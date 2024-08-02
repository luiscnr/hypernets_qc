import os
from datetime import datetime as dt
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
from plot_multiple import PlotMultiple
from flag_manager import Flags
from plot_options import PlotOptions
from plot_spectra import PlotSpectra

# import sys
# code_home = os.path.dirname(os.path.dirname(__init__.__file__))
# sys.path.append(code_home)

class HYPERNETS_DAY_FILE_LAND():

    def __init__(self, file_nc, path_images):
        self.VALID = True
        if not os.path.exists(file_nc):
            self.VALID = False
        if file_nc is None:
            self.VALID = False
        self.file_nc = file_nc
        self.path_images = path_images
        if self.path_images is None:
            self.path_images = os.path.dirname(self.file_nc)
        self.path_images_date = None
        self.sequences = []
        self.format_img = '.png'
        if self.VALID:
            self.sequences = self.get_sequences()
            self.isequence = 0
            self.valid_sequences = None

        self.ref_wl = 800
        self.ref_wl_idx = -1
        self.rho = r'ρ$_w$'

        self.flag_builder = None
        self.sequences_no_data = None

        self.flags_rgb = ['target_rad1', 'target_rad2', 'target_rad3', 'sky_irr1', 'sky_irr2', 'sun']

        self.context = {
            'bad_pointing_threshold_zenith': 3,
            'bad_pointing_threshold_azimuth': 3,
            'plot_polar_min': 0,
            'plot_polar_max': 0.8,
            'legendfontsize': 8,
            'fontsize': 14,
            'ylim_irradiance': None,
            'ylim_reflectance': None,
            'linestyles': None
        }

    def get_sequences(self):
        from netCDF4 import Dataset
        import numpy as np
        sequences = []
        dataset = Dataset(self.file_nc)
        seq_var = dataset.variables['sequence_ref']
        for idx in range(len(seq_var)):
            val = seq_var[idx]
            if np.ma.is_masked(val):
                sequences.append(None)
                continue
            time = dt.utcfromtimestamp(float(val))
            sequences.append(time.strftime('%Y%m%dT%H%M'))
        dataset.close()
        return sequences

    def set_path_images_date(self, site, date_here):
        folder_date = os.path.join(self.path_images, site, date_here.strftime('%Y'), date_here.strftime('%m'),
                                   date_here.strftime('%d'))
        if os.path.exists(folder_date):
            self.path_images_date = folder_date
        else:
            self.path_images_date = None

    def save_report_image(self, site, delete_images, overwrite):
        print(f'[INFO] Sequence {self.isequence}: SEQ{self.sequences[self.isequence]}')
        if self.sequences[self.isequence] is None:
            return
        file_out = os.path.join(os.path.dirname(self.file_nc),
                                f'{site}_{self.sequences[self.isequence]}_Report{self.format_img}')
        if os.path.exists(file_out) and not overwrite:
            return
        dir_img = os.path.join(os.path.dirname(self.file_nc), 'IMG')
        if not os.path.exists(dir_img):
            os.mkdir(dir_img)

        file_img = os.path.join(dir_img, f'CameraImages_{self.isequence}_all{self.format_img}')
        if not os.path.isfile(file_img) or (os.path.isfile(file_img) and overwrite):
            self.save_img_files(True)

        file_middle = os.path.join(dir_img, f'Angles_{self.isequence}_all{self.format_img}')
        if not os.path.isfile(file_middle) or (os.path.isfile(file_middle) and overwrite):
            self.save_middle_panel(file_middle)

        file_spectra = os.path.join(dir_img, f'Spectra_{self.isequence}_all{self.format_img}')
        if not os.path.exists(file_spectra) or (os.path.isfile(file_spectra) and overwrite):
            self.save_spectra_files(file_spectra)

        pm = PlotMultiple()
        nrow = 3
        ncol = 1
        pm.start_multiple_plot_advanced(nrow, ncol, 10, 18, 0.1, 0.1, True)
        pm.get_axes(0, 0).set_title(self.get_title(site), fontsize=20)
        pm.plot_image(file_img, 0, 0)
        flag_list = self.get_flags_sequence()
        if len(flag_list) == 0:
            str_flag_list = 'NO FLAGGED'
        else:
            str_list = ','.join(flag_list)
            str_flag_list = f'FLAGS:{str_list}'
        pm.get_axes(1, 0).set_title(str_flag_list, fontsize=12)
        pm.plot_image(file_middle, 1, 0)
        # info_l2 = self.get_info_l2()
        info_l2 = ''
        pm.get_axes(2, 0).set_title(info_l2, fontsize=12)
        pm.plot_image(file_spectra, 2, 0)

        pm.save_fig(file_out)
        pm.close_plot()

        # if delete_images:
        #     for name in os.listdir(dir_img):
        #         file_here = os.path.join(dir_img, name)
        #         os.remove(file_here)
        #     os.rmdir(dir_img)

    def get_flags_sequence(self):
        from netCDF4 import Dataset
        import numpy as np
        dataset = Dataset(self.file_nc)
        flag_value_series = np.uint64(dataset.variables['l2_quality_flag'][self.isequence])
        all_flag_values = [np.uint64(x) for x in dataset.variables['l2_quality_flag'].flag_masks.split(',')]
        all_flag_meaninigs = dataset.variables['l2_quality_flag'].flag_meanings
        cflags = Flags(all_flag_values, all_flag_meaninigs)
        final_list = {}
        for flag_value in flag_value_series:
            list, mask = cflags.Decode(flag_value)
            for l in list:
                if l not in final_list.keys():
                    final_list[l] = 1
                else:
                    final_list[l] = final_list[l] + 1
        dataset.close()
        list = []
        for fl in final_list:
            list.append(f'{fl}(x{final_list[fl]})')

        return list

    def get_title(self, site):
        date_time_here = dt.strptime(self.sequences[self.isequence], '%Y%m%dT%H%M')
        date_str = date_time_here.strftime('%Y-%m-%d')
        time_str = date_time_here.strftime('%H:%M')
        title = f'{site} {self.sequences[self.isequence]} - {date_str} {time_str} - {self.isequence + 1}/{len(self.sequences)}'
        return title

    def save_spectra_files(self, file_spectra):
        dir_img = os.path.join(os.path.dirname(self.file_nc), 'IMG')
        if not os.path.exists(dir_img):
            os.mkdir(dir_img)
        from netCDF4 import Dataset
        dataset = Dataset(self.file_nc)
        vza = np.round(dataset.variables['l2_viewing_zenith_angle'][0, :])
        vaa = np.round(dataset.variables['l2_viewing_azimuth_angle'][0, :])

        index_1 = np.logical_and(vza == 30, vaa == 113).nonzero()[0].tolist()
        index_2 = np.logical_and(vza == 0, vaa == 98).nonzero()[0].tolist()
        index_3 = np.logical_and(vza == 30, vaa == 293).nonzero()[0].tolist()
        dataset.close()

        file_1 = os.path.join(dir_img, f'Reflectance_{self.isequence}_1{self.format_img}')
        self.plot_reflectance_impl(file_1, index_1, None, 'DEFAULT')
        file_2 = os.path.join(dir_img, f'Reflectance_{self.isequence}_2{self.format_img}')
        self.plot_reflectance_impl(file_2, index_2, None, 'DEFAULT')
        file_3 = os.path.join(dir_img, f'Reflectance_{self.isequence}_3{self.format_img}')
        self.plot_reflectance_impl(file_3, index_3, None, 'DEFAULT')

        pm = PlotMultiple()
        nrow = 1
        ncol = 2
        file_12 = os.path.join(dir_img, f'Reflectance_{self.isequence}_12{self.format_img}')
        pm.start_multiple_plot_advanced(nrow, ncol, 10, 5, 0.35, 0, True)
        pm.plot_image(file_1, 0, 0)
        pm.plot_image(file_2, 0, 1)
        pm.save_fig(file_12)
        pm.close_plot()

        pmf = PlotMultiple()
        nrow = 2
        ncol = 1
        pmf.start_multiple_plot_advanced(nrow, ncol, 10, 8, 0.1, 0, True)
        pmf.plot_image(file_12, 0, 0)
        pmf.plot_image(file_3, 1, 0)
        pmf.save_fig(file_spectra)
        pmf.close_plot()

    def save_middle_panel(self, file_middle):
        dir_img = os.path.join(os.path.dirname(self.file_nc), 'IMG')
        if not os.path.exists(dir_img):
            os.mkdir(dir_img)

        ##irradiance
        file_irradiance = os.path.join(dir_img, f'Irradiance_{self.isequence}{self.format_img}')
        self.plot_irradiance_impl(file_irradiance, None)

        ##polar
        file_polar = os.path.join(dir_img, f'Polar_{self.isequence}{self.format_img}')
        pm = PlotMultiple()
        nrow = 1
        ncol = 2
        pm.start_multiple_plot_polar(nrow, ncol, 10, 5, 0.35, 0, True)
        self.ref_wl = 680
        ax_here = pm.get_axes(0, 0)
        self.plot_polar_files_impl(None, ax_here)
        self.ref_wl = 800
        ax_here = pm.get_axes(0, 1)
        self.plot_polar_files_impl(None, ax_here)
        pm.save_fig(file_polar)
        pm.close_plot()

        pm_middle = PlotMultiple()
        pm_middle.start_multiple_plot_advanced(2, 1, 10, 8, 0.1, 0, True)
        pm_middle.plot_image(file_irradiance, 0, 0)
        pm_middle.plot_image(file_polar, 1, 0)
        pm_middle.save_fig(file_middle)
        pm_middle.close_plot()

    def plot_irradiance_impl(self, file_out, series_id):
        from netCDF4 import Dataset
        dataset = Dataset(self.file_nc)
        irradiance = dataset.variables['l1_irradiance'][self.isequence]
        xdata = dataset.variables['wavelength'][:]

        ##if series_id is None, plot first and last spectra
        if series_id is None:
            series_id = [0, irradiance.shape[0] - 1]
        labels = None
        if len(series_id) > 1:
            labels = []
            times = dataset.variables['l2_acquisition_time'][self.isequence]
            for id in series_id:
                time_here = dt.utcfromtimestamp(times[id])
                labels.append(time_here.strftime('%Y-%m-%d %H:%M:%S'))
        fig1, ax1 = plt.subplots(figsize=(10, 5))
        if labels is None:
            ydata = irradiance[series_id[0], :]
            ax1.plot(xdata, ydata, alpha=0.3)
        else:
            ydata = irradiance[series_id, :]
            for idx in range(len(series_id)):
                ax1.plot(xdata, ydata[idx], label=labels[idx], alpha=0.3)
            ax1.legend(fontsize=self.context['legendfontsize'])
        ax1.set_xlabel("Wavelength (nm)", fontsize=self.context['fontsize'])
        ax1.set_ylabel(r"Irradiance ($mW\ nm^{-1}\ m^{-2}$)", fontsize=self.context['fontsize'])
        ylim = self.context['ylim_irradiance']
        if ylim is not None:
            ax1.set_ylim(ylim)
        else:
            ymax = np.percentile(ydata, 95) * 1.2
            if np.isfinite(ymax):
                ax1.set_ylim([0, ymax])
        fig1.savefig(file_out, dpi=300, bbox_inches="tight")
        plt.close(fig1)
        dataset.close()

    def plot_reflectance_impl(self, file_out, series_id, type_labels, title):
        from netCDF4 import Dataset
        dataset = Dataset(self.file_nc)
        reflectance = dataset.variables['l2_reflectance'][self.isequence]
        xdata = dataset.variables['wavelength'][:]
        labels = None
        if len(series_id) > 1 and type_labels is not None:
            # to implement
            labels = []
        if title is not None and title == 'DEFAULT':
            if labels is None:  ##len(series_id)==1, only one series
                iseries = series_id[0]
                vza = np.round(dataset.variables['l2_viewing_zenith_angle'][self.isequence, iseries])
                vaa = np.round(dataset.variables['l2_viewing_azimuth_angle'][self.isequence, iseries])
                time = dataset.variables['l2_acquisition_time'][self.isequence, iseries]
                time_here = dt.utcfromtimestamp(time)
                title = f'{time_here.strftime("%Y-%m-%d %H:%M:%S")}(vza={vza:.0f},vaa={vaa:.0f})'
        dataset.close()
        # fig1, ax1 = plt.subplots(figsize=(10, 5))
        fig1, ax1 = plt.subplots()
        if labels is None:
            ydata = reflectance[series_id[0], :]
            ax1.plot(xdata, ydata, alpha=0.3, linewidth=0, marker='o', color='orangered', markeredgewidth=0)
        else:
            linestyles = self.context['linestyles']
            if linestyles is not None and len(linestyles) != len(labels):
                linestyles = None
            ydata = reflectance[series_id, :]
            for idx in range(len(labels)):
                if linestyles is None:
                    ax1.plot(xdata, ydata[idx], label=labels[idx], alpha=0.3)
                else:
                    ax1.plot(xdata, ydata[idx], label=labels[idx], ls=linestyles[idx], alpha=0.3)
            if len(labels) > 10:
                ax1.legend(
                    bbox_to_anchor=(1.04, 1),
                    loc="upper left",
                    fontsize=self.context['legendfontsize'],
                    ncol=2,
                )
            else:
                ax1.legend(loc="lower right", fontsize=self.context['legendfontsize'])
        ax1.set_xlabel("Wavelength (nm)", fontsize=self.context['fontsize'])
        ax1.set_ylabel(r"Reflectance", fontsize=self.context['fontsize'])
        if title is not None:
            ax1.set_title(title, fontsize=self.context['fontsize'] + 2)
        ylim = self.context['ylim_reflectance']
        if ylim is not None:
            ax1.set_ylim(ylim)
        else:
            ymax = np.percentile(ydata, 95) * 1.2
            if np.isfinite(ymax):
                ax1.set_ylim([0, ymax])
        if file_out is not None:
            fig1.savefig(file_out, bbox_inches="tight", dpi=300)
        plt.close(fig1)

    def plot_polar_files_impl(self, file_out, ax):
        from netCDF4 import Dataset
        dataset = Dataset(self.file_nc)
        saa = np.mean(dataset.variables['l2_solar_azimuth_angle'][self.isequence, :] % 360)
        sza = np.mean(dataset.variables['l2_solar_zenith_angle'][self.isequence, :])
        vaa = dataset.variables['l2_viewing_azimuth_angle'][self.isequence, :] % 360

        vza = dataset.variables['l2_viewing_zenith_angle'][self.isequence, :]
        refl = dataset.variables['l2_reflectance'][self.isequence, :, :]

        vaa_grid = np.arange(8, 368, 15)
        vza_grid = np.array([0, 5, 10, 20, 30, 40, 50, 60])
        raa_grid = vaa_grid - saa

        self.ref_wl_idx = np.argmin(np.abs(self.ref_wl - dataset.variables['wavelength'][:]))
        print(f'[INFO][Polar plot] Wavelength indces for {self.ref_wl}: {self.ref_wl_idx}')

        vaa_mesh, vza_mesh = np.meshgrid(np.radians(vaa_grid), vza_grid)

        refl_2d = np.zeros((len(vaa_grid), len(vza_grid)))
        for i in range(len(vaa_grid)):
            for ii in range(len(vza_grid)):
                id_series = np.where(
                    (np.abs(vaa - vaa_grid[i]) < self.context["bad_pointing_threshold_azimuth"]) &
                    (np.abs(vza - vza_grid[ii]) < self.context["bad_pointing_threshold_zenith"])
                )[0]
                if len(id_series) > 0:
                    if len(id_series) == 1:
                        refl_2d[i, ii] = np.abs(refl[id_series, self.ref_wl_idx])
                    else:
                        refl_2d[i, ii] = np.mean(np.abs(refl[id_series, self.ref_wl_idx]))

        refl_2d[refl_2d == 0] = np.nan
        dataset.close()

        ##plotting
        fig = None
        if ax is None:
            fig = plt.figure()
            ax = plt.subplot(1, 1, 1, projection="polar")

        ax.set_theta_direction(-1)
        ax.set_theta_offset(np.pi / 2.0)
        im = ax.pcolormesh(
            vaa_mesh,
            vza_mesh,
            refl_2d.T,
            shading="auto",
            cmap=plt.get_cmap("jet"),
            vmin=self.context["plot_polar_min"],
            vmax=self.context["plot_polar_max"],
        )

        ax.plot(np.radians(saa), sza, color="k", ls="none", marker="o")
        if fig is not None:
            cbar = fig.colorbar(im)
            cbar.set_label("reflectance at %s nm" % self.ref_wl, rotation=270, labelpad=15)
        else:
            cbar = plt.colorbar(im, pad=0.15)
            cbar.set_label("reflectance at %s nm" % self.ref_wl, rotation=270, labelpad=15)

        if fig is not None:
            if file_out is not None:
                fig.savefig(file_out, dpi=300)
            plt.close(fig)

    def save_img_files(self, multiple_plot):
        # flags = ['sky_irr_1', 'sky_rad_1', 'water_rad', 'sky_rad_2', 'sky_irr_2', 'sun']
        dir_img = os.path.join(os.path.dirname(self.file_nc), 'IMG')
        if not os.path.exists(dir_img):
            os.mkdir(dir_img)
        file_out_base = os.path.join(dir_img, f'CameraImages_{self.isequence}')
        if multiple_plot:
            pm = PlotMultiple()
            nrow = 2
            ncol = 3
            pm.start_multiple_plot_advanced(nrow, ncol, 10, 7.0, 0.15, 0.15, True)
        index_row = 0
        index_col = 0
        for flag in self.flags_rgb:
            file_img, title = self.get_img_file(flag)

            if index_col == ncol:
                index_col = 0
                index_row = index_row + 1

            if multiple_plot:
                if file_img is not None:
                    pm.plot_image_hypernets(file_img, index_row, index_col, title)
                else:

                    pm.plot_blank_with_title(index_row, index_col, title)

            index_col = index_col + 1
        if multiple_plot:
            pm.save_fig(f'{file_out_base}_all{self.format_img}')
            pm.close_plot()

    def get_img_file(self, flag):
        from netCDF4 import Dataset
        import numpy as np
        if flag.startswith('pictures_'):
            name_var = flag
        else:
            name_var = f'pictures_{flag}'
        dataset = Dataset(self.file_nc)
        if name_var not in dataset.variables:
            return None

        var = dataset.variables[name_var]
        prefix = var.prefix
        suffix = var.suffix
        seq_here = self.sequences[self.isequence]

        val = var[self.isequence]
        file_img = None
        if not np.ma.is_masked(val):
            time = dt.utcfromtimestamp(float(val))
            time_str = time.strftime('%Y%m%dT%H%M')
            name_file_img = f'{prefix}_{seq_here}_{time_str}_{suffix}'
            if self.path_images_date is not None:
                file_img = os.path.join(self.path_images_date, name_file_img)
                if not os.path.exists(file_img):
                    file_img = None

        title = f'{flag}(VZA={var.oza};VAA={var.oaa})'
        dataset.close()
        return file_img, title

    def plot_from_options_impl(self, options_figure):
        if not options_figure['apply']:
            return None
        if options_figure['type'] == 'spectraplot' and options_figure['type_rrs'] == 'user_defined':
            print(f'[INFO] Plotting spectra plot...')
            self.plot_spectra_plot_from_options(options_figure)
        if options_figure['type'] == 'timeseries':
            print(f'[INFO] Plotting time series plot...')
            self.plot_time_series_from_options(options_figure)
        if options_figure['type'] == 'sequence':
            print(f'[INFO] Plotting sequence plot...')
            daily_sequence_summary = self.plot_sequence_plot_from_options(options_figure)
            return daily_sequence_summary
        if options_figure['type'] == 'flagplot' and options_figure['type_flagplot'] == 'comparison':
            print(f'[INFO] Plotting flag plot...')
            self.plot_flag_plot_comparison(options_figure)
        if options_figure['type'] == 'angleplot':
            print(f'[INFO] Plotting angle plot...')
            self.plot_angle_plot_from_options(options_figure)

        return None

    def plot_sequence_plot_from_options(self, options_figure):
        from netCDF4 import Dataset
        import numpy as np
        import pytz
        import pandas as pd
        import seaborn as sns
        from matplotlib import pyplot as plt
        dataset = Dataset(self.file_nc)
        time_array = dataset.variables['l2_acquisition_time'][:]
        time_array = np.ma.masked_values(time_array, 0)  ##solving a problem find 0n 30/05/2023, it shouln't be happen

        start_time_real = dt.utcfromtimestamp(np.min(time_array[:]))
        end_time_real = dt.utcfromtimestamp(np.max(time_array[:]))
        if start_time_real.strftime('%Y%m%d') != end_time_real.strftime('%Y%m%d'):
            dataset.close()
            print('[ERROR] Plot is only created for a single day')
            return

        use_default_flag = True
        if options_figure['flagBy'] is not None:
            options_figure['flagType'] = 'flag'
            options_figure = self.check_gs_options_impl(options_figure, 'flagBy', 'flagType', 'flagValues')
            # noptions = len(options_figure[options_figure['flagBy']]['flag_values'])
            if options_figure['flagBy'] in options_figure.keys():
                use_default_flag = False
                flag_meanings = options_figure[options_figure['flagBy']]['flag_meanings']
                flag_values = options_figure[options_figure['flagBy']]['flag_values']
                flag_array = options_figure[options_figure['flagBy']]['flag_array']

        if use_default_flag:
            dataset.close()
            print(f'[ERROR] Defaults options for land are not set')
            return
        time_fix_axis = self.get_fix_axis_time(options_figure, start_time_real, end_time_real)
        ntime = len(time_fix_axis)
        time_fix_min_max = np.zeros((ntime, 2))
        seconds_ref = self.get_time_interval_seconds(options_figure['frequency'], 'minutes')
        time_fix_axis_ts = np.array([x.replace(tzinfo=pytz.utc).timestamp() for x in time_fix_axis]).astype(np.float64)
        time_fix_min_max[:, 0] = time_fix_axis_ts - seconds_ref
        time_fix_min_max[:, 1] = time_fix_axis_ts + seconds_ref
        time_ticks = []
        for itime in range(ntime):
            htick = time_fix_axis[itime].strftime('%H')
            mtick = time_fix_axis[itime].strftime('%M')
            time_ticks.append(f'{htick}:{mtick}')
        nseries = flag_array.shape[1]

        vza = np.round(dataset.variables['l2_viewing_zenith_angle'][0])
        vaa = np.round(dataset.variables['l2_viewing_azimuth_angle'][0])
        series_ticks = []
        for idx in range(nseries):
            series_ticks.append(f'{vza[idx]:.0f},{vaa[idx]:.0f}')



        daily_summary_sequences = {
            'NTotal': 0,
            'NAvailable': 0,
            'start_time': time_fix_axis[0].strftime('%Y-%m-%d %H:%M'),
            'end_time': time_fix_axis[-1].strftime('%Y-%m-%d %H:%M'),
            'expected_sequences': ntime
        }
        legend_values = flag_meanings
        if options_figure['legendValues'] is not None:
            legend_values = options_figure['legendValues']
        for f in legend_values:
            daily_summary_sequences[f] = {
                'nsequences': 0,
                'nsequences_completed': 0
            }

        data = pd.DataFrame(index=series_ticks, columns=time_ticks).astype(np.float64)
        # data[:] = 0
        for itime in range(ntime):
            htick = time_fix_axis[itime].strftime('%H')
            mtick = time_fix_axis[itime].strftime('%M')
            time_tick = f'{htick}:{mtick}'

            time_valid = np.logical_and(time_array >= time_fix_min_max[itime, 0],
                                        time_array < time_fix_min_max[itime, 1])

            if np.count_nonzero(time_valid) == nseries:
                daily_summary_sequences['NAvailable'] = daily_summary_sequences['NAvailable'] + 1
                flag_array_series = flag_array[time_valid]

                for iseries in range(nseries):
                    try:
                        fvalue = flag_array_series[iseries]
                        index_flag = flag_values.index(int(fvalue))
                        # flag_meaning = legend_values[index_flag]
                        # daily_summary_sequences[flag_meaning]['nseries'] = daily_summary_sequences[flag_meaning]['nseries']+1
                        data.loc[series_ticks[iseries]].at[time_tick] = index_flag + 1
                    except:
                        pass

        for flag_meaning, flag_value in zip(flag_meanings, flag_values):
            for time_tick in time_ticks:
                y_array = data.loc[:, time_tick]
                nvalues = np.count_nonzero(y_array == flag_value + 1)
                if nvalues > 0:
                    daily_summary_sequences[flag_meaning]['nsequences'] = daily_summary_sequences[flag_meaning][
                                                                              'nsequences'] + 1
                if nvalues == nseries:
                    daily_summary_sequences[flag_meaning]['nsequences_completed'] = \
                        daily_summary_sequences[flag_meaning]['nsequences_completed'] + 1

        if self.sequences_no_data is not None and len(self.sequences_no_data) > 0:
            daily_summary_sequences['NTotal'] = daily_summary_sequences['NAvailable'] + len(self.sequences_no_data)
            for seq in self.sequences_no_data:
                time_stamp_seq = dt.strptime(seq[3:], '%Y%m%dT%H%M').replace(tzinfo=pytz.utc).timestamp()
                iwhere = np.where(
                    np.logical_and(time_stamp_seq >= time_fix_min_max[:, 0], time_stamp_seq < time_fix_min_max[:, 1]))
                if len(iwhere[0]) > 0:
                    index = iwhere[0][0]
                    date_seq = dt.utcfromtimestamp(time_fix_axis_ts[index])
                    data.loc[:, date_seq.strftime('%H:%M')] = 0
        else:
            daily_summary_sequences['NTotal'] = daily_summary_sequences['NAvailable']

        print('---------------->',len(data.index))
        print(data.columns)

        plt.Figure()
        if options_figure['color'] is not None:
            colors = ['gainsboro'] + options_figure['color']
        else:
            colors = ['gainsboro', 'red', 'cyan', 'green', 'magenta']
        vmax = len(colors) - 1
        ax = sns.heatmap(data, vmin=0, vmax=vmax, cmap=colors, linewidths=0.5, linecolor='gray')
        plt.yticks(rotation='horizontal')
        colorbar = ax.collections[0].colorbar
        ticks = np.arange(0.5, vmax, 0.75).tolist()
        if options_figure['legendTicks'] is not None:
            ticks = options_figure['legendTicks']
        colorbar.set_ticks(ticks)

        legend_values = ['NO-L2'] + legend_values
        colorbar.set_ticklabels(legend_values, rotation=90)
        colorbar.ax.tick_params(size=0)
        plt.xlabel('Time')
        plt.ylabel('Series(vza,vaa)')

        if options_figure['title'] is not None:
            title = options_figure['title']
            title = title.replace('$DATE$', start_time_real.strftime('%Y-%m-%d'))
            plt.title(title, fontsize=options_figure['fontsizetitle'])

        plt.gcf().tight_layout()
        file_out = options_figure['file_out']
        if file_out is not None:
            plt.savefig(file_out, dpi=300)

        plt.close()

        return daily_summary_sequences

    def get_fix_axis_time(self, options_figure, start_time_real, end_time_real):
        start_date = options_figure['start_date']
        end_date = options_figure['end_date']
        if start_date is None:
            start_date = start_time_real.strftime('%Y-%m-%d')
        if end_date is None:
            end_date = end_time_real.strftime('%Y-%m-%d')
        start_time = options_figure['start_time']
        end_time = options_figure['end_time']
        if start_time is None:
            start_time = start_time_real.strftime('%H:%M')
        if end_time is None:
            end_time = end_time_real.strftime('%H:%M')
        start_date_time = dt.strptime(f'{start_date}T{start_time}', '%Y-%m-%dT%H:%M')
        end_date_time = dt.strptime(f'{end_date}T{end_time}', '%Y-%m-%dT%H:%M')
        print(f'[INFO] Fix time axis from {start_date_time} to {end_date_time}')
        frequency = float(options_figure['frequency'])
        print(f'[INFO] ->Frequency: {frequency} {options_figure["frequency_units"]}')
        if frequency == -999.0:
            return None

        poptions = PlotOptions(None, None)
        return poptions.get_fix_time_axis(frequency, options_figure['frequency_units'], start_date_time, end_date_time)

    def plot_flag_plot_comparison(self, options_figure):


        dataset = Dataset(self.file_nc)
        qf_array = dataset.variables['l2_quality_flag'][:]
        nsequences = qf_array.shape[0]
        nseries = qf_array.shape[1]

        qf_flags_meanings = dataset.variables['l2_quality_flag'].flag_meanings
        qf_flags_list = qf_flags_meanings.split(' ')
        qf_flags_values = [np.uint64(x) for x in dataset.variables['l2_quality_flag'].flag_masks.split(',')]
        ylabel = ['ALLVALID','VALID'] + qf_flags_list
        yarray = np.zeros((len(ylabel)))
        fcheck = Flags(qf_flags_values, qf_flags_meanings)
        nflag_total = 0

        for index, value in enumerate(qf_flags_values):
            mask = fcheck.Mask(qf_array, [qf_flags_list[index]])
            mask = np.ma.filled(mask, 0)

            nflag_bysequence = np.count_nonzero(mask > 0, axis=1)
            nflag_bysequence[nflag_bysequence > 0] = 1
            nflag = np.sum(nflag_bysequence)
            if nflag > 0:
                yarray[index + 1] = nflag
                nflag_total = nflag_total + nflag

        valid = qf_array == 0
        nvalid_bysequence = np.count_nonzero(valid, axis=1)
        nvalid_bysequence_all = nvalid_bysequence.copy()
        nvalid_bysequence_all[nvalid_bysequence == nseries] = 1
        nvalid_bysequence_all[nvalid_bysequence < nseries] = 0
        yarray[0] = np.sum(nvalid_bysequence_all)
        nvalid_bysequence[nvalid_bysequence>0]=1
        yarray[1] = np.sum(nvalid_bysequence)


        xarray = np.arange(len(yarray))
        plt.barh(xarray, yarray)
        xarray = np.append(xarray, xarray[-1] + 1)
        ylabel.append('')
        plt.yticks(xarray - 0.5, ylabel, fontsize=8, verticalalignment='bottom')
        if nsequences <= 5:
            xticks = [0, 1, 2, 3, 4, 5]
        elif nsequences >= 5 and nsequences < 10:
            xticks = [0, 2, 4, 6, 8, 10]
        else:
            xmax = 5 * (np.floor(nsequences / 5) + 1)
            xticks = np.arange(0, xmax + 1, 5).tolist()
        plt.xticks(xticks, [f'{x:.0f}' for x in xticks])
        plt.xlabel(f'Number of sequences')
        plt.grid()
        plt.gcf().tight_layout()
        file_out = options_figure['file_out']
        if file_out is not None:
            plt.savefig(file_out, dpi=300)

        dataset.close()

    def plot_time_series_from_options(self, options_figure):
        from netCDF4 import Dataset
        import numpy as np
        import pytz

        dataset = Dataset(self.file_nc)
        time_var = options_figure['time_var']
        avg_vars = options_figure['avg_var']

        if time_var is None or avg_vars is None:
            print(f'[ERROR] time_var and avg_var should be defined in the configuration file')
            return
        if time_var not in dataset.variables:
            print(f'[ERROR] {time_var} is not defined in {self.mrfile.file_path}')
            return
        else:
            nsequences = dataset.variables[time_var].shape[0]
            nseries = dataset.variables[time_var].shape[1]

        is_spectral = []
        for avg_var in avg_vars:
            if avg_var not in dataset.variables:
                print(f'[ERROR] {avg_var} is not defined in {self.mrfile.file_path}')
                return
            else:
                sh = dataset.variables[avg_var].shape
                if len(sh)==3:
                    is_spectral.append(True)
                elif len(sh) == 2:
                    is_spectral.append(False)


        options_figure['wlref_index'] = -1
        if options_figure['wlref'] != -999.0 and is_spectral.count(True) > 0:
            wavelength = dataset.variables['wavelength'][:]
            index_ref = np.argmin(np.abs(options_figure['wlref'] - wavelength))
            if abs(wavelength[index_ref] - options_figure['wlref']) <= 2:
                options_figure['wlref_index'] = index_ref

        options_figure = self.check_gs_options_impl(options_figure, 'groupBy', 'groupType', 'groupValues')
        time_array = dataset.variables[time_var][:]
        time_array = np.ma.filled(time_array.astype(np.float64), -999.0)



        ngroups, groupValues, groupArray, str_legend = self.check_group_and_legend(options_figure)

        if options_figure['legend_values'] is not None and len(options_figure['legend_values']) == ngroups:
            str_legend = options_figure['legend_values']

        start_time_real = dt.utcfromtimestamp(np.min(time_array[time_array != -999.0]))
        end_time_real = dt.utcfromtimestamp(np.max(time_array[time_array != -999.0]))



        # time_fix_axis = None
        time_fix_min_max = None
        if options_figure['type_time_axis'] == 'fix':
            time_fix_axis = self.get_fix_axis_time(options_figure, start_time_real, end_time_real)
            ntime = len(time_fix_axis)
            time_fix_min_max = np.zeros((ntime, 2))
            seconds_ref = self.get_time_interval_seconds(options_figure['frequency'], options_figure['frequency_units'])
            time_fix_axis_ts = np.array([x.replace(tzinfo=pytz.utc).timestamp() for x in time_fix_axis]).astype(
                np.float64)
            time_fix_min_max[:, 0] = time_fix_axis_ts - seconds_ref
            time_fix_min_max[:, 1] = time_fix_axis_ts + seconds_ref
            xarray = np.arange(ntime)

            xvalues = xarray.tolist()
            xticks = [x.strftime('%H:%M') for x in time_fix_axis]
            if options_figure['xticks_range'] > 1:
                if options_figure['xticks_labels_range'] > 1:
                    indices_new = list(range(0, ntime, options_figure['xticks_labels_range']))
                    for idx in range(ntime):
                        if idx not in indices_new:
                            xticks[idx] = ''
                indices_new = list(range(0, ntime, options_figure['xticks_range']))
                xvalues = [xvalues[index] for index in indices_new]
                xticks = [xticks[index] for index in indices_new]



        ##PLOTTING
        pspectra = PlotSpectra()
        pspectra.xdata = xarray

        if ngroups == 1:  ##no groups,we could use multiple variables
            pass
        elif ngroups > 1:  ##multiple groups, only one variable, optional multiple series
            avg_var = avg_vars[0]
            avg_array = dataset.variables[avg_var][:]
            avg_array = np.ma.filled(avg_array.astype(np.float64), -999.0)
            if is_spectral[0]:
                index_ref = options_figure['wlref_index']
                avg_array = avg_array[:,:,index_ref]
                avg_array = np.squeeze(avg_array)

            print(f'[INFO] NSequences:  {nsequences}. NSeries: {nseries}. Variable: {avg_var}. Array shape: {avg_array.shape}')
            markers = options_figure['marker']
            nseries_ref = len(options_figure['iseriesref'])

            handles = []
            str_legend_valid = []
            for icolor, gvalue in enumerate(groupValues):
                color = options_figure['color'][icolor]
                legend_added = False
                for idx in range(nsequences):
                    for idx_series,iseries in enumerate(options_figure['iseriesref']):

                        xdata, ydata = self.get_data_nospectral_fix_time(time_fix_min_max[idx, 0],
                                                                     time_fix_min_max[idx, 1], time_array,
                                                                     avg_array, idx,iseries, groupArray, gvalue)
                        if xdata is not None and ydata is not None:
                            marker_here = markers[0]
                            if len(markers)==nseries_ref:
                                marker_here = markers[idx_series]
                            h = pspectra.plot_single_marker(xdata, ydata, marker_here, 6, color, None, 0)
                            if not legend_added:
                                handles.append(h[0])
                                str_legend_valid.append(icolor)
                                legend_added = True

            pspectra.set_grid()
            if options_figure['xlabel'] is not None:
                pspectra.set_xaxis_title(options_figure['xlabel'])
            if options_figure['ylabel'] is not None:
                pspectra.set_yaxis_title(options_figure['ylabel'])

            ##legend
            if len(str_legend) > 0 and len(str_legend_valid) > 0:
                str_legend = [str_legend[int(x)] for x in str_legend_valid]
                pspectra.legend_options['loc'] = 'lower center'
                pspectra.legend_options['bbox_to_anchor'] = (0.5, -0.25)
                pspectra.legend_options['ncols'] = len(str_legend)
                if len(handles) == 0:
                    pspectra.set_legend(str_legend)
                elif len(handles) == len(str_legend):
                    pspectra.set_legend_h(handles, str_legend)

            ##y-range
            pspectra.set_y_range(options_figure['y_min'], options_figure['y_max'])

            pspectra.set_xticks(xvalues, xticks, 0, 10)
            pspectra.set_tigth_layout()

            file_out = options_figure['file_out']
            if file_out is not None:
                pspectra.save_fig(file_out)
            pspectra.close_plot()
        dataset.close()


    def plot_angle_plot_from_options(self, options_figure):
        from netCDF4 import Dataset
        import numpy as np
        options_figure = self.check_gs_options_impl(options_figure, 'groupBy', 'groupType', 'groupValues')
        dataset = Dataset(self.file_nc)
        if options_figure['angle_var'] not in dataset.variables:
            print(f'[ERROR] Angle var: {options_figure["angle_var"]} is not available')
            dataset.close()
            return
        if options_figure['avg_var'] not in dataset.variables:
            print(f'[ERROR] Angle var: {options_figure["avg_var"]} is not available')
            dataset.close()
            return

        angle_array = dataset.variables[options_figure['angle_var']][:]
        angle_array = np.ma.filled(angle_array.astype(np.float64), -999.0)


        avg_array = dataset.variables[options_figure['avg_var']][:]
        avg_array = np.ma.filled(avg_array, -999.0)
        is_spectral = False
        if len(avg_array.shape) == 3:is_spectral = True

        if is_spectral:
            index_ref = -1
            if options_figure['wlref'] != -999.0:
                wavelength = dataset.variables['wavelength'][:]
                index_ref = np.argmin(np.abs(options_figure['wlref'] - wavelength))
                if abs(wavelength[index_ref] - options_figure['wlref']) > 2:
                    index_ref = -1
            if index_ref == -1:
                print(f'[ERROR] index_ref for wavelength {options_figure["wlref"]} is not valid')
                return
            avg_array = np.squeeze(avg_array[:, :, index_ref])


        ngroups, groupValues, groupArray, str_legend = self.check_group_and_legend(options_figure)

        ##PLOTTING
        fig, ax_here = plt.subplots(figsize=(3, 3), subplot_kw={'projection': 'polar'})
        ax_here = plt.gca()
        ax_here.set_rscale(options_figure['scale'])
        ax_here.set_rlim(options_figure['rlim'])
        ax_here.set_rticks(options_figure['rticks'])
        ax_here.set_theta_zero_location(options_figure['theta_zero_location'])
        ax_here.set_theta_direction(options_figure['theta_direction'])

        valid_array = np.logical_and(angle_array != -999.0, avg_array != -999.0)

        for idx_series, iseries in enumerate(options_figure['iseriesref']):
            marker_series = None
            if len(options_figure['marker_series'])==len(options_figure['iseriesref']) and len(options_figure['iseriesref'])>1:
                marker_series = options_figure['marker_series'][idx_series]
            angle_array_series = angle_array[:,iseries]
            avg_array_series = avg_array[:,iseries]
            valid_array_series = valid_array[:,iseries]
            groupArray_series = groupArray[:, iseries]
            angle_array_series = angle_array_series[valid_array_series]
            avg_array_series = avg_array_series[valid_array_series]
            if ngroups > 1 and groupArray is not None:
                groupArray_series = groupArray_series[valid_array_series]
            if options_figure['min_data'] is not None:
                avg_array_series[avg_array_series < options_figure['min_data']] = options_figure['min_data'] + 1e-6

            if ngroups > 1 and groupArray is not None:
                colors = options_figure['color']
                if len(colors) != ngroups:
                    import MDB_reader.MDBPlotDefaults as pdefaults
                    colors = pdefaults.get_color_list(ngroups)
                point_size = options_figure['point_size']
                point_marker = options_figure['point_marker']
                if len(point_size) != ngroups:
                    point_size = [point_size[0]] * ngroups
                if len(point_marker) != ngroups:
                    point_marker = [point_marker[0]] * ngroups
                for ig, g in enumerate(groupValues):
                    valid_g = groupArray_series == g
                    if np.count_nonzero(valid_g) > 0:
                        angle_array_rad = (angle_array_series[valid_g] * np.pi) / 180
                        avg_array_here = avg_array_series[valid_g]
                        marker_here = marker_series if marker_series is not None else point_marker[ig]
                        ax_here.scatter(angle_array_rad, avg_array_here, marker=marker_here, s=point_size[ig],color=colors[ig])
            else:
                angle_array_rad = (angle_array_series * np.pi) / 180
                marker_here = marker_series if marker_series is not None else options_figure['point_marker'][0]
                ax_here.scatter(angle_array_rad, avg_array_series, marker=marker_here,s=options_figure['point_size'][0], color=options_figure['color'][0])

            if options_figure['title'] is not None:
                ax_here.set_title(options_figure['title'])

            ax_here.set_thetamin(options_figure['theta_min'])
            ax_here.set_thetamax(options_figure['theta_max'])
            ax_here.set_xticks((np.pi * np.array(options_figure['xticks']) / 180))
            if options_figure['rlabel_position'] is not None:
                if options_figure['rlabel_position'] == 'AUTO':
                    rpos = self.get_rlabel_position(angle_array)
                else:
                    rpos = float(options_figure['rlabel_position'])
                ax_here.set_rlabel_position(rpos)

            ax_here.tick_params(axis='x', labelsize=options_figure['label_size'])
            ax_here.tick_params(axis='y', labelsize=options_figure['label_size'])
            plt.tight_layout()
            file_out = options_figure['file_out']
            if file_out is not None:
                plt.savefig(file_out, dpi=300)

        dataset.close()

    def plot_spectra_plot_from_options(self, options_figure):
        from netCDF4 import Dataset
        import numpy as np
        plot_spectra = True
        if options_figure['plot_spectra'][0].lower() == 'none':
            plot_spectra = False
        plot_stats = options_figure['plot_stats']
        if not plot_spectra and not plot_stats:
            print('[WARNING] Please active plot_spectra or plot_stats')
            return

        ##start potential virtual flags
        options_figure = self.check_gs_options_impl(options_figure, 'groupBy', 'groupType', 'groupValues')

        ##DATA SELECTION
        dataset = Dataset(self.file_nc)
        wavelength = np.array(dataset.variables[options_figure['wl_variable']])
        y_variable = options_figure['y_variable']
        var_y = dataset.variables[y_variable]
        # ndim = len(var_y.shape)
        if '_FillValue' in var_y.ncattrs():
            fill_value = var_y._FillValue
        else:
            from netCDF4 import default_fillvals
            fill_value = default_fillvals[var_y.dtype.str[1:]]

        iseries = options_figure['iseriesref']
        spectra = np.squeeze(np.array(dataset.variables[y_variable][:, iseries, :]))

        ##Checking spectra with fill values
        spectra_check = np.where(spectra == fill_value, 1, 0)
        spectra_check = np.sum(spectra_check, axis=1)

        spectra = spectra[spectra_check == 0, :]

        dataset.close()
        ##DATA SELECTION

        ##CHECKING GROUP AND LEGEND
        ngroup, groupValues, groupArray, str_legend = self.check_group_and_legend(options_figure)
        groupArray = groupArray[:,iseries]
        handles = []
        str_legend_valid = []
        if ngroup > 1:
            groupArray = groupArray[spectra_check == 0]


        ##PLOTTING
        pspectra = PlotSpectra()
        pspectra.xdata = wavelength

        # single spectra plotting
        if plot_spectra:
            line_color = options_figure['color']
            marker = options_figure['marker']
            marker_size = options_figure['markersize']
            line_type = options_figure['linestyle']
            line_size = options_figure['linewidth']
            if ngroup > 1:
                if len(line_color) != ngroup:  ##assing default colors
                    line_color = [line_color[0]] * ngroup
                    for idx in range(ngroup):
                        line_color[idx] = self.get_color_default(idx, 0, ngroup)
                if len(marker) != ngroup:
                    marker = [marker[0]] * ngroup
                if len(marker_size) != ngroup:
                    marker_size = [marker_size[0]] * ngroup
                if len(line_type) != ngroup:
                    line_type = [line_type[0]] * ngroup
                if len(line_size) != ngroup:
                    line_size = [line_size[0]] * ngroup

                for idx in range(ngroup):
                    val = groupValues[idx]
                    hline = pspectra.plot_single_line(spectra[groupArray == val, :].transpose(), line_color[idx],
                                                      line_type[idx], line_size[idx], marker[idx], marker_size[idx])
                    if len(hline) > 0:
                        handles.append(hline[0])
                        str_legend_valid.append(idx)
            else:
                pspectra.plot_single_line(spectra.tranpose(), line_color[0], line_type[0], line_size[0], marker[0],
                                          marker_size[0])

        if plot_stats:
            line_color = options_figure['color']
            if ngroup > 1:
                if len(line_color) != ngroup:  ##assing default colors
                    line_color = [line_color[0]] * ngroup
                    for idx in range(ngroup):
                        line_color[idx] = self.get_color_default(idx, 0, ngroup)
                for idx in range(ngroup):
                    val = groupValues[idx]
                    spectra_here = spectra[groupArray == val, :]
                    if len(spectra_here) > 0:
                        stats_here = self.get_spectra_stats(spectra_here)
                        pspectra.stats_style['fill']['color'] = line_color[idx]
                        pspectra.stats_style['central']['color'] = line_color[idx]
                        hline = pspectra.plot_stats(stats_here, None, None)
                        handles.append(hline[0])
                        str_legend_valid.append(idx)
            else:
                stats = self.get_spectra_stats(spectra)
                pspectra.plot_stats(stats, None, None)

        ##legend
        if len(str_legend) > 0 and len(str_legend_valid) > 0:
            str_legend = [str_legend[int(x)] for x in str_legend_valid]
            pspectra.legend_options = pspectra.legend_options_bottom
            pspectra.legend_options['ncols'] = len(str_legend)

            if len(handles) == 0:
                pspectra.set_legend(str_legend)
            elif len(handles) == len(str_legend):
                pspectra.set_legend_h(handles, str_legend)

        ##y-range
        pspectra.set_y_range(options_figure['y_min'], options_figure['y_max'])

        ##label and titles
        if options_figure['xlabel'] is not None:
            pspectra.set_xaxis_title(options_figure['xlabel'])
        if options_figure['ylabel'] is not None:
            pspectra.set_yaxis_title(options_figure['ylabel'])
        if options_figure['title'] is not None:
            pspectra.set_title(options_figure['title'])
        pspectra.set_grid_horizontal()
        # saveing to file
        if options_figure['file_out'] is not None:
            file_out = options_figure['file_out']
            pspectra.save_fig(file_out)

        pspectra.close_plot()


    def save_report_summary_image(self, site, date_here, dir_img_summary, daily_sequences_summary):

        file_out = os.path.join(os.path.dirname(self.file_nc),
                                f'{site}_{date_here.strftime("%Y%m%d")}_DailySummary{self.format_img}')
        print(f'[INFO] Output file: {file_out}')


        ##TIME SERIES
        # start_multiple_plot_advanced(self, nrow, ncol, xfigsize, yfigsize, wspace, hspace, frameon)
        file_out_ts = os.path.join(dir_img_summary, f'TS{self.format_img}')
        pmts = PlotMultiple()
        pmts.start_multiple_plot_advanced(2, 1, 3, 3.5, 0, 0, True)
        pmts.plot_image(os.path.join(dir_img_summary, 'time_series_r680.tif'), 0, 0)
        pmts.plot_image(os.path.join(dir_img_summary, 'time_series_r800.tif'), 1, 0)
        pmts.save_fig(file_out_ts)
        pmts.close_plot()

        ##TOP PANEL
        file_out_top = os.path.join(dir_img_summary, f'TOP{self.format_img}')
        pmtop = PlotMultiple()
        pmtop.start_multiple_plot_advanced(1, 3, 10, 3.5, 0, 0, True)
        pmtop.plot_image(os.path.join(dir_img_summary, 'sequence_info.tif'), 0, 0)
        pmtop.plot_image(os.path.join(dir_img_summary, 'flag_plot.tif'), 0, 1)
        pmtop.plot_image(file_out_ts, 0, 2)
        pmtop.set_text(-1250, 50, f'DAILY SUMMARY REPORT - {date_here.strftime("%Y-%m-%d")}')
        # daily_sequences_summary = self.get_sequence_info()
        line = f'Total sequences: {daily_sequences_summary["NTotal"]}/{daily_sequences_summary["expected_sequences"]}. Processed to L2: {daily_sequences_summary["NAvailable"]}.'
        skip = ['NTotal','NAvailable','start_time','end_time','expected_sequences']
        for key in daily_sequences_summary.keys():
            if key in skip: continue
            nsequences = daily_sequences_summary[key]['nsequences']
            nsequences_completed = daily_sequences_summary[key]['nsequences_completed']
            if nsequences==0: continue
            #if daily_sequences_summary[key]==0: continue

            line = f'{line} {key}: {nsequences}({nsequences_completed}). '

        pmtop.set_text_size(-1750, 850, line, 8)
        pmtop.save_fig(file_out_top)
        pmtop.close_plot()

        ##MIDDLE(ANGLE-PANEL)
        file_out_middle = os.path.join(dir_img_summary, f'MIDDLE{self.format_img}')
        pmmiddle = PlotMultiple()
        pmmiddle.start_multiple_plot_advanced(2, 2, 10, 6, 0.35, 0, True)
        pmmiddle.plot_image(os.path.join(dir_img_summary, 'angle_680_sza.tif'), 0, 0)
        pmmiddle.plot_image(os.path.join(dir_img_summary, 'angle_680_saa.tif'), 0, 1)
        pmmiddle.plot_image(os.path.join(dir_img_summary, 'angle_800_sza.tif'), 1, 0)
        pmmiddle.plot_image(os.path.join(dir_img_summary, 'angle_800_saa.tif'), 1, 1)

        pmmiddle.save_fig(file_out_middle)
        pmmiddle.close_plot()


        ##IRRADIANCE PANEL
        file_out_irr = os.path.join(dir_img_summary, f'IRR{self.format_img}')
        pirr = PlotMultiple()
        pirr.start_multiple_plot_advanced(1, 2, 10, 3, 0.35, 0, True)
        pirr.plot_image(os.path.join(dir_img_summary, 'irradiance_1.tif'), 0, 0)
        pirr.plot_image(os.path.join(dir_img_summary, 'irradiance_2.tif'), 0, 1)
        pirr.save_fig(file_out_irr)
        pirr.close_plot()

        #Reflectance panel
        file_out_ref = os.path.join(dir_img_summary, f'REF{self.format_img}')
        pref = PlotMultiple()
        pref.start_multiple_plot_advanced(1, 3, 10, 3, 0.35, 0, True)
        pref.plot_image(os.path.join(dir_img_summary, 'reflectance_1.tif'), 0, 0)
        pref.plot_image(os.path.join(dir_img_summary, 'reflectance_2.tif'), 0, 1)
        pref.plot_image(os.path.join(dir_img_summary, 'reflectance_3.tif'), 0, 2)
        pref.save_fig(file_out_ref)
        pref.close_plot()

        ##DOWN(SPECTRA-PANEL)
        file_out_down = os.path.join(dir_img_summary, f'DOWN{self.format_img}')
        pdown = PlotMultiple()
        pdown.start_multiple_plot_advanced(2, 1, 10, 6, 0, 0, True)
        pdown.plot_image(file_out_irr, 0, 0)
        pdown.plot_image(file_out_ref, 1, 0)
        pdown.save_fig(file_out_down)
        pdown.close_plot()

        ##final plot
        pm = PlotMultiple()
        pm.start_multiple_plot_advanced(3, 1, 10, 18, 0, 0, True)
        pm.plot_image(file_out_top, 0, 0)
        pm.plot_image(file_out_middle, 1, 0)
        pm.plot_image(file_out_down, 2, 0)
        pm.save_fig(file_out)
        pm.close_plot()

        ##remove blank space
        import matplotlib.image as mpimg
        import numpy as np
        image = mpimg.imread(file_out)
        image_new = np.concatenate([image[0:901, :, :], image[1200:image.shape[0], :, :]])
        height, width, nbands = image_new.shape
        figsize = width / float(300), height / float(300)
        plt.close()
        fig = plt.figure(figsize=figsize)
        ax = fig.add_axes([0, 0, 1, 1])
        ax.imshow(image_new, interpolation='nearest')
        ax.axis(False)
        fig.savefig(file_out, dpi=300, bbox_inches="tight")
        plt.close(fig)

        # ##delete intermeddiate figures
        # for name in os.listdir(dir_img_summary):
        #     os.remove(os.path.join(dir_img_summary, name))
        #
        # os.rmdir(dir_img_summary)

    def get_color_default(self, value, min, max):
        nvalues = (max - min) + 1
        colors_default = ['Blue', 'Red', 'Green', 'm', 'Cyan', 'Orange', 'Yellow']
        if nvalues < 6:
            index = value - min
            return colors_default[index]
        import matplotlib as mpl
        cm = mpl.colormaps['jet']
        return cm((value - min) / (max - min))

    def get_rlabel_position(self, angle):
        median_op_angle = np.median(angle) + 180
        if median_op_angle > 360:
            median_op_angle = median_op_angle - 360
        angle_labels = np.arange(0, 361, 45)
        angle_label = angle_labels[int(np.argmin(np.abs(angle_labels - median_op_angle)))]
        if angle_label == 270:
            angle_label = 225
        if angle_label == 90:
            angle_label = 135
        return angle_label

    def get_data_nospectral_fix_time(self,time_min, time_max, time_array, var_array, output_value,iseries,groupArray, groupValue):


        if groupArray is not None:
            valid_time = np.logical_and(np.logical_and(time_array >= time_min, time_array < time_max),
                                        groupArray == groupValue)
        else:
            valid_time = np.logical_and(time_array >= time_min, time_array < time_max)


        var_array_series =var_array[:,iseries]
        valid_time_series = valid_time[:,iseries]

        nvalid = np.count_nonzero(valid_time_series)

        if nvalid == 0:
            return None, None

        xdata = np.array([output_value] * nvalid).astype(np.float64)
        ydata = var_array_series[valid_time_series]

        return xdata, ydata


    def check_group_and_legend(self, options_figure):

        ngroup = 1
        groupValues = None
        groupArray = None
        if 'groupValues' in options_figure.keys():
            groupValues = options_figure['groupValues']
        if groupValues is not None:
            ngroup = len(groupValues)



        str_legend = []

        if ngroup > 1 and options_figure['legend']:
            str_legend = self.get_str_legend(options_figure)



        if ngroup > 1:

            groupArray, all_flag_values, all_flag_meanings = self.get_gs_array(options_figure,
                                                                               options_figure['groupBy'],
                                                                               options_figure['groupType'])
            if 'color' in options_figure.keys():
                if len(options_figure['color']) != ngroup:
                    import MDB_reader.MDBPlotDefaults as default
                    options_figure['color'] = default.get_color_list(ngroup)

        return ngroup, groupValues, groupArray, str_legend

    def get_gs_array(self, options_figure, by, type):
        from netCDF4 import Dataset
        import numpy as np
        dataset = Dataset(self.file_nc)
        if type == 'float':
            array_flag = np.array(dataset.variables[by])
            all_flag_values = None
            all_flag_meanings = None
        else:
            if by in dataset.variables:
                array_flag = np.array(dataset.variables[by][:])
                all_flag_values = dataset.variables[by].flag_values
                all_flag_meanings = dataset.variables[by].flag_meanings.split(' ')
            else:  ##previously built as a virtual flag in method check_gs_options_impl, line 904
                array_flag = options_figure[by]['flag_array']
                all_flag_values = options_figure[by]['flag_values']
                all_flag_meanings = options_figure[by]['flag_meanings']
        dataset.close()

        return array_flag, all_flag_values, all_flag_meanings


    def get_str_legend(self, options):
        if options['legend_values'] is not None:
            return options['legend_values']
        str_legend = []
        groupValues = options['groupValues']
        groupValues = np.unique(np.array(groupValues)).tolist()
        if groupValues is not None:
            ngroup = len(groupValues)
            if ngroup > 1:
                if options['groupType'] == 'float' or options['groupType'] == 'wavelength':
                    for g in groupValues:
                        str_legend.append(f'{g:.2f}')
                if options['groupType'] == 'flag':
                    flag_name = options['groupBy']
                    str_legend = self.get_flag_list(groupValues, options[flag_name]['flag_values'],
                                                    options[flag_name]['flag_meanings'])
                    if 'FUB' in str_legend:
                        index = str_legend.index('FUB')
                        if index >= 0:
                            str_legend[index] = 'S3 FUB-CSIRO'
                    if 'STANDARD' in str_legend:
                        index = str_legend.index('STANDARD')
                        if index >= 0:
                            str_legend[index] = 'WFR'
                    if 'POLYMER' in str_legend:
                        index = str_legend.index('POLYMER')
                        if index >= 0:
                            str_legend[index] = 'CMEMS-OLCI'
                    if 'CCIALL' in str_legend:
                        index = str_legend.index('CCIALL')
                        if index >= 0:
                            str_legend[index] = 'OC-CCI v.6 (complete time series)'
                    if 'CCI' in str_legend:
                        index = str_legend.index('CCI')
                        if index >= 0:
                            str_legend[index] = 'OC-CCI v.6 (OLCI period)'

        return str_legend

    def get_flag_list(self, values, allValues, allFlags):
        import numpy as np
        flag_list = []
        for val in values:
            if val == -1:
                flag_list.append('GLOBAL')
            indext = np.where(np.array(allValues) == val)
            index = indext[0]
            if len(index) == 1:
                indexf = index[0]
                flag_list.append(allFlags[indexf])
        return flag_list


    def get_time_interval_seconds(self, frequency, units):
        interval = frequency / 2
        if units == 'minutes':
            interval = interval * 60
        elif units == 'hours':
            interval = interval * 60 * 60
        elif units == 'days':
            interval = interval * 60 * 60 * 24
        elif units == 'months':
            interval = interval * 60 * 60 * 24 * 30
        elif units == 'years':
            interval = interval * 60 * 60 * 24 * 365
        return interval

    def check_gs_options_impl(self, options_figure, by, type, values):
        from netCDF4 import Dataset
        import numpy as np
        var_group_name = options_figure[by]
        if var_group_name is None:
            return options_figure

        dataset = Dataset(self.file_nc)
        if options_figure[type] == 'flag':
            if var_group_name in dataset.variables:
                flag_values = dataset.variables[var_group_name].flag_values
                flag_meanings_list = dataset.variables[var_group_name].flag_meanings.split(' ')
                flag_meanings = [x.strip() for x in flag_meanings_list]
                options_figure[var_group_name] = {
                    'flag_values': flag_values,
                    'flag_meanings': flag_meanings
                }
            else:  ##virtual flag
                virtual_flags_options = self.flag_builder.get_virtual_flags_options()
                if virtual_flags_options[var_group_name]['type'] == 'ranges':
                    array, flag_meanings, flag_values = self.flag_builder.create_flag_array_ranges_v2(
                        virtual_flags_options[var_group_name])
                elif virtual_flags_options[var_group_name]['type'] == 'flag':
                    array, dims, flag_meanings, flag_values = self.flag_builder.create_flag_array_flag(
                        virtual_flags_options[var_group_name])

                options_figure[var_group_name] = {
                    'flag_values': flag_values,
                    'flag_meanings': flag_meanings,
                    'flag_array': array
                }

            if options_figure[values] is None:
                options_figure[values] = flag_values
            else:
                flag_list_config = options_figure[values]
                flag_values_config = []
                for flag_config in flag_list_config:
                    if flag_config.strip() == 'GLOBAL':
                        flag_values_config.append(-1)
                        continue
                    try:
                        iflag = flag_meanings.index(flag_config.strip())
                        flag_values_config.append(flag_values[iflag])
                    except:
                        print(f'[WARNING] Flag {flag_config.strip()} is not in the list')
                        return None
                options_figure[values] = flag_values_config

        if options_figure[type] == 'float':
            if var_group_name not in dataset.variables:
                return None
            all_group_values = np.unique(np.array(dataset.variables[var_group_name]))
            if options_figure[values] is None:
                options_figure[values] = list(all_group_values)
            else:
                group_values_given = options_figure[values]
                group_values = []
                for val in group_values_given:
                    imin = np.argmin(np.abs(val - all_group_values))
                    if abs(val - all_group_values[imin]) < 0.1:
                        group_values.append((all_group_values[imin]))
                    else:
                        print(f'[WARNING] Value {val} is not in the variable {var_group_name}')
                        return None
                options_figure[values] = group_values
        dataset.close()
        return options_figure

    def save_report_image_only_pictures(self, site, delete_images, overwrite, seq, files_img):
        print(f'[INFO] Sequence {seq} (No Level-2 data available)')
        seq_time_str = seq[3:]
        seq_time = dt.strptime(seq_time_str, '%Y%m%dT%H%M')
        file_out = os.path.join(os.path.dirname(self.file_nc), f'{site}_{seq_time_str}_Report{self.format_img}')
        if os.path.exists(file_out) and not overwrite:
            return
        names_img = self.flags_rgb
        names_img_files = {}
        for ref in files_img:
            name_img = files_img[ref]['name_img']
            if name_img in names_img:
                names_img_files[name_img] = files_img[ref]['file_img']

        dir_img = os.path.join(os.path.dirname(self.file_nc), 'IMG')
        if not os.path.exists(dir_img):
            os.mkdir(dir_img)

        # file_out_base = os.path.join(dir_img, f'CameraImages_{self.isequence}')

        pm = PlotMultiple()
        nrow = 2
        ncol = 3
        pm.start_multiple_plot_advanced(nrow, ncol, 10, 7.0, 0.02, 0.15, True)
        index_row = 0
        index_col = 0
        for name_img in names_img_files:
            file_img = names_img_files[name_img]
            title = name_img
            if index_col == ncol:
                index_col = 0
                index_row = index_row + 1

            if file_img is not None:
                pm.plot_image_hypernets(file_img, index_row, index_col, title)
            else:
                pm.plot_blank_with_title(index_row, index_col, title)

            index_col = index_col + 1

        date_str = seq_time.strftime('%Y-%m-%d')
        time_str = seq_time.strftime('%H:%M')
        title = f'{site} {seq_time_str} - {date_str} {time_str} - No L2 Data'
        pm.fig.suptitle(title)
        line = f'ANOMALY: ?'
        pm.fig.text(0.20, 0.05, line)
        # pm.get_axes(0,0).set_title(title)
        pm.save_fig(file_out)
        pm.close_plot()