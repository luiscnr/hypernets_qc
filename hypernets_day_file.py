import os, configparser, pytz
from datetime import datetime as dt
from datetime import timedelta
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import matplotlib as mpl
import numpy as np
import pandas as pd
import seaborn as sns
import statistics as st
from netCDF4 import Dataset
from netCDF4 import default_fillvals
from math import floor, ceil
import plot_defaults as pdefaults
from plot_multiple import PlotMultiple
from plot_options import PlotOptions
from plot_spectra import PlotSpectra
from flag_manager import Flags


class HYPERNETS_DAY_FILE:

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

    def get_sequences(self):

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

    def get_metadata_files(self, input_path):
        metadata_files = []
        dataset = Dataset(self.file_nc)
        seq_var = dataset.variables['sequence_ref']
        metadata_files_keys = {}
        for name in os.listdir(input_path):
            if name.endswith('_metadata.txt'):
                time_metadata = dt.strptime(name[3:name.find('_metadata.txt')], '%Y%m%dT%H%M%S')
                metadata_files_keys[time_metadata.strftime('%Y%m%dT%H%M')] = os.path.join(input_path, name)

        for idx in range(len(seq_var)):
            val = seq_var[idx]
            if np.ma.is_masked(val):
                metadata_files.append(None)
                continue
            time = dt.utcfromtimestamp(float(val))
            time_str = time.strftime('%Y%m%dT%H%M')
            if time_str in metadata_files_keys:
                metadata_files.append(metadata_files_keys[time_str])
            else:
                metadata_files.append(None)

        dataset.close()
        return metadata_files

    def get_angles_from_metadata_files(self, metadata_files):
        if metadata_files is None:
            return None

        offset_tilt = 56
        nseries = len(metadata_files)
        array_angles = np.zeros((nseries, 3)).astype(np.float64)  # paa, vaa, vza
        array_angles[:] = np.nan
        for idx, file in enumerate(metadata_files):
            if file is not None:
                options = configparser.ConfigParser()
                options.read(file)
                if options.has_option('01_008_0270_2_0040', 'pt_ref'):
                    val = options['01_008_0270_2_0040']['pt_ref']
                    paa_ref = float(val.split(';')[0].strip())
                    vza_ref = float(val.split(';')[1].strip())
                    array_angles[idx, 0] = self.normalizeddeg(paa_ref, 0, 360)
                    array_angles[idx, 1] = self.normalizeddeg(array_angles[idx, 0] + 180, 0, 360)
                    array_angles[idx, 2] = self.normalizeddeg(vza_ref + float(offset_tilt), 0, 360)

        return array_angles

    def normalizeddeg(self, num, lower=0.0, upper=360.0, b=False):
        res = num
        if not b:
            if lower >= upper:
                raise ValueError(
                    "Invalid lower and upper limits: (%s, %s)" % (lower, upper)
                )

            res = num
            if num > upper or num == lower:
                num = lower + abs(num + upper) % (abs(lower) + abs(upper))
            if num < lower or num == upper:
                num = upper - abs(num - lower) % (abs(lower) + abs(upper))

            res = lower if res == upper else num
        else:
            total_length = abs(lower) + abs(upper)
            if num < -total_length:
                num += ceil(num / (-2 * total_length)) * 2 * total_length
            if num > total_length:
                num -= floor(num / (2 * total_length)) * 2 * total_length
            if num > upper:
                num = total_length - num
            if num < lower:
                num = -total_length - num

            res = num * 1.0  # Make all numbers float, to be consistent

        return res

    ##methods for multiple dates. start_date and end_date are dt objects, start_time and end_time with format %H:%M
    ##these parameters are retrieved from options_figures using check_dates_times
    def get_sequences_range(self, start_date, end_date, start_time, end_time):
        if start_date is None and end_date is None and start_time is None and end_time is None:
            return None, None
        if self.sequences is None:
            self.get_sequences()
        sequences_here = []
        sequences_indices = None
        for iseq in range(len(self.sequences)):
            seq = self.sequences[iseq]
            if seq is None:
                continue
            time_seq = dt.strptime(seq, '%Y%m%dT%H%M')
            if start_time is not None or end_time is not None:
                date_seq_str = time_seq.strftime('%Y-%m-%d')
            if start_date is not None:
                if time_seq < start_date:
                    continue
            if end_date is not None:
                end_date = end_date.replace(hour=23, minute=59, second=59)
                if time_seq > end_date:
                    continue
            if start_time is not None:
                start_time_ref = dt.strptime(f'{date_seq_str}T{start_time}', '%Y-%m-%dT%H:%M')
                if time_seq < start_time_ref:
                    continue
            if end_time is not None:
                end_time_ref = dt.strptime(f'{date_seq_str}T{end_time}', '%Y-%m-%dT%H:%M')
                if time_seq > end_time_ref:
                    continue
            sequences_here.append(seq)
            if sequences_indices is None:
                sequences_indices = [iseq]
            else:
                sequences_indices.append(iseq)

        return sequences_here, sequences_indices

    ##method only if file is for a specific data
    def get_sequences_interval(self, start_time, end_time):
        if self.sequences is None:
            self.get_sequences()
        sequences_here = []
        sequences_indices = None
        if len(self.sequences) == 0:
            return sequences_here
        for iseq in range(len(self.sequences)):
            seq = self.sequences[iseq]
            if seq is None:
                continue
            time_seq = dt.strptime(seq, '%Y%m%dT%H%M')
            if start_time <= time_seq <= end_time:
                sequences_here.append(seq)
                if sequences_indices is None:
                    sequences_indices = [iseq]
                else:
                    sequences_indices.append(iseq)

        return sequences_here, sequences_indices

    def get_report_files_interval(self, sequences_here, site, start_time, end_time):
        if sequences_here is None:
            sequences_here, range_here = self.get_sequences_interval(start_time, end_time)
        report_files = []
        if len(sequences_here) == 0:
            return report_files
        for seq in sequences_here:
            file_out = os.path.join(os.path.dirname(self.file_nc), f'{site}_{seq}_Report{self.format_img}')
            if os.path.exists(file_out):
                report_files.append(file_out)
        return report_files

    def get_ref_wl_idx(self):

        dataset = Dataset(self.file_nc)
        wavelength = np.array(dataset.variables['wavelength'][:])
        self.ref_wl_idx = np.argmin(np.abs(wavelength - self.ref_wl))
        dataset.close()

    def reflectance_ref(self, apply_nosc):

        if self.ref_wl_idx == -1:
            self.get_ref_wl_idx()
        variable = 'l2_reflectance'
        if apply_nosc:
            variable = 'l2_reflectance_nosc'
        dataset = Dataset(self.file_nc)
        reflectance = np.array(dataset.variables[variable][:, self.ref_wl_idx])
        dataset.close()
        return reflectance

    def reflectance_ref_l1(self, apply_nosc):

        if self.ref_wl_idx == -1:
            self.get_ref_wl_idx()
        variable = 'l1_reflectance'
        if apply_nosc:
            variable = 'l1_reflectance_nosc'
        dataset = Dataset(self.file_nc)
        reflectance = np.array(dataset.variables[variable][self.isequence, self.ref_wl_idx, :])
        dataset.close()
        return reflectance

    def angle_rad_level1(self, flag):

        if flag == 'sza':
            variable = 'l1_solar_zenith_angle'
        elif flag == 'saa':
            variable = 'l1_solar_azimuth_angle'
        elif flag == 'paa':
            variable = 'l1_pointing_azimuth_angle'
        dataset = Dataset(self.file_nc)
        angle = np.array(dataset.variables[variable][self.isequence, :])
        dataset.close()
        angle = (angle * np.pi) / 180
        return angle

    def angle_rad(self, flag):

        if flag == 'sza':
            variable = 'l2_solar_zenith_angle'
        elif flag == 'saa':
            variable = 'l2_solar_azimuth_angle'
        elif flag == 'paa':
            variable = 'l2_pointing_azimuth_angle'
        dataset = Dataset(self.file_nc)
        angle = np.array(dataset.variables[variable][:])
        dataset.close()
        median_op_angle = np.median(angle) + 180
        if median_op_angle > 360:
            median_op_angle = median_op_angle - 360
        angle_labels = np.arange(0, 361, 45)
        angle_label = angle_labels[int(np.argmin(np.abs(angle_labels - median_op_angle)))]
        if angle_label == 270:
            angle_label = 225
        if angle_label == 90:
            angle_label = 135
        angle = (angle * np.pi) / 180
        return angle, angle_label

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

    def get_valid_flags(self):

        dataset = Dataset(self.file_nc)
        valid_sequences = dataset.variables['l2_quality_flag'][:].astype(np.float64)
        valid_sequences = np.ma.filled(valid_sequences, -999.0)
        valid_sequences[valid_sequences > 0] = 1
        epsilon_array = dataset.variables['l2_epsilon'][:].astype(np.float64)
        epsilon_array = np.ma.filled(epsilon_array, -999.0)
        valid_sequences[np.logical_and(valid_sequences == 0, epsilon_array < (-0.05))] = 2
        valid_sequences[np.logical_and(valid_sequences == 0, epsilon_array >= 0.05)] = 3
        self.valid_sequences = valid_sequences
        dataset.close()

    def get_flags_sequence(self):

        dataset = Dataset(self.file_nc)
        flag_value = np.uint64(dataset.variables['l2_quality_flag'][self.isequence])
        all_flag_values = [np.uint64(x) for x in dataset.variables['l2_quality_flag'].flag_masks.split(',')]
        all_flag_meaninigs = dataset.variables['l2_quality_flag'].flag_meanings
        cflags = Flags(all_flag_values, all_flag_meaninigs)
        list, mask = cflags.Decode(flag_value)
        epsilon_here = dataset.variables['l2_epsilon'][self.isequence]
        if not np.ma.is_masked(epsilon_here):
            if epsilon_here < (-0.05): list.append('ENEG')
            if epsilon_here >= 0.05: list.append('EHIGH')
        dataset.close()
        return list

    def get_info_l2(self):

        dataset = Dataset(self.file_nc)
        epsilon = dataset.variables['l2_epsilon'][self.isequence]
        rho = dataset.variables['l2_rhof'][self.isequence]
        raa = dataset.variables['l1_rhof_raa'][self.isequence, 0]
        sza = dataset.variables['l1_rhof_sza'][self.isequence, 0]
        vza = dataset.variables['l1_rhof_vza'][self.isequence, 0]
        ws = dataset.variables['l1_rhof_wind'][self.isequence, 0]
        dataset.close()
        str = f'epsilon={epsilon:.4f};rho={rho:.4f}(raa={raa:.1f};sza={sza:.1f};vza={vza:.1f};ws={ws:.2f})'
        return str

    def get_title(self, site):
        date_time_here = dt.strptime(self.sequences[self.isequence], '%Y%m%dT%H%M')
        date_str = date_time_here.strftime('%Y-%m-%d')
        time_str = date_time_here.strftime('%H:%M')
        title = f'{site} {self.sequences[self.isequence]} - {date_str} {time_str} - {self.isequence + 1}/{len(self.sequences)}'
        return title

    def get_array_variable(self, name_var):

        dataset = Dataset(self.file_nc)
        array = np.array(dataset.variables[name_var][:])
        dataset.close()
        return array

    def creating_copy_with_new_angles(self, array_angles):

        file_nc_old = self.file_nc[:-3] + '_old.nc'
        if os.path.exists(file_nc_old):
            print('[WARNING] Angle correction has already been done. Skipping...')
            return
        os.rename(self.file_nc, file_nc_old)
        file_nc_new = self.file_nc
        # file_nc_old = self.file_nc
        # file_nc_new = self.file_nc[:-3] + '_new.nc'

        input_dataset = Dataset(file_nc_old)
        ncout = Dataset(file_nc_new, 'w', format='NETCDF4')

        # copy attribtues
        ncout.setncatts(input_dataset.__dict__)

        # copy dimensions
        for name, dimension in input_dataset.dimensions.items():
            ncout.createDimension(
                name, (len(dimension) if not dimension.isunlimited() else None))

        nseries = array_angles.shape[0]
        for name, variable in input_dataset.variables.items():
            fill_value = None
            if '_FillValue' in list(variable.ncattrs()):
                fill_value = variable._FillValue

            ncout.createVariable(name, variable.datatype, variable.dimensions, fill_value=fill_value, zlib=True,
                                 complevel=6)

            # copy variable attributes all at once via dictionary
            ncout[name].setncatts(input_dataset[name].__dict__)

            if name == 'l2_pointing_azimuth_angle':
                ncout[name][:] = array_angles[:, 0]
            elif name == 'l2_viewing_azimuth_angle':
                ncout[name][:] = array_angles[:, 1]
            elif name == 'l2_viewing_zenith_angle':
                ncout[name][:] = array_angles[:, 2]
            elif name == 'l1_pointing_azimuth_angle':
                nscan = input_dataset[name].shape[1]
                array_new = array_angles[:, 0].reshape(nseries, 1).repeat(nscan, axis=1)
                ncout[name][:] = array_new[:]
            elif name == 'l1_viewing_azimuth_angle':
                nscan = input_dataset[name].shape[1]
                array_new = array_angles[:, 1].reshape(nseries, 1).repeat(nscan, axis=1)
                ncout[name][:] = array_new[:]
            elif name == 'l1_viewing_zenith_angle':
                nscan = input_dataset[name].shape[1]
                array_new = array_angles[:, 2].reshape(nseries, 1).repeat(nscan, axis=1)
                ncout[name][:] = array_new[:]
            else:
                ncout[name][:] = input_dataset[name][:]
        input_dataset.close()

        ncout.close()

    # flag: name_variable, sky_irr_1, sky_irr_2, sky_rad_1, sky_rad_1, water_rad, sun
    def get_img_file(self, flag):

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

        #check if exist seq folder (RBINS)
        time_seq_here = dt.strptime(seq_here,'%Y%m%dT%H%M')
        path_images_here = self.path_images_date
        for name in os.listdir(self.path_images_date):
            if name.startswith('SEQ') and os.path.isdir(os.path.join(self.path_images_date,name,'image')):
                time_seq_folder = dt.strptime(name[3:],'%Y%m%dT%H%M%S')
                if abs(time_seq_folder-time_seq_here).total_seconds()<150.0:
                    path_images_here = os.path.join(self.path_images_date,name,'image')


        val = var[self.isequence]
        file_img = None
        if not np.ma.is_masked(val):
            time = dt.utcfromtimestamp(float(val))
            time_str = time.strftime('%Y%m%dT%H%M')
            name_file_img = f'{prefix}_{seq_here}_{time_str}_{suffix}'
            if self.path_images_date is not None:
                file_img = os.path.join(path_images_here, name_file_img)
                if not os.path.exists(file_img):
                    file_img = None

        title = f'{flag}(az={var.oaa};zn={var.oza})'
        dataset.close()
        return file_img, title

    def set_path_images_date(self, site, date_here):
        folder_date = os.path.join(self.path_images, site, date_here.strftime('%Y'), date_here.strftime('%m'),
                                   date_here.strftime('%d'))
        if os.path.exists(folder_date):
            self.path_images_date = folder_date
        else:
            self.path_images_date = None

    def get_water_images(self, site, date_here, time_min, time_max, interval_minutes):

        if time_min is None:
            time_min = '0600'
        if time_max is None:
            time_max = '1700'
        if interval_minutes is None:
            interval_minutes = 20
        self.set_path_images_date(site, date_here)

        date_here_str = date_here.strftime('%Y%m%d')
        date_here_min = dt.strptime(f'{date_here_str}T{time_min}', '%Y%m%dT%H%M')
        date_here_max = dt.strptime(f'{date_here_str}T{time_max}', '%Y%m%dT%H%M')
        date_here = date_here_min
        water_images = {}
        while date_here <= date_here_max:
            date_here_str = date_here.strftime('%Y%m%dT%H%M')
            water_images[date_here_str] = {
                'file_img': None,
                'title': None
            }
            date_here = date_here + timedelta(minutes=interval_minutes)

        dataset = Dataset(self.file_nc)
        img_var = dataset.variables['pictures_water_rad']
        for idx in range(len(img_var)):

            val = img_var[idx]
            if np.ma.is_masked(val):
                continue
            self.isequence = idx
            file_img, title = self.get_img_file('water_rad')
            seq = self.sequences[self.isequence]
            if seq in water_images.keys():
                water_images[seq]['file_img'] = file_img
                paa = float(dataset.variables['l2_pointing_azimuth_angle'][idx])
                sza = float(dataset.variables['l2_solar_zenith_angle'][idx])
                saa = float(dataset.variables['l2_solar_azimuth_angle'][idx])
                title = f'{seq} (sza={sza:.1f};paa={paa:.1f};saa={saa:.1f})'
                title = f'{seq} (paa={paa:.1f})'
                water_images[seq]['title'] = title
        dataset.close()

        return water_images

    def plot_water_images(self, wimages):
        file_out_base = os.path.join(os.path.dirname(self.file_nc), f'DayComparison_202304023_20230425')
        nimages = len(wimages)
        pm = PlotMultiple()
        nrow = nimages
        ncol = 2
        pm.start_multiple_plot_advanced(nrow, ncol, 5, 12, 0, 0.35, True)
        index_row = 0
        for wimage in wimages:
            pm.plot_image_title(wimages[wimage]['file_1'], index_row, 0, wimages[wimage]['title_1'])
            pm.plot_image_title(wimages[wimage]['file_2'], index_row, 1, wimages[wimage]['title_2'])
            index_row = index_row + 1
        pm.save_fig(f'{file_out_base}.{self.format_img}')
        pm.close_plot()

    def save_img_files(self, multiple_plot):
        flags = ['sky_irr_1', 'sky_rad_1', 'water_rad', 'sky_rad_2', 'sky_irr_2', 'sun']
        dir_img = os.path.join(os.path.dirname(self.file_nc), 'IMG')
        if not os.path.exists(dir_img):
            os.mkdir(dir_img)
        file_out_base = os.path.join(dir_img, f'CameraImages_{self.isequence}')
        if multiple_plot:
            pm = PlotMultiple()
            nrow = 2
            ncol = 3
            pm.start_multiple_plot_advanced(nrow, ncol, 10, 7.0, 0.1, 0.15, True)
        index_row = 0
        index_col = 0
        for flag in flags:
            file_img, title = self.get_img_file(flag)

            if index_col == ncol:
                index_col = 0
                index_row = index_row + 1

            if multiple_plot:
                if file_img is not None:
                    # pm.plot_image_title(file_img,index_row,index_col,title)
                    pm.plot_image_hypernets(file_img, index_row, index_col, title)
                else:
                    pm.plot_blank_with_title(index_row, index_col, title[:9])

            index_col = index_col + 1
        if multiple_plot:
            pm.save_fig(f'{file_out_base}_all{self.format_img}')
            pm.close_plot()

    def save_spectra_files(self, multiple_plot):
        dir_img = os.path.join(os.path.dirname(self.file_nc), 'IMG')
        if not os.path.exists(dir_img):
            os.mkdir(dir_img)
        file_out_base = os.path.join(dir_img, f'Spectra_{self.isequence}')
        flags = ['irradiance', 'downwelling_radiance', 'upwelling_radiance', 'water_leaving_radiance',
                 'reflectance_nosc', 'reflectance']
        if multiple_plot:
            pm = PlotMultiple()
            nrow = 2
            ncol = 3
            pm.start_multiple_plot_advanced(nrow, ncol, 10, 5.5, 0.25, 0.30, True)
        index_row = 0
        index_col = 0
        for flag in flags:
            if multiple_plot:
                if index_col == ncol:
                    index_col = 0
                    index_row = index_row + 1
                ax_here = pm.get_axes(index_row, index_col)
                self.plot_spectra(flag, ax_here)
                index_col = index_col + 1

        if multiple_plot:
            pm.save_fig(f'{file_out_base}_all{self.format_img}')
            pm.close_plot()

    def save_angle_files(self, multiple_plot):
        dir_img = os.path.join(os.path.dirname(self.file_nc), 'IMG')
        if not os.path.exists(dir_img):
            os.mkdir(dir_img)
        file_out_base = os.path.join(dir_img, f'Angles_{self.isequence}')
        flags = ['sza_nosc', 'saa_nosc', 'paa_nosc', 'sza', 'saa', 'paa']

        if multiple_plot:
            pm = PlotMultiple()
            nrow = 2
            ncol = 3
            pm.start_multiple_plot_polar(nrow, ncol, 10, 6, 0.25, 0.40, True)
        index_row = 0
        index_col = 0
        for flag in flags:
            if multiple_plot:
                if index_col == ncol:
                    index_col = 0
                    index_row = index_row + 1
                ax_here = pm.get_axes(index_row, index_col)
                self.plot_angle(flag, ax_here)
                index_col = index_col + 1
        if multiple_plot:
            pm.save_fig(f'{file_out_base}_all{self.format_img}')
            pm.close_plot()

    def save_report_summary_image(self, site, date_here, dir_img_summary, daily_sequences_summary):

        file_out = os.path.join(os.path.dirname(self.file_nc),
                                f'{site}_{date_here.strftime("%Y%m%d")}_DailySummary{self.format_img}')
        print(f'[INFO] Output file: {file_out}')


        ##TIME SERIES
        # start_multiple_plot_advanced(self, nrow, ncol, xfigsize, yfigsize, wspace, hspace, frameon)
        file_out_ts = os.path.join(dir_img_summary, f'TS{self.format_img}')
        pmts = PlotMultiple()
        pmts.start_multiple_plot_advanced(2, 1, 3, 3.5, 0, 0, True)
        pmts.plot_image(os.path.join(dir_img_summary, 'time_series_epsilon.tif'), 0, 0)
        pmts.plot_image(os.path.join(dir_img_summary, 'time_series_r800_nosc.tif'), 1, 0)
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

        line = f'Total sequences: {daily_sequences_summary["NTotal"]}/{daily_sequences_summary["expected_sequences"]}. Processed to L2: {daily_sequences_summary["NAvailable"]}.'
        skip = ['NTotal', 'NAvailable', 'start_time', 'end_time', 'expected_sequences']
        for key in daily_sequences_summary.keys():
            if key in skip: continue
            if daily_sequences_summary[key] == 0: continue
            line = f'{line} {key}: {daily_sequences_summary[key]}. '

        pmtop.set_text_size(-1750, 850, line, 8)
        pmtop.save_fig(file_out_top)
        pmtop.close_plot()

        ##MIDDLE(ANGLE-PANEL)
        file_out_middle = os.path.join(dir_img_summary, f'MIDDLE{self.format_img}')
        pmmiddle = PlotMultiple()
        pmmiddle.start_multiple_plot_advanced(2, 3, 10, 6, 0, 0, True)
        pmmiddle.plot_image(os.path.join(dir_img_summary, 'angle_800_nosc_sza.tif'), 0, 0)
        pmmiddle.plot_image(os.path.join(dir_img_summary, 'angle_800_nosc_saa.tif'), 0, 1)
        pmmiddle.plot_image(os.path.join(dir_img_summary, 'angle_800_nosc_paa.tif'), 0, 2)
        pmmiddle.plot_image(os.path.join(dir_img_summary, 'angle_800_nosc_sza.tif'), 1, 0)
        pmmiddle.plot_image(os.path.join(dir_img_summary, 'angle_800_nosc_saa.tif'), 1, 1)
        pmmiddle.plot_image(os.path.join(dir_img_summary, 'angle_800_nosc_paa.tif'), 1, 2)
        pmmiddle.save_fig(file_out_middle)
        pmmiddle.close_plot()

        ##DOWN(SPECTRA-PANEL)
        file_out_down = os.path.join(dir_img_summary, f'DOWN{self.format_img}')
        pdown = PlotMultiple()
        pdown.start_multiple_plot_advanced(2, 3, 10, 6, 0, 0, True)
        pdown.plot_image(os.path.join(dir_img_summary, 'downwelling_irradiance.tif'), 0, 0)
        pdown.plot_image(os.path.join(dir_img_summary, 'downwelling_radiance.tif'), 0, 1)
        pdown.plot_image(os.path.join(dir_img_summary, 'upwelling_radiance.tif'), 0, 2)
        pdown.plot_image(os.path.join(dir_img_summary, 'water_leaving_radiance.tif'), 1, 0)
        pdown.plot_image(os.path.join(dir_img_summary, 'reflectance_nosc.tif'), 1, 1)
        pdown.plot_image(os.path.join(dir_img_summary, 'reflectance.tif'), 1, 2)
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

        ##delete intermeddiate figures
        for name in os.listdir(dir_img_summary):
            os.remove(os.path.join(dir_img_summary, name))

        os.rmdir(dir_img_summary)

    def save_report_image_only_pictures(self,site,delete_images,overwrite,seq,files_img):
        print(f'[INFO] Sequence {seq} (No Level-2 data available)')
        seq_time_str = seq[3:]
        seq_time = dt.strptime(seq_time_str,'%Y%m%dT%H%M')
        file_out = os.path.join(os.path.dirname(self.file_nc),f'{site}_{seq_time_str}_Report{self.format_img}')
        if os.path.exists(file_out) and not overwrite:
            return
        names_img = ['sky_irr_1', 'sky_rad_1', 'water_rad', 'sky_rad_2', 'sky_irr_2', 'sun']
        names_img_files = {}
        for ref in files_img:
            name_img = files_img[ref]['name_img']
            if name_img in names_img:
                names_img_files[name_img] = files_img[ref]['file_img']

        dir_img = os.path.join(os.path.dirname(self.file_nc), 'IMG')
        if not os.path.exists(dir_img):
            os.mkdir(dir_img)

        pm = PlotMultiple()
        nrow = 2
        ncol = 3
        pm.start_multiple_plot_advanced(nrow, ncol, 10, 7.0, 0.02,0.15, True)
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
        line = f'ANOMALY: '
        pm.fig.text(0.20,0.05,line)
        pm.save_fig(file_out)
        pm.close_plot()

        if delete_images:
            for name in os.listdir(dir_img):
                file_here = os.path.join(dir_img, name)
                os.remove(file_here)
            os.rmdir(dir_img)


    def save_report_image(self, site, delete_images, overwrite):
        print(f'[INFO] Sequence {self.isequence}: SEQ{self.sequences[self.isequence]}')
        if self.sequences[self.isequence] is None:
            return
        file_out = os.path.join(os.path.dirname(self.file_nc),
                                f'{site}_{self.sequences[self.isequence]}_Report{self.format_img}')
        if os.path.exists(file_out) and not overwrite:
            print(f'[WARNING] File {file_out} already exists. Skipping...')
            return
        dir_img = os.path.join(os.path.dirname(self.file_nc), 'IMG')
        if not os.path.exists(dir_img):
            os.mkdir(dir_img)
        file_angle = os.path.join(dir_img, f'Angles_{self.isequence}_all{self.format_img}')
        if not os.path.isfile(file_angle) or (os.path.isfile(file_angle) and overwrite):
            self.save_angle_files(True)
        file_img = os.path.join(dir_img, f'CameraImages_{self.isequence}_all{self.format_img}')
        if not os.path.isfile(file_img) or (os.path.isfile(file_img) and overwrite):
            self.save_img_files(True)
        file_spectra = os.path.join(dir_img, f'Spectra_{self.isequence}_all{self.format_img}')
        if not os.path.exists(file_spectra) or (os.path.isfile(file_spectra) and overwrite):
            self.save_spectra_files(True)
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
        pm.plot_image(file_angle, 1, 0)
        info_l2 = self.get_info_l2()
        pm.get_axes(2, 0).set_title(info_l2, fontsize=12)
        pm.plot_image(file_spectra, 2, 0)

        pm.save_fig(file_out)
        pm.close_plot()

        if delete_images:
            for name in os.listdir(dir_img):
                file_here = os.path.join(dir_img, name)
                os.remove(file_here)
            os.rmdir(dir_img)

    def plot_angle(self, flag, ax_here):

        angle_flag = flag.split('_')[0]
        apply_nosc = False
        if flag.endswith('_nosc'):
            apply_nosc = True

        reflectance_ref = self.reflectance_ref(apply_nosc)
        reflectance_ref_l1 = self.reflectance_ref_l1(apply_nosc)
        angle_rad, angle_label = self.angle_rad(angle_flag)
        angle_rad_l1 = self.angle_rad_level1(angle_flag)

        reflectance_ref[reflectance_ref < 1e-5] = 1.1e-5
        reflectance_ref_l1[reflectance_ref_l1 < 1e-5] = 1.1e-5

        if self.valid_sequences is None:
            self.get_valid_flags()

        ax_here.set_rscale('log')
        ax_here.set_rlim((1e-5, 1))
        ax_here.set_rticks([1e-4, 1e-3, 1e-2, 1e-1, 1])
        ax_here.set_theta_zero_location("N")
        ax_here.set_theta_direction(-1)

        ##valid
        ax_here.scatter(angle_rad[self.valid_sequences == 0], reflectance_ref[self.valid_sequences == 0], marker='o',
                        s=8, color='green')
        ##qc_flagedd
        ax_here.scatter(angle_rad[self.valid_sequences == 1], reflectance_ref[self.valid_sequences == 1], marker='o',
                        s=8, color='red')

        ##ENEG
        ax_here.scatter(angle_rad[self.valid_sequences == 2], reflectance_ref[self.valid_sequences == 2], marker='o',
                        s=8, color='cyan')

        ##EHIGH
        ax_here.scatter(angle_rad[self.valid_sequences == 3], reflectance_ref[self.valid_sequences == 3], marker='o',
                        s=8, color='magenta')

        ax_here.scatter(angle_rad_l1[:], reflectance_ref_l1[:], marker='s', s=8, color='gray')
        ax_here.scatter(angle_rad[self.isequence], reflectance_ref[self.isequence], marker='s', s=12, color='blue')

        val = (angle_rad[self.isequence] * 180) / np.pi
        if apply_nosc:
            title = f'{self.rho}({self.ref_wl})/NOSC vs.{angle_flag} ({val:.2f}°)'
        else:
            title = f'{self.rho}({self.ref_wl}) vs.{angle_flag} ({val:.2f}°)'
        ax_here.set_title(title, fontsize=10)

        if angle_flag == 'sza':
            ax_here.set_thetamin(0)
            ax_here.set_thetamax(90)
            ax_here.set_xticks((np.pi * np.array([0, 15, 30, 45, 60, 75, 90]) / 180))
            ax_here.set_rticks([1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1])
        if angle_flag == 'saa' or angle_flag == 'paa':
            ax_here.set_rlabel_position(angle_label)

        ax_here.tick_params(axis='x', labelsize=10)
        ax_here.tick_params(axis='y', labelsize=10)

    def plot_spectra(self, flag, ax_here):

        info = {
            'irradiance': {
                'ylabel': 'Ed',
                'title': 'Downwelling Irradiance'
            },
            'downwelling_radiance': {
                'ylabel': 'Li',
                'title': 'Downwelling radiance'
            },
            'upwelling_radiance': {
                'ylabel': 'Lt',
                'title': 'Upwelling radiance'
            },
            'water_leaving_radiance': {
                'ylabel': 'Lw',
                'title': 'Water-leaving radiance'
            },
            'reflectance': {
                'ylabel': r'ρ$_w$',
                'title': 'Reflectance'
            },
            'reflectance_nosc': {
                'ylabel': r'ρ$_w$',
                'title': 'Reflectance Nosc'
            },
        }
        l1_variable = f'l1_{flag}'
        l2_variable = f'l2_{flag}'
        dataset = Dataset(self.file_nc)

        spectra = np.array(dataset.variables[l1_variable][self.isequence, :, :]).transpose()
        spectra_l2 = None
        if l2_variable in dataset.variables:
            spectra_l2 = np.array(dataset.variables[l2_variable][self.isequence, :]).transpose()
        wavelength = np.array(dataset.variables['wavelength'])
        for ispectra in range(spectra.shape[0]):
            ax_here.plot(wavelength, spectra[ispectra, :], color='gray', linewidth=0.5)
        if spectra_l2 is not None:
            ax_here.plot(wavelength, spectra_l2, color='black', linewidth=0.25)

        ax_here.set_xlabel('Wavelength(nm)', fontsize=7)
        ax_here.set_ylabel(info[flag]['ylabel'], fontsize=7)
        ax_here.tick_params(axis='x', labelsize=7)
        ax_here.tick_params(axis='y', labelsize=7)
        ax_here.set_title(info[flag]['title'], fontsize=7)
        ax_here.grid(which='major', color='lightgray', linestyle='--', axis='y')
        dataset.close()

    def check_img_files(self):
        flags = ['sky_irr_1', 'sky_irr_2', 'sky_rad_1', 'sky_rad_1', 'water_rad', 'sun']
        for flag in flags:
            file_img, title = self.get_img_file(flag)
            print(f'[INFO] Sequence: {self.sequences[self.isequence]}')
            if file_img is None:
                print(f'[INFO] {flag}-> N/Av')
            else:
                print(f'[INFO] {flag}->{file_img}->{os.path.exists(file_img)}')

    def start_dataset_w(self, file_out):

        dataset_w = Dataset(file_out, 'w', format='NETCDF4')
        dataset_r = Dataset(self.file_nc)
        # copy attributes
        dataset_w.setncatts(dataset_r.__dict__)

        # copy dimensions
        for name, dimension in dataset_r.dimensions.items():
            dataset_w.createDimension(name, (len(dimension) if not dimension.isunlimited() else None))

        # copy all file data except for the excluded
        for name, variable in dataset_r.variables.items():
            dataset_w.createVariable(name, variable.datatype, variable.dimensions)
            # copy variable attributes all at once via dictionary
            dataset_w[name].setncatts(dataset_r[name].__dict__)

            if name == 'wavelength' or name == 'bandwidth':
                dataset_w[name][:] = dataset_r[name][:]

        dataset_r.close()
        return dataset_w

    def add_flag_variables(self, dataset_w, flags):
        for flag in flags:
            var_name = flag
            if not var_name.startswith('flag_'):
                var_name = f'flag_{flag}'
            var = dataset_w.createVariable(var_name, 'i8', ('series',))
            var[:] = 0
            var.flag_values = flags[flag]['values']
            var.flag_meanings = ' '.join(flags[flag]['meanings'])
        return dataset_w

    def set_data_dataset_w(self, dataset_w, sindices, index_w):

        dini = index_w
        dfin = index_w + len(sindices)
        dataset_r = Dataset(self.file_nc)
        for variable in dataset_w.variables:
            if variable == 'wavelength' or variable == 'bandwidth':
                continue
            if variable not in dataset_r.variables:  ##for flag variables
                continue
            ndim = len(dataset_w[variable].shape)
            if ndim == 1:
                dataset_w[variable][dini:dfin] = dataset_r[variable][sindices]
            elif ndim == 2:
                dataset_w[variable][dini:dfin, :] = dataset_r[variable][sindices, :]
            elif ndim == 3:
                dataset_w[variable][dini:dfin, :, :] = dataset_r[variable][sindices, :, :]

        dataset_r.close()
        return dataset_w

    def set_data_flag(self, dataset_w, index_w, flag, flag_array):
        dini = index_w
        dfin = index_w + len(flag_array)
        var_name = flag
        if not var_name.startswith('flag_'):
            var_name = f'flag_{flag}'
        if var_name in dataset_w.variables:
            dataset_w[var_name][dini:dfin] = flag_array[:]
        return dataset_w

    def get_csv_col_names(self):

        col_names = ['sequence_ref']
        dataset_r = Dataset(self.file_nc)
        for var in dataset_r.variables:
            if var.startswith('l2'):
                ndim = len(dataset_r.variables[var].shape)
                if ndim == 1:
                    col_names.append(var[3:])

        col_names = col_names + ['rhow_nosc_800', 'rhow_800', 'rhow_nosc_350_450', 'rhow_350_450', 'isequence',
                                 'file_nc']

        dataset_r.close()

        return col_names

    def get_dataframe_lines(self, sindices, col_names):

        if self.sequences is None:
            self.sequences = self.get_sequences()
        dataset_r = Dataset(self.file_nc)
        data = {}
        for c in col_names:
            if c == 'file_nc':
                data[c] = [self.file_nc] * len(sindices)
            elif c == 'isequence':
                data[c] = sindices
            elif c == 'sequence_ref':
                data[c] = [self.sequences[idx] for idx in sindices]
            elif c.startswith('rhow'):
                lc = c.split('_')
                if lc[1] == 'nosc':
                    wl_ini = float(lc[2])
                    variable = 'l2_reflectance_nosc'
                else:
                    wl_ini = float(lc[1])
                    variable = 'l2_reflectance'
                wl_fin = float(lc[-1])
                wls = np.array(dataset_r.variables['wavelength'][:])
                i_ini = np.argmin(np.abs(wl_ini - wls))
                i_fin = np.argmin(np.abs(wl_fin - wls)) + 1
                reflectance = dataset_r.variables[variable][sindices, i_ini:i_fin]
                data[c] = np.mean(reflectance, axis=1)
            else:
                variable = f'l2_{c}'
                data[c] = np.array(dataset_r.variables[variable][sindices])

        dataset_r.close()

        df = pd.DataFrame(data)
        return df

    def plot_from_options_impl(self, options_figure):
        if not options_figure['apply']:
            return
        if options_figure['type'] == 'spectraplot' and options_figure['type_rrs'] == 'user_defined':
            print(f'[INFO] Spectra plot')
            self.plot_spectra_plot_from_options(options_figure)
        if options_figure['type'] == 'timeseries':
            print(f'[INFO] Time series plot')
            self.plot_time_series_from_options(options_figure)
        if options_figure['type'] == 'sequence':
            print(f'[INFO] Sequence plot')
            daily_sequence_summary =  self.plot_sequence_plot_from_options(options_figure)
            return daily_sequence_summary
        if options_figure['type'] == 'flagplot' and options_figure['type_flagplot'] == 'comparison':
            print(f'[INFO] Flag plot')
            self.plot_flag_plot_comparison(options_figure)
        if options_figure['type'] == 'angleplot':
            print(f'[INFO] Angle plot')
            self.plot_angle_plot_from_options(options_figure)

        return None
    def get_sequence_info(self):

        dataset = Dataset(self.file_nc)
        qf_array = dataset.variables['l2_quality_flag'][:]
        epsilon_array = dataset.variables['l2_epsilon'][:]
        ndata = len(qf_array)
        qf_flags_meanings = dataset.variables['l2_quality_flag'].flag_meanings
        qf_flags_list = qf_flags_meanings.split(' ')
        qf_flags_values = [np.uint64(x) for x in dataset.variables['l2_quality_flag'].flag_masks.split(',')]
        fcheck = Flags(qf_flags_values, qf_flags_meanings)
        nflag_total = 0
        for index, value in enumerate(qf_flags_values):
            mask = fcheck.Mask(qf_array, [qf_flags_list[index]])
            mask = np.ma.filled(mask, 0)
            nflag = np.count_nonzero(mask > 0)
            if nflag > 0:
                nflag_total = nflag_total + nflag

        eneg = np.logical_and(qf_array == 0, epsilon_array < (-0.05))
        ehigh = np.logical_and(qf_array == 0, epsilon_array >= 0.05)
        valid = np.logical_and(qf_array == 0, np.logical_and(epsilon_array >= (-0.05), epsilon_array < 0.05))
        nvalid = int(np.count_nonzero(valid))
        nepsilon = int(np.count_nonzero(eneg) + np.count_nonzero(ehigh))
        navailable = nvalid + nepsilon + nflag_total
        daily_sequences_summary = {
            'Total': ndata,
            'Available': navailable,
            'Valid': nvalid,
            'QFlagged': nflag_total,
            'EFlagged': nepsilon
        }
        dataset.close()
        return daily_sequences_summary

    def plot_flag_plot_comparison(self, options_figure):

        dataset = Dataset(self.file_nc)
        qf_array = dataset.variables['l2_quality_flag'][:]
        epsilon_array = dataset.variables['l2_epsilon'][:]
        ndata = len(qf_array)
        qf_flags_meanings = dataset.variables['l2_quality_flag'].flag_meanings
        qf_flags_list = qf_flags_meanings.split(' ')
        qf_flags_values = [np.uint64(x) for x in dataset.variables['l2_quality_flag'].flag_masks.split(',')]
        ylabel = ['VALID'] + qf_flags_list + ['ENEG', 'EHIGH']
        yarray = np.zeros((len(ylabel)))
        fcheck = Flags(qf_flags_values, qf_flags_meanings)
        nflag_total = 0
        for index, value in enumerate(qf_flags_values):
            mask = fcheck.Mask(qf_array, [qf_flags_list[index]])
            mask = np.ma.filled(mask, 0)
            nflag = np.count_nonzero(mask > 0)
            if nflag > 0:
                yarray[index + 1] = nflag
                nflag_total = nflag_total + nflag

        eneg = np.logical_and(qf_array == 0, epsilon_array < (-0.05))
        ehigh = np.logical_and(qf_array == 0, epsilon_array >= 0.05)
        valid = np.logical_and(qf_array == 0, np.logical_and(epsilon_array >= (-0.05), epsilon_array < 0.05))
        yarray[-2] = np.count_nonzero(eneg)
        yarray[-1] = np.count_nonzero(ehigh)
        yarray[0] = np.count_nonzero(valid)
        xarray = np.arange(len(yarray))
        plt.barh(xarray, yarray)
        xarray = np.append(xarray, xarray[-1] + 1)
        ylabel.append('')
        plt.yticks(xarray - 0.5, ylabel, fontsize=8, verticalalignment='bottom')
        if ndata <= 5:
            xticks = [0, 1, 2, 3, 4, 5]
        elif ndata >= 5 and ndata < 10:
            xticks = [0, 2, 4, 6, 8, 10]
        else:
            xmax = 5 * (np.floor(ndata / 5) + 1)
            xticks = np.arange(0, xmax + 1, 5).tolist()
        plt.xticks(xticks, [f'{x:.0f}' for x in xticks])
        plt.xlabel(f'Number of sequences')
        plt.grid()
        plt.gcf().tight_layout()
        file_out = options_figure['file_out']
        if file_out is not None:
            plt.savefig(file_out, dpi=300)

        dataset.close()

    def plot_sequence_plot_from_options(self, options_figure):

        dataset = Dataset(self.file_nc)
        time_array = dataset.variables['l2_acquisition_time'][:]
        time_array = np.ma.masked_values(time_array, 0)  ##solving a problem find 0n 30/05/2023, it shouln't be happen
        start_time_real = dt.utcfromtimestamp(np.min(time_array))
        end_time_real = dt.utcfromtimestamp(np.max(time_array))
        if start_time_real.strftime('%Y%m%d') != end_time_real.strftime('%Y%m%d'):
            dataset.close()
            print('[ERROR] Plot is only created for a single day')
            return

        # print(options_figure)

        use_default_flag = True
        if options_figure['flagBy'] is not None:
            # print('===========================================================================================================')
            options_figure['flagType'] = 'flag'
            options_figure = self.check_gs_options_impl(options_figure, 'flagBy', 'flagType', 'flagValues')
            if options_figure['flagBy'] in options_figure.keys() and len(
                    options_figure[options_figure['flagBy']]['flag_values']) == 4:
                use_default_flag = False
                flag_meanings = options_figure[options_figure['flagBy']]['flag_meanings']
                flag_values = options_figure[options_figure['flagBy']]['flag_values']
                flag_array = options_figure[options_figure['flagBy']]['flag_array']
                # print(options_figure[options_figure['flagBy']])
            # print('********************************************************')
        if use_default_flag:
            qf_array = dataset.variables['l2_quality_flag'][:]
            epsilon_array = dataset.variables['l2_epsilon'][:]
            flag_meanings = ['FLAGGED', 'ENEG', 'EHIGH', 'VALID']

        time_fix_axis = self.get_fix_axis_time(options_figure, start_time_real, end_time_real)
        ntime = len(time_fix_axis)
        time_fix_min_max = np.zeros((ntime, 2))
        seconds_ref = self.get_time_interval_seconds(options_figure['frequency'], 'minutes')
        time_fix_axis_ts = np.array([x.replace(tzinfo=pytz.utc).timestamp() for x in time_fix_axis]).astype(np.float64)
        time_fix_min_max[:, 0] = time_fix_axis_ts - seconds_ref
        time_fix_min_max[:, 1] = time_fix_axis_ts + seconds_ref
        xarray = np.arange(ntime)
        yarray = np.zeros(xarray.shape)
        hours_ticks = []
        minutes_ticks = []
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
            daily_summary_sequences[f] = 0

        for itime in range(ntime):
            time_valid = np.logical_and(time_array >= time_fix_min_max[itime, 0],
                                        time_array < time_fix_min_max[itime, 1])
            htick = time_fix_axis[itime].strftime('%H')
            mtick = time_fix_axis[itime].strftime('%M')

            if htick not in hours_ticks: hours_ticks.append(htick)
            if mtick not in minutes_ticks: minutes_ticks.append(mtick)

            if np.count_nonzero(time_valid) == 1:
                daily_summary_sequences['NAvailable'] = daily_summary_sequences['NAvailable'] + 1
                # print('time_valid_good-->',time_valid,time_valid.shape,)
                if use_default_flag:
                    qf_value = qf_array[time_valid][0]
                    epsilon_value = epsilon_array[time_valid][0]
                    if qf_value == 0:
                        if epsilon_value < (-0.05):
                            yarray[itime] = 2
                            daily_summary_sequences['ENEG'] = daily_summary_sequences['ENEG'] + 1
                        elif (-0.05) <= epsilon_value < 0.05:
                            daily_summary_sequences['VALID'] = daily_summary_sequences['VALID'] + 1
                            yarray[itime] = 3
                        elif epsilon_value > 0.05:
                            daily_summary_sequences['EHIGH'] = daily_summary_sequences['EHIGH'] + 1
                            yarray[itime] = 4
                    else:
                        daily_summary_sequences['FLAGGED'] = daily_summary_sequences['FLAGGED']+1
                        yarray[itime] = 1
                else:
                    try:
                        fvalue = flag_array[time_valid]
                        index_flag = flag_values.index(int(fvalue))
                        flag_meaning = legend_values[index_flag]
                        daily_summary_sequences[flag_meaning] = daily_summary_sequences[flag_meaning]+1
                        yarray[itime] = index_flag + 1
                    except:
                        pass

        hours_ticks.reverse()
        plt.Figure()
        data = pd.DataFrame(index=hours_ticks, columns=minutes_ticks).astype(np.float64)
        yarray[yarray == 0] = np.nan
        for itime, tf in enumerate(time_fix_axis):
            data.loc[tf.strftime('%H')].at[tf.strftime('%M')] = yarray[itime]


        if self.sequences_no_data is not None and len(self.sequences_no_data) > 0:
            daily_summary_sequences['NTotal'] = daily_summary_sequences['NAvailable'] + len(self.sequences_no_data)
            for seq in self.sequences_no_data:
                time_stamp_seq = dt.strptime(seq[3:], '%Y%m%dT%H%M').replace(tzinfo=pytz.utc).timestamp()
                iwhere = np.where(np.logical_and(time_stamp_seq>=time_fix_min_max[:,0],time_stamp_seq<time_fix_min_max[:,1]))
                if len(iwhere[0])>0:
                    index = iwhere[0][0]
                    date_seq = dt.utcfromtimestamp(time_fix_axis_ts[index])
                    data.loc[date_seq.strftime('%H')].at[date_seq.strftime('%M')] = 0
        else:
            daily_summary_sequences['NTotal'] = daily_summary_sequences['NAvailable']

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
        plt.xlabel('Minutes')
        plt.ylabel('Hours')

        if options_figure['title'] is not None:
            title = options_figure['title']
            title = title.replace('$DATE$', start_time_real.strftime('%Y-%m-%d'))
            plt.title(title, fontsize=options_figure['fontsizetitle'])

        plt.gcf().tight_layout()
        file_out = options_figure['file_out']
        if file_out is not None:
            plt.savefig(file_out, dpi=300)

        plt.close()


        dataset.close()

        return daily_summary_sequences

    def plot_angle_plot_from_options(self, options_figure):

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
        nseries = angle_array.shape[0]
        nscan = 1
        if options_figure['angle_var'].startswith('l1'):
            nscan = angle_array.shape[2]
            angle_array = self.reduce_l1_dimensions(angle_array)

        avg_array = dataset.variables[options_figure['avg_var']][:]
        avg_array = np.ma.filled(avg_array, -999.0)
        if options_figure['avg_var'].startswith('l1'):
            if len(avg_array.shape) == 2: is_spectral = False
            if len(avg_array.shape) == 3: is_spectral = True

        if options_figure['avg_var'].startswith('l2'):
            if len(avg_array.shape) == 1: is_spectral = False
            if len(avg_array.shape) == 2: is_spectral = True

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
            if len(avg_array.shape) == 3: avg_array = avg_array[:, index_ref, :]
            if len(avg_array.shape) == 2: avg_array = avg_array[:, index_ref]

        if nscan > 1 and len(avg_array.shape) == 1: avg_array = self.multiply_array_by_scan(avg_array, nseries, nscan)
        if nscan > 1 and len(avg_array.shape) == 2: avg_array = self.reduce_l1_dimensions(avg_array)

        ngroups, groupValues, groupArray, str_legend = self.check_group_and_legend(options_figure, nseries, nscan)

        ##PLOTTING
        fig, ax_here = plt.subplots(figsize=(3, 3), subplot_kw={'projection': 'polar'})
        ax_here = plt.gca()
        ax_here.set_rscale(options_figure['scale'])
        ax_here.set_rlim(options_figure['rlim'])
        ax_here.set_rticks(options_figure['rticks'])
        ax_here.set_theta_zero_location(options_figure['theta_zero_location'])
        ax_here.set_theta_direction(options_figure['theta_direction'])

        valid_array = np.logical_and(angle_array != -999.0, avg_array != -999.0)
        angle_array = angle_array[valid_array]
        avg_array = avg_array[valid_array]
        if ngroups > 1 and groupArray is not None:
            groupArray = groupArray[valid_array]
        if options_figure['min_data'] is not None:
            avg_array[avg_array < options_figure['min_data']] = options_figure['min_data'] + 1e-6

        if ngroups > 1 and groupArray is not None:
            colors = options_figure['color']
            if len(colors) != ngroups:
                colors = pdefaults.get_color_list(ngroups)
            point_size = options_figure['point_size']
            point_marker = options_figure['point_marker']
            if len(point_size) != ngroups:
                point_size = [point_size[0]] * ngroups
            if len(point_marker) != ngroups:
                point_marker = [point_marker[0]] * ngroups
            for ig, g in enumerate(groupValues):
                valid_g = groupArray == g
                if np.count_nonzero(valid_g) > 0:
                    angle_array_rad = (angle_array[valid_g] * np.pi) / 180
                    avg_array_here = avg_array[valid_g]
                    ax_here.scatter(angle_array_rad, avg_array_here, marker=point_marker[ig], s=point_size[ig],
                                    color=colors[ig])
        else:
            angle_array_rad = (angle_array * np.pi) / 180
            ax_here.scatter(angle_array_rad, avg_array, marker=options_figure['point_marker'][0],
                            s=options_figure['point_size'][0], color=options_figure['color'][0])

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

    def reduce_l1_dimensions(self, array):
        if len(array.shape) == 3:  ##spectral variables
            nspectra_total = array.shape[2] * array.shape[0]
            new_shape = (nspectra_total, array.shape[1])
            array_new = np.zeros(new_shape)
            for iscan in range(array.shape[2]):
                iini = iscan * array.shape[0]
                ifin = iini + array.shape[0]
                array_new[iini:ifin, :] = array[:, :, iscan]
            return array_new

        if len(array.shape) == 2:  ##no spectral variables
            ntotal = array.shape[0] * array.shape[1]
            new_shape = (ntotal,)
            array_new = np.zeros(new_shape)
            for iscan in range(array.shape[1]):
                iini = iscan * array.shape[0]
                ifin = iini + array.shape[0]
                array_new[iini:ifin] = array[:, iscan]
            return array_new

    def multiply_by_scan(self, nseries, nscan, sequences_here, sequence_indices):
        nhere = len(sequences_here)
        ntotal = nhere * nscan
        sequence_here_total = [None] * ntotal
        sequence_indices_total = [None] * ntotal
        for iscan in range(nscan):
            iini = iscan * nhere
            ifin = iini + nhere
            sequence_here_total[iini:ifin] = sequences_here[:]
            sequence_indices_total[iini:ifin] = sequence_indices[:] + (iscan * nseries)

        return sequence_here_total, sequence_indices_total

    def multiply_array_by_scan(self, array, nseries, nscan):
        if array.shape[0] != nseries:
            return None
        ntotal = nseries * nscan
        array_new = np.zeros((ntotal,), dtype=array.dtype)
        for iscan in range(nscan):
            iini = iscan * nseries
            ifin = iini + nseries
            array_new[iini:ifin] = array[:]
        return array_new

    ##spectra without _fILLValue
    def get_spectra_stats(self, spectra_good):

        spectra_avg = np.mean(spectra_good, axis=0)
        spectra_std = np.std(spectra_good, axis=0)
        indices_max = np.argmax(spectra_good, axis=0)
        imax = st.mode(indices_max)
        spectra_max_real = spectra_good[imax, :]
        spectra_max = np.max(spectra_good, axis=0)

        indices_min = np.argmin(spectra_good, axis=0)
        imin = st.mode(indices_min)
        spectra_mim_real = spectra_good[imin, :]
        spectra_min = np.min(spectra_good, axis=0)

        spectra_median = np.median(spectra_good, axis=0)

        spectra_p25 = np.percentile(spectra_good, 25, axis=0)
        spectra_p75 = np.percentile(spectra_good, 75, axis=0)

        spectra_stats = {
            'avg': spectra_avg,
            'std': spectra_std,
            'spectra_min_real': spectra_mim_real,
            'spectra_max_real': spectra_max_real,
            'spectra_min': spectra_min,
            'spectra_max': spectra_max,
            'median': spectra_median,
            'p25': spectra_p25,
            'p75': spectra_p75
        }
        return spectra_stats

    def plot_time_series_from_options(self, options_figure):

        dataset = Dataset(self.file_nc)
        time_var = options_figure['time_var']
        avg_vars = options_figure['avg_var']
        nseries = 0
        nscan = 1
        if time_var is None or avg_vars is None:
            print(f'[ERROR] time_var and avg_var should be defined in the configuration file')
            return
        if time_var not in dataset.variables:
            print(f'[ERROR] {time_var} is not defined in {self.file_nc}')
            return
        else:
            nseries = dataset.variables[time_var].shape[0]

        is_spectral = []
        for avg_var in avg_vars:
            if avg_var not in dataset.variables:
                print(f'[ERROR] {avg_var} is not defined in {self.file_nc}')
                return
            else:
                sh = dataset.variables[avg_var].shape
                if avg_var.startswith('l1'):
                    nscan = sh[-1]
                    if len(sh) == 2:
                        is_spectral.append(False)
                    elif len(sh) == 3:
                        is_spectral.append(True)
                if avg_var.startswith('l2'):
                    if len(sh) == 1:
                        is_spectral.append(False)
                    elif len(sh) == 2:
                        is_spectral.append(True)

        options_figure['wlref_index'] = -1
        if options_figure['wlref'] != -999.0 and is_spectral.count(True) > 0:
            wavelength = dataset.variables['wavelength'][:]
            index_ref = np.argmin(np.abs(options_figure['wlref'] - wavelength))
            if abs(wavelength[index_ref] - options_figure['wlref']) <= 2:
                options_figure['wlref_index'] = index_ref

        options_figure = self.check_gs_options_impl(options_figure, 'groupBy', 'groupType', 'groupValues')
        time_array = dataset.variables[time_var][:]
        time_array = np.ma.filled(time_array.astype(np.float64), -999.0)

        if nscan > 1:
            time_array = self.reduce_l1_dimensions(time_array)

        ngroups, groupValues, groupArray, str_legend = self.check_group_and_legend(options_figure, nseries, nscan)
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
        elif ngroups > 1:  ##multiple groups, only one variable
            avg_var = avg_vars[0]
            avg_array = dataset.variables[avg_var][:]
            avg_array = np.ma.filled(avg_array.astype(np.float64), -999.0)
            if is_spectral[0]:
                index_ref = options_figure['wlref_index']
                if len(avg_array.shape) == 3:
                    avg_array = avg_array[:, index_ref, :]
                elif len(avg_array.shape) == 2:
                    avg_array = avg_array[:, index_ref]
                avg_array = np.squeeze(avg_array)
            print(f'[INFO] NSeries:  {nseries} NScans: {nscan} Variablel: {avg_var} Array shape: {avg_array.shape}')
            if nscan > 1 and nseries > 1 and len(avg_array.shape) == 1:
                avg_array = self.multiply_array_by_scan(avg_array, nseries, nscan)
            if nscan > 1 and len(avg_array.shape) == 2:
                avg_array = self.reduce_l1_dimensions(avg_array)
            handles = []
            str_legend_valid = []
            for icolor, gvalue in enumerate(groupValues):
                color = options_figure['color'][icolor]
                legend_added = False
                for idx in range(ntime):

                    xdata, ydata = self.get_data_nospectral_fix_time(options_figure, time_fix_min_max[idx, 0],
                                                                     time_fix_min_max[idx, 1], time_array,
                                                                     avg_array, idx, groupArray, gvalue)
                    if xdata is not None and ydata is not None:
                        h = pspectra.plot_single_marker(xdata, ydata, 'o', 6, color, None, 0)
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

    def get_data_nospectral_fix_time(self, options_figure, time_min, time_max, time_array, var_array, output_value,
                                     groupArray, groupValue):

        if groupArray is not None:
            valid_time = np.logical_and(np.logical_and(time_array >= time_min, time_array < time_max),
                                        groupArray == groupValue)
        else:
            valid_time = np.logical_and(time_array >= time_min, time_array < time_max)
        nvalid = np.count_nonzero(valid_time)
        if nvalid == 0:
            return None, None

        xdata = np.array([output_value] * nvalid).astype(np.float64)
        ydata = var_array[valid_time]

        return xdata, ydata

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

    def plot_spectra_plot_from_options(self, options_figure):

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
        ndim = len(var_y.shape)
        if '_FillValue' in var_y.ncattrs():
            fill_value = var_y._FillValue
        else:
            fill_value = default_fillvals[var_y.dtype.str[1:]]
        nscan = 1
        nseries = 0
        if ndim == 3:  # l1 variable
            spectra = np.array(dataset.variables[y_variable][:, :, :])
            nseries = spectra.shape[0]
            nscan = spectra.shape[2]
            spectra = self.reduce_l1_dimensions(spectra)
        elif ndim == 2:
            spectra = np.array(dataset.variables[y_variable][:, :])
            nseries = spectra.shape[0]

        ##Checking spectra with fill values
        spectra_check = np.where(spectra == fill_value, 1, 0)
        spectra_check = np.sum(spectra_check, axis=1)

        sequences_here, sequence_indices = self.get_sequences_range(options_figure['start_date'],
                                                                    options_figure['end_date'],
                                                                    options_figure['start_time'],
                                                                    options_figure['end_time'])
        if ndim == 3 and sequence_indices is not None:
            sequences_here, sequence_indices = self.multiply_by_scan(nseries, nscan, sequences_here, sequence_indices)

        if sequence_indices is not None:
            spectra_check[sequence_indices] = spectra_check[sequence_indices] + 100
            spectra_check = np.where(spectra_check == 100, 0, 1)

        spectra = spectra[spectra_check == 0, :]

        dataset.close()
        ##DATA SELECTION

        ##CHECKING GROUP AND LEGEND
        ngroup, groupValues, groupArray, str_legend = self.check_group_and_legend(options_figure, nseries, nscan)
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

    def get_gs_array(self, options_figure, by, type):

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

    def check_group_and_legend(self, options_figure, nseries, nscan):

        ngroup = 1
        groupValues = None
        groupArray = None
        if 'groupValues' in options_figure.keys():
            groupValues = options_figure['groupValues']
        if groupValues is not None:
            ngroup = len(groupValues)
        str_legend = []
        handles = []
        str_legend_valid = []
        if ngroup > 1 and options_figure['legend']:
            str_legend = self.get_str_legend(options_figure)

        if ngroup > 1:
            groupArray, all_flag_values, all_flag_meanings = self.get_gs_array(options_figure,
                                                                               options_figure['groupBy'],
                                                                               options_figure['groupType'])
            if nscan > 1:
                if len(groupArray.shape) == 1 and groupArray.shape[0] == nseries:
                    groupArray = self.multiply_array_by_scan(groupArray, nseries, nscan)
                if len(groupArray.shape) == 2 and groupArray.shape[0] == nseries and groupArray.shape[1] == nscan:
                    groupArray = self.reduce_l1_dimensions(groupArray)

            if 'color' in options_figure.keys():
                if len(options_figure['color']) != ngroup:
                    options_figure['color'] = pdefaults.get_color_list(ngroup)

        return ngroup, groupValues, groupArray, str_legend

    def check_gs_options_impl(self, options_figure, by, type, values):

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

                array, flag_meanings, flag_values = self.flag_builder.create_flag_array_ranges_v2(
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

    def get_color_default(self, value, min, max):
        nvalues = (max - min) + 1
        colors_default = ['Blue', 'Red', 'Green', 'm', 'Cyan', 'Orange', 'Yellow']
        if nvalues < 6:
            index = value - min
            return colors_default[index]

        cm = mpl.colormaps['jet']
        return cm((value - min) / (max - min))
