import numpy as np
from matplotlib import pyplot as plt


class PlotSpectra():

    def __init__(self):
        self.start_plot()
        self.xdata = []

        line_style_default = {
            'color': 'b',
            'marker': None,
            'linestyle': 'solid',
            'linewidth': 1,
            'markersize': 25
        }
        self.legend_options = {
            'loc': 'upper left',
            'bbox_to_anchor': (1.0, 1.0),
            'framealpha': 0.8,
            'ncols': 1,
            'markerscale': 1
        }
        self.legend_options_bottom = {
            'loc': 'lower center',
            'bbox_to_anchor': (0.5, -0.25),
            'framealpha': 0.8,
            'ncols': 1,
            'markerscale': 1
        }

        self.line_style_default = line_style_default.copy()
        self.spectra_style = line_style_default.copy()
        self.spectra_style['color'] = 'gray'

        self.stats_style = {
            'central': line_style_default.copy(),
            'dispersion': line_style_default.copy(),
            'minmax': line_style_default.copy(),
            'fill': {
                'color': 'gray',
                'alpha': 0.2
            }
        }
        self.stats_style['central']['color'] = 'black'
        self.stats_style['central']['linewidth'] = 1

        self.stats_style['dispersion']['color'] = 'black'
        self.stats_style['dispersion']['linewidth'] = 0
        self.stats_style['dispersion']['linestyle'] = 'dashed'

        self.stats_style['minmax']['color'] = 'black'
        self.stats_style['minmax']['linewidth'] = 0
        self.stats_style['minmax']['linestyle'] = 'dashed'

        self.stats_plot = {
            'central': 'avg',
            'dispersion': 'std',
            'factor': 1,
        }

    def start_plot(self):
        plt.close()
        plt.figure()

    def close_plot(self):
        plt.close()

    def save_plot(self, file_out):
        plt.savefig(file_out, dpi=300)

    def set_color(self,color):
        self.stats_style['central']['color'] = color
        self.stats_style['dispersion']['color'] = color
        self.stats_style['minmax']['color'] = color
        self.stats_style['fill']['color'] = color

    def set_iqr_as_stats_plot(self):
        self.stats_plot['central'] = 'median'
        self.stats_plot['dispersion'] = 'iqr'

    def set_std_as_stats_plot(self, factor):
        self.stats_plot['central'] = 'avg'
        self.stats_plot['dispersion'] = 'std'
        self.stats_plot['factor'] = factor

    def plot_data(self, ydata, style):
        h = plt.plot(self.xdata, ydata,
                     color=style['color'],
                     linestyle=style['linestyle'],
                     linewidth=style['linewidth'],
                     marker=style['marker'],
                     markersize=style['markersize'])
        return h

    def plot_multiple_data(self,ydata,style):
        if len(ydata.shape)==2 and ydata.shape[0]==len(self.xdata):
            xdata_multiple = np.zeros(ydata.shape)
            xdata_multiple = xdata_multiple + np.array(self.xdata)
            h = plt.plot(xdata_multiple, ydata,
                         color=style['color'],
                         linestyle=style['linestyle'],
                         linewidth=style['linewidth'],
                         marker=style['marker'],
                         markersize=style['markersize'])
            return h

    def plot_single_line(self, ydata, line_color, line_type, line_width, marker, marker_size):
        style = self.line_style_default
        if line_color is not None:
            style['color'] = line_color
        if line_type is not None:
            style['linestyle'] = line_type
        if line_width is not None:
            style['linewidth'] = line_width
        if marker is not None:
            style['marker'] = marker
        if marker_size is not None:
            style['markersize'] = marker_size

        h = self.plot_data(ydata, style)

        return h

    def plot_single_marker(self, xpoint, ypoint, marker, marker_size, color, edge_color, edge_width):
        h = plt.plot(xpoint, ypoint,
                 color=color,
                 linewidth=0,
                 marker=marker,
                 markersize=marker_size,
                 mec=edge_color,
                 mew=edge_width)
        return h

    def plot_single_bar_series(self, ydata, color, width, offset,linewidth):

        h = plt.bar(self.xdata + offset, ydata, width, color=color, linewidth=linewidth, edgecolor='k')
        return h

    def set_legend(self, str_legend):
        print(self.legend_options)
        plt.legend(str_legend, loc=self.legend_options['loc'], bbox_to_anchor=self.legend_options['bbox_to_anchor'],
                   framealpha=self.legend_options['framealpha'], ncol=self.legend_options['ncols'])

    def set_legend_h(self, handles, str_legend):
        plt.legend(handles, str_legend, loc=self.legend_options['loc'],
                   bbox_to_anchor=self.legend_options['bbox_to_anchor'], framealpha=self.legend_options['framealpha'],
                   ncol=self.legend_options['ncols'],markerscale=self.legend_options['markerscale'])

    def set_title(self, title):
        plt.title(title)
    def set_title_size(self, title, fontsize):
        plt.title(title,fontsize=fontsize)


    def set_xticks(self, xticks, xtickvalues, rotation, fontsize):
        if rotation is None:
            rotation = 0
        if rotation < 0 or rotation > 90:
            rotation = 0
        if xtickvalues is None:
            xtickvalues = xticks
        if fontsize is None:
            fontsize = 9
        plt.xticks(xticks, xtickvalues, rotation=rotation, fontsize=fontsize)

    def set_xticks_minor(self, xticks, xtickvalues, rotation, fontsize):
        if rotation is None:
            rotation = 0
        if rotation < 0 or rotation > 90:
            rotation = 0
        if xtickvalues is None:
            xtickvalues = xticks
        if fontsize is None:
            fontsize = 9
        plt.xticks(xticks, xtickvalues, minor=True, rotation=rotation, fontsize=fontsize)
        # plt.xticks([],minor = True)

    def set_yticks(self, yticks, ytickvalues, rotation, fontsize):
        if rotation is None:
            rotation = 0
        if rotation < 0 or rotation > 90:
            rotation = 0
        if ytickvalues is None:
            ytickvalues = yticks
        if fontsize is None:
            fontsize = 9
        if yticks is None:
            plt.yticks(rotation=rotation, fontsize=fontsize)
        else:
            plt.yticks(yticks, ytickvalues, rotation=rotation, fontsize=fontsize)

    def get_xticks(self):
        locs, labels = plt.xticks()
        return locs, labels

    def get_yticks(self):
        locs, labels = plt.yticks()
        return locs

    def set_xaxis_title(self, xaxis_title):
        plt.xlabel(xaxis_title, fontsize=14)

    def set_yaxis_title(self, yaxis_title):
        plt.ylabel(yaxis_title, fontsize=14)

    def set_xaxis_title_f(self, xaxis_title, fontsize):
        plt.xlabel(xaxis_title, fontsize=fontsize)

    def set_yaxis_title_f(self, yaxis_title, fontsize):
        plt.ylabel(yaxis_title, fontsize=fontsize)

    def save_fig(self, file_out):
        # plt.savefig(file_out, dpi=300)
        if file_out.endswith('.tif'):
            plt.savefig(file_out, dpi=300, bbox_inches='tight', pil_kwargs={"compression": "tiff_lzw"})
        else:
            plt.savefig(file_out, dpi=300, bbox_inches='tight')

    def set_equal_apect(self):
        plt.gca().set_aspect('equal', adjustable='box')

    def set_grid(self):
        # plt.grid(b=True, which='major', color='gray', linestyle='--')
        plt.grid(which='major', color='gray', linestyle='--', axis='both')

    def set_grid_horizontal(self):
        plt.grid(which='major',color='lightgray',linestyle='--',axis='y')

    def set_tigth_layout(self):
        plt.gcf().tight_layout()

    def set_y_range(self, ymin, ymax):
        plt.ylim(ymin, ymax)

    def get_y_range(self):
        ymin, ymax = plt.ylim()
        return ymin, ymax

    def plot_text(self, xpos, ypos, str):
        htext = plt.text(xpos, ypos, str)
        return htext

    def plot_multiple_spectra(self, wavelength, spectra, stats, wlmin, wlmax):
        imin, imax = self.get_imin_imax_from_wavelength(wavelength, wlmin, wlmax)

        # print(wlmin,wlmax,imin,imax)
        self.xdata = wavelength[imin:imax]
        if spectra is not None:
            spectra = spectra[:, imin:imax]
            # for ispectra in range(0,spectra.shape[0]):
            #     self.plot_data(spectra[ispectra,:],self.spectra_style)
            # print(len(self.xdata))
            #
            # print(spectra.transpose().shape)
            self.plot_data(spectra.transpose(), self.spectra_style)

        if stats is None:
            return imin, imax

        self.plot_data(stats['avg'][imin:imax], self.stats_style['central'])
        y1 = stats['avg'][imin:imax] - stats['std'][imin:imax]
        y2 = stats['avg'][imin:imax] + stats['std'][imin:imax]
        self.plot_data(y1, self.stats_style['dispersion'])
        self.plot_data(y2, self.stats_style['dispersion'])
        plt.fill_between(self.xdata, y1, y2, color=self.stats_style['fill']['color'],
                         alpha=self.stats_style['fill']['alpha'])
        self.plot_data(stats['spectra_min'][imin:imax], self.stats_style['minmax'])
        self.plot_data(stats['spectra_max'][imin:imax], self.stats_style['minmax'])

        ymin, ymax = self.get_ymin_ymax_from_stats(stats, imin, imax)
        self.set_y_range(ymin, ymax)

        return imin, imax

    def plot_stats(self, stats, imin, imax):

        if imin is None:
            imin = 0
        if imax is None:
            imax = len(stats['avg'])

        stat_central = self.stats_plot['central']  # avg or median
        stat_dispersion = self.stats_plot['dispersion']  # std o iqr
        # h = self.plot_data(stats['avg'][imin:imax], self.stats_style['central'])
        h = self.plot_data(stats[stat_central][imin:imax], self.stats_style['central'])

        if stat_dispersion == 'std':
            factor = self.stats_plot['factor']
            y1 = stats['avg'][imin:imax] - (factor * stats['std'][imin:imax])
            y2 = stats['avg'][imin:imax] + (factor * stats['std'][imin:imax])
        elif stat_dispersion == 'iqr':
            y1 = stats['p25'][imin:imax]
            y2 = stats['p75'][imin:imax]

        if self.stats_style['dispersion']['linewidth'] > 0:
            self.plot_data(y1, self.stats_style['dispersion'])
            self.plot_data(y2, self.stats_style['dispersion'])

        if self.stats_style['minmax']['linewidth'] > 0:
            self.plot_data(stats['spectra_min'][imin:imax], self.stats_style['minmax'])
            self.plot_data(stats['spectra_max'][imin:imax], self.stats_style['minmax'])

        if self.stats_style['fill']['color'] is not None:
            plt.fill_between(self.xdata, y1, y2, facecolor=self.stats_style['fill']['color'],alpha=self.stats_style['fill']['alpha'])

        return h

    def plot_iqr_basic(self, y1, y2, color):
        plt.fill_between(self.xdata, y1, y2, facecolor=color, alpha=0.5)

    def get_ymin_ymax_from_stats(self, stats, imin, imax):
        y1 = stats['avg'][imin:imax] - (2 * stats['std'][imin:imax])
        y2 = stats['avg'][imin:imax] + (2 * stats['std'][imin:imax])
        ymin = np.min(stats['spectra_min'])
        yminstd = np.min(y1)
        if yminstd > ymin:
            ymin = yminstd
        ymax = np.max(stats['spectra_max'])
        ymaxstd = np.max(y2)
        if ymaxstd < ymax:
            ymax = ymaxstd
        return ymin, ymax

    def get_imin_imax_from_wavelength(self, wavelength, wlmin, wlmax):
        imin = 0
        imax = len(wavelength) - 1
        if wlmin is not None and wavelength[imin] < wlmin < wavelength[imax]:
            imin = np.argmin(np.abs(wavelength - wlmin))
        if wlmax is not None and wavelength[imin] < wlmax < wavelength[imax]:
            imax = np.argmin(np.abs(wavelength - wlmax))
        imax = imax + 1
        return imin, imax
