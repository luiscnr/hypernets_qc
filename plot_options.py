import os.path
import MDBPlotDefaults as defaults


class PlotOptions:

    def __init__(self, options, config_file):

        if options is None and config_file is not None and os.path.exists(config_file):
            import configparser
            options = configparser.ConfigParser()
            options.read(config_file)

        self.options = options

        ##GLOBAL OPTIONS
        self.global_options = {}
        # self.mu_valid_variable = 'mu_valid'
        # self.format_image = 'png'
        # self.image_resolution = 300
        # self.output_path = None

        self.valid_stats = {
            'N': 0,
            'NMU': 0,
            'slope': 0.0,
            'intercept': 0.0,
            'PCC(r)': 0.0,
            'p_value': 0.0,
            'std_err': 0.0,
            'RMSD': 0.0,
            'RPD': 0.0,
            'APD': 0.0,
            'BIAS': 0.0,
            'DETER(r2)': 0.0,
            'SLOPE': 0.0,
            'OFFSET': 0.0,
            'XAVG': 0.0,
            'YAVG': 0.0,
            'CPRMSE': 0.0,
            'MAE': 0.0
        }

    def set_global_options(self):
        section = 'GLOBAL_OPTIONS'
        self.global_options = {}
        for goption in defaults.global_options:
            default = defaults.global_options[goption]['default']
            type = defaults.global_options[goption]['type']
            potential_values = None
            if 'values' in defaults.global_options[goption].keys():
                potential_values = defaults.global_options[goption]['values']
            self.global_options[goption] = self.get_value_param(section, goption, default, type, potential_values)
            # if type=='str' and 'values' in defaults.global_options[goption].keys():
            #     values = defaults.global_options[goption]['values']
            #     if not self.global_options[goption] in values:
            #         print(f'[ERROR] [{section}] {self.global_options[goption]} is not a valid  value for {goption}. Valid values: {values} ')

        # self.output_path = self.get_value_param(section, 'output_path', self.output_path, 'directory')
        # self.mu_valid_variable = self.get_value_param(section, 'mu_valid_variable', self.mu_valid_variable, 'str')

    def get_list_figures(self):
        sections = self.options.sections()
        list_figures = []
        for s in sections:
            apply = self.get_value_param(s, 'apply', False, 'boolean', None)
            if apply:
                list_figures.append(s)
        return list_figures

    def get_options(self, section):
        options_out = {'apply': self.get_value_param(section, 'apply', False, 'boolean', None)}
        if not options_out['apply']:
            return None
        options_out['type'] = self.get_value_param(section, 'type', None, 'str', defaults.type_list)
        if options_out['type'] is None:
            return None
        options_out['name'] = section
        if self.global_options['output_path'] is not None:
            name_default = options_out['name'] + '.' + self.global_options['fig_extension']
            file_out_default = os.path.join(self.global_options['output_path'], name_default)
        else:
            file_out_default = None
        options_out['file_out'] = self.get_value_param(section, 'file_out', file_out_default, 'str', None)

        # options_out['multiple_plot'] = self.get_value_param(section, 'multiple_plot', None, 'str')
        # if options_out['type'] == 'csvtable':
        #     options_out = self.get_options_csv(section,options_out)
        # if options_out['type'] == 'scatterplot':
        #     print(f'[INFO] Plot type: scatterplot')
        #     options_out = self.get_options_scatterplot(section, options_out)
        #
        # if options_out['type'] == 'spectraplot':
        #     print(f'[INFO] Plot type: spectraplot')
        #     options_out = self.get_options_spectraplot(section, options_out)

        if options_out['type'].startswith('statstable'):
            print(f'[INFO] Plot type: statstable')
            options_out = self.get_options_csv_statstable(section, options_out)
        else:
            print(f'[INFO] Plot type: {options_out["type"]}')
            options_out = self.get_options_impl(options_out['type'], section, options_out)

        # if options_out['type'] == 'histogram':
        #     print(f'[INFO] Plot type: histogram')
        #     options_out = self.get_options_histogram(section,options_out)
        #
        # if options_out['type'] == 'timeseries':
        #     print(f'[INFO] Plot type: timeseries')
        #     options_out = self.get_options_timeseries(section, options_out)
        #
        # if options_out['type'] == 'sequence':
        #     print(f'[INFO] Plot type: sequence')
        #     options_out = self.get_options_sequence(section, options_out)

        # if options_out['type'] == 'statswlplot':
        #     options_out = self.get_select_options(section, options_out)
        #     options_out = self.get_options_statswlplot(section, options_out)
        # if options_out['type'] == 'spectraplot':
        #     options_out['type_rrs'] = self.get_value_param(section, 'type_rrs', 'ins', 'str')
        #     if options_out['type_rrs'].startswith('flag'):
        #         options_out = self.get_group_options(section, options_out)
        #     options_out = self.get_options_spectraplot(section, options_out)
        #     options_out = self.get_select_options(section, options_out)
        # if options_out['type'] == 'flagplot':
        #      options_out = self.get_options_flag(section, options_out)
        if options_out is not None:
            for option in options_out:
                print(option, '->', options_out[option])

        return options_out

    def get_options_impl(self, type, section, options_out):
        doptions = None
        if type == 'scatterplot':
            doptions = defaults.get_options_scatterplots()
        if type == 'spectraplot':
            doptions = defaults.get_options_spectraplots()
        if type == 'timeseries':
            doptions = defaults.get_options_timeseries()
        if type == 'sequence':
            doptions = defaults.get_options_sequence()
        if type == 'histogram':
            doptions = defaults.get_options_histogram()
        if type == 'flagplot':
            doptions = defaults.get_options_flag_plot()
        if type == 'angleplot':
            doptions = defaults.get_options_angleplot()

        if doptions is None:
            print(f'[ERROR] Options for plot {type} are not implemented yet')
            return None
        for option in doptions:
            pvalues = None
            if 'values' in doptions[option]:
                pvalues = doptions[option]['values']
            options_out[option] = self.get_value_param(section, option, doptions[option]['default'],
                                                       doptions[option]['type'], pvalues)
        if type == 'scatterplot':
            options_out = self.get_options_scatterplot(section, options_out)

        return options_out

    def get_fix_time_axis(self, frequency, units, time_start, time_stop):

        from datetime import timedelta
        time_instants = []
        time_here = time_start
        if units == 'months':
            frequency = frequency * 30
        if units == 'years':
            frequency = frequency * 365
        while time_here <= time_stop:
            time_instants.append(time_here)
            if units == 'minutes':
                time_here = time_here + timedelta(minutes=frequency)
            elif units == 'hours':
                time_here = time_here + timedelta(hours=frequency)
            elif units == 'days' or units == 'months' or units == 'years':
                time_here = time_here + timedelta(days=frequency)
            else:
                break

        return time_instants

    def get_options_csv_statstable(self, section, options_out):
        options_out['xvar'] = self.get_value_param(section, 'xvar', 'mu_ins_rrs', 'str')
        options_out['yvar'] = self.get_value_param(section, 'yvar', 'mu_sat_rrs', 'str')
        options_out['params'] = self.get_value_param(section, 'params', self.valid_stats.keys(), 'strlist')
        options_out['log_scale'] = self.get_value_param(section, 'log_scale', False, 'boolean')
        options_out['use_rhow'] = self.get_value_param(section, 'use_rhow', False, 'boolean')
        options_out['flag'] = self.get_value_param(section, 'flag', 'GLOBAL', 'str')
        if self.output_path is not None:
            name_default = options_out['name'] + '.csv'
            file_out_default = os.path.join(self.output_path, name_default)
        options_out['file_out'] = self.get_value_param(section, 'file_out', file_out_default, 'str')
        return options_out

    def get_options_scatterplot(self, section, options_out):

        if options_out['type_scatterplot'] == 'chla' or options_out['type_scatterplot'] == 'kd':  ##CHANGE DEFAULTS
            if self.get_value(section, 'units') is None:
                options_out['units'] = defaults.units_default[options_out['type_scatterplot']]
            if self.get_value(section, 'scale_factor') is None:
                options_out['scale_factor'] = 1
            if self.get_value(section, 'log_scale') is None:
                options_out['log_scale'] = True
        if options_out['type_scatterplot'] == 'rrs' and options_out['use_rhow']:
            if self.get_value(section, 'units') is None:
                options_out['units'] = defaults.units_default['rhow']
            if self.get_value(section, 'scale_factor') is None:
                options_out['scale_factor'] = 1
        if options_out['xlabel'] is None:
            options_out['xlabel'] = defaults.get_label_scatterplot('x', options_out['type_scatterplot'],
                                                                   options_out['use_rhow'], options_out['scale_factor'])
        if options_out['ylabel'] is None:
            options_out['ylabel'] = defaults.get_label_scatterplot('y', options_out['type_scatterplot'],
                                                                   options_out['use_rhow'], options_out['scale_factor'])

        options_out['scale_factor_str'] = defaults.get_scale_factor_str(options_out['scale_factor'])

        if options_out['groupBy'] is not None:
            if self.get_value(section, 'edgecolor') is None:
                options_out['edgecolor'] = ['black']
            if self.get_value(section, 'linewidth') is None:
                options_out['linewidth'] = [0.25]
        density_plots = False
        if options_out['groupBy'] is None and options_out['apply_density']:
            if not options_out['selectByWavelength'] and not options_out['apply_wavelength_color']:
                density_plots = True
            if options_out['selectByWavelength']:
                density_plots = True
        if density_plots:
            if self.get_value(section, 'edgecolor') is None:
                options_out['edgecolor'] = ['blue']
            if self.get_value(section, 'linewidth') is None:
                options_out['linewidth'] = [0]

        if options_out['include_stats']:
            stat_list = options_out['stat_list']
            for stat in stat_list:
                if stat.upper() == 'WL':
                    continue
                elif stat.upper() == 'EQUATION':
                    type_regression = options_out['type_regression'].upper()

                    val_format_slope = self.get_value(section, f'SLOPE_{type_regression}_FORMAT')
                    if val_format_slope is None:
                        val_format_slope = defaults.valid_stats[f'SLOPE_{type_regression}']['format']
                    options_out[f'SLOPE_{type_regression}_FORMAT'] = val_format_slope.strip()

                    val_format_offset = self.get_value(section, f'OFFSET_{type_regression}_FORMAT')
                    if val_format_offset is None:
                        val_format_offset = defaults.valid_stats[f'OFFSET_{type_regression}']['format']
                    options_out[f'OFFSET_{type_regression}_FORMAT'] = val_format_offset.strip()
                else:
                    val_format = self.get_value(section, f'{stat}_FORMAT')
                    val_nameplot = self.get_value(section, f'{stat}_NAMEPLOT')
                    if val_format is None:
                        val_format = defaults.valid_stats[stat.upper()]['format']
                    if val_nameplot is None:
                        val_nameplot = stat.upper()
                        if 'name_plot' in defaults.valid_stats[stat.upper()]:
                            val_nameplot = defaults.valid_stats[stat.upper()]['name_plot']

                    options_out[f'{stat.upper()}_FORMAT'] = val_format.strip()
                    options_out[f'{stat.upper()}_NAMEPLOT'] = val_nameplot.strip()

        if options_out['log_scale']:
            if options_out['min_xy'] is None:
                options_out['min_xy'] = 0.1
                if options_out['type_scatterplot'] == 'kd':
                    options_out['min_xy'] = 0.01
            if options_out['max_xy'] is None:
                options_out['max_xy'] = 100
                if options_out['type_scatterplot'] == 'kd':
                    options_out['max_xy'] = 10

        return options_out

    def get_options_scatterplot_deprecated(self, section, options_out):

        options_out['type_scatterplot'] = self.get_value_param(section, 'type_scatterplot', 'rrs', 'str')

        options_out['xvar'] = self.get_value_param(section, 'xvar', 'mu_ins_rrs', 'str')
        options_out['yvar'] = self.get_value_param(section, 'yvar', 'mu_sat_rrs', 'str')

        options_out['legend'] = self.get_value_param(section, 'legend', True, 'boolean')
        options_out['legend_values'] = self.get_value_param(section, 'legend_values', None, 'strlist')
        options_out['include_stats'] = self.get_value_param(section, 'include_stats', False, 'boolean')

        options_out['regression_line_groups'] = self.get_value_param(section, 'regression_line_groups', False,
                                                                     'boolean')
        options_out['apply_density'] = self.get_value_param(section, 'apply_density', True, 'boolean')
        options_out['title'] = self.get_value_param(section, 'title', None, 'str')
        # print(self.global_options)
        if self.global_options['output_path'] is not None:
            name_default = options_out['name'] + '.' + self.global_options['fig_extension']
            file_out_default = os.path.join(self.global_options['output_path'], name_default)
            # name_default = options_out['name'] + '.' + self.format_image
            # file_out_default = os.path.join(self.output_path, name_default)
        options_out['file_out'] = self.get_value_param(section, 'file_out', file_out_default, 'str')
        options_out['log_scale'] = self.get_value_param(section, 'log_scale', False, 'boolean')
        options_out['use_rhow'] = self.get_value_param(section, 'use_rhow', False, 'boolean')
        options_out['min_xy'] = self.get_value_param(section, 'min_xy', None, 'float')
        options_out['max_xy'] = self.get_value_param(section, 'max_xy', None, 'float')
        options_out['ticks'] = self.get_value_param(section, 'ticks', None, 'floatlist')
        options_out['fontsizeaxis'] = self.get_value_param(section, 'fontsizeaxis', 12, 'float')
        options_out['fontsizelabels'] = self.get_value_param(section, 'fontsizelabels', 12, 'float')
        options_out['fontsizetitle'] = self.get_value_param(section, 'fontsizetitle', 12, 'float')
        options_out['fontsizestats'] = self.get_value_param(section, 'fontsizestats', 12, 'float')
        sfdefault = None
        unitsdefault = None
        if options_out['type_scatterplot'] == 'rrs':
            unitsdefault = r'sr$^-$$^1$'
            sfdefault = 1000
        if options_out['type_scatterplot'] == 'chla':
            unitsdefault = r'mg m$^-$$^3$'
            sfdefault = 1
            options_out['log_scale'] = True
        if options_out['type_scatterplot'] == 'kd':
            unitsdefault = r'm$^-$$^1$'
            sfdefault = 1
            options_out['log_scale'] = True
        xlabeldefault = defaults.xlabel_default
        ylabeldefault = defaults.ylabel_default
        options_out['scale_factor'] = self.get_value_param(section, 'scale_factor', sfdefault, 'float')
        options_out['xlabel'] = self.get_value_param(section, 'xlabel', xlabeldefault, 'str')
        options_out['ylabel'] = self.get_value_param(section, 'ylabel', ylabeldefault, 'str')
        options_out['units'] = self.get_value_param(section, 'units', unitsdefault, 'str')
        options_out['identity_line'] = self.get_value_param(section, 'identity_line', True, 'boolean')
        options_out['regression_line'] = self.get_value_param(section, 'regression_line', True, 'boolean')
        # marker, markersize, color, edgecolor, linewidth
        # o, 25, 'black', None, None
        options_out['marker'] = self.get_value_param(section, 'marker', ['o'], 'strlist')
        options_out['markersize'] = self.get_value_param(section, 'markersize', [25], 'intlist')
        options_out['color'] = self.get_value_param(section, 'color', ['black'], 'strlist')

        edgeColorDefault = None
        lineWidthDefault = None
        # if options_out['groupBy'] is not None:
        # if options_out['groupType'] == 'rrs':
        #     edgeColorDefault = 'gray'
        #     lineWidthDefault = 1.5
        # if options_out['groupType'] == 'flag':
        #     edgeColorDefault = 'black'
        #     lineWidthDefault = 0.25
        options_out['edgecolor'] = self.get_value_param(section, 'edgecolor', [edgeColorDefault], 'strlist')
        options_out['linewidth'] = self.get_value_param(section, 'linewidth', [lineWidthDefault], 'floatlist')

        options_out['xfigsize'] = self.get_value_param(section, 'xfigsize', 7, 'float')
        options_out['yfigsize'] = self.get_value_param(section, 'yfigsize', 7, 'float')
        options_out['widthspace'] = self.get_value_param(section, 'widthspace', 0.1, 'float')
        options_out['heightspace'] = self.get_value_param(section, 'heightspace', 0.1, 'float')
        options_out['stat_list'] = self.get_value_param(section, 'stat_list', None, 'strlist')
        options_out['stats_xpos'] = self.get_value_param(section, 'stats_xpos', 0.05, 'float')
        options_out['stats_ypos'] = self.get_value_param(section, 'stats_ypos', 0.70, 'float')
        options_out['individual_axis'] = self.get_value_param(section, 'individual_axis', False, 'boolean')
        return options_out

    def get_value(self, section, key):
        value = None
        if self.options.has_option(section, key):
            value = self.options[section][key]
        return value

    def get_value_param(self, section, key, default, type, potential_values):
        value = self.get_value(section, key)
        if value is None:
            return default
        if type == 'str':
            if potential_values is None:
                return value.strip()
            else:
                if value.strip().lower() in potential_values:
                    return value
                else:
                    print(
                        f'[ERROR] [{section}] {value} is not a valid  value for {key}. Valid values: {potential_values} ')
                    return default
        if type == 'file':
            if not os.path.exists(value.strip()):
                return default
            else:
                return value.strip()
        if type == 'directory':
            directory = value.strip()
            if not os.path.isdir(directory):
                try:
                    os.mkdir(directory)
                    return directory
                except:
                    return default
            else:
                return directory
        if type == 'int':
            return int(value)
        if type == 'float':
            return float(value)
        if type == 'boolean':
            if value == '1' or value.upper() == 'TRUE':
                return True
            elif value == '0' or value.upper() == 'FALSE':
                return False
            else:
                return True
        if type == 'rrslist':
            list_str = value.split(',')
            list = []
            for vals in list_str:
                vals = vals.strip().replace('.', '_')
                list.append(f'RRS{vals}')
            return list
        if type == 'strlist':
            list_str = value.split(',')
            list = []
            for vals in list_str:
                list.append(vals.strip())
            return list
        if type == 'floatlist':
            list_str = value.split(',')
            list = []
            for vals in list_str:
                vals = vals.strip()
                list.append(float(vals))
            return list
        if type == 'floattuple':
            list_str = value.split(',')
            list = []
            for vals in list_str:
                vals = vals.strip()
                list.append(float(vals))
            return tuple(list)
        if type == 'intlist':
            list_str = value.split(',')
            list = []
            for vals in list_str:
                vals = vals.strip()
                list.append(int(vals))
            return list

        if type == 'linestyle':
            list_str = value.split(',')
            if len(list_str) != 5:
                print(
                    f'[WARNING] {section}-{key} is not valid, it should be defined as a line style with 5 comma-separated values:')
                print(f'[WARNING] line_color,marker,marker_size,line_style,line_size')
                print(f'[WARNING] Using default line style for {section}-{key}')
                return default
            try:
                marker_size = float(list_str[2].strip())
                line_size = float(list_str[4].strip())
                marker = list_str[1].strip()
                if marker.lower() == 'none':
                    marker = None
                style = {
                    'color': list_str[0].strip(),
                    'marker': marker,
                    'markersize': marker_size,
                    'linestyle': list_str[3].strip(),
                    'linewidth': line_size
                }
                return style

            except:
                print(f'[WARNING] {section}-{key} is not valid line style, using default style')
                return default

        if type == 'fillstyle':
            list_str = value.split(',')
            if len(list_str) != 2:
                print(
                    f'[WARNING] {section}-{key} is not valid, it should be defined as a fill style with 2 comma-separated values:')
                print(f'[WARNING] color, alpha')
                print(f'[WARNING] Using default fill style for {section}-{key}')
                return default
            try:
                alpha = float(list_str[1].strip())
                style = {
                    'color': list_str[0].strip(),
                    'alpha': alpha
                }
                return style
            except:
                print(f'[WARNING] {section}-{key} is not valid fill style, using default style')
                return default

        if type == 'date':
            from datetime import datetime as dt
            try:
                date = dt.strptime(value.strip(), '%Y-%m-%d')
                return date
            except:
                return default
        if type == 'time':
            from datetime import datetime as dt
            val = value.strip()
            try:
                val_check = dt.now().strftime('%Y-%m-%d')
                val_check = f'{val_check}T{val}'
                dt.strptime(val_check, '%Y-%m-%dT%H:%M')
                return val
            except:
                return default
