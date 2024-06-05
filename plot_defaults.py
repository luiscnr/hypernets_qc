import math
import matplotlib as mpl

color_dict = dict({ \
    '400.00': 'LightBlue', \
    '412.50': 'DeepSkyBlue', \
    '442.50': 'DodgerBlue', \
    '490.00': 'Blue', \
    '510.00': 'ForestGreen', \
    '560.00': 'Green', \
    '620.00': 'LightCoral', \
    '665.00': 'Red', \
    '673.75': 'Crimson', \
    '681.25': 'FireBrick', \
    '708.75': 'LightGray', \
    '753.75': 'Silver', \
    '778.75': 'DarkGray', \
    '865.00': 'Gray', \
    '885.00': 'DimGray', \
    '1020.50': 'Pink'})

# xlabel_default = {
#     'rrs': r'In situ R$_r$$_s$ (10$^-$$^3$ sr$^-$$^1$)',
#     'chla': r'In situ chl-a (mg m$^-$$^3$)',
#     'kd': r'In situ Kd (m$^-$$^1$)',
#     'rhow': r'In situ ρ$_w$'
# }
# ylabel_default = {
#     'rrs': r'Satellite R$_r$$_s$ (10$^-$$^3$ sr$^-$$^1$)',
#     'chla': r'Satellite chl-a (mg m$^-$$^3$)',
#     'kd': r'Satellite Kd (m$^-$$^1$)',
#     'rhow': r'Satellite ρ$_w$'
# }

ylabel_rrs = r'R$_r$$_s$ (sr$^-$$^1$)'
ylabel_rrs_scaled = r'R$_r$$_s$ (10$^-$$^3$ sr$^-$$^1$)'

label_insitu_default = 'In situ Rrs ()'
xlabel_wl_default = 'Wavelength (nm)'

colors_default = ['Blue', 'Red', 'Green', 'm', 'Cyan', 'Orange', 'Yellow']
fontsizetitle_default = 12
fontsizeaxis_default = 12
fontsizelabels_default = 12
fontsizestats_default = 12
units_default = {
    'rrs': r'sr$^-$$^1$',
    'rhow': '',
    'chla': r'mg m$^-$$^3$',
    'kd': r'm$^-$$^1$'
}

marker_default = 'o'
marker_size_default = 25
marker_color_default = 'black'
edge_color_default = 'gray'
line_width_default = 0.5
regression_line_style_default = {
    'color': 'black',
    'marker': None,
    'markersize': 0,
    'linestyle': 'solid',
    'linewidth': 1
}
identity_line_style_default = {
    'color': 'black',
    'marker': None,
    'markersize': 0,
    'linestyle': 'dashed',
    'linewidth': 0.75
}

line_style_default = {
    'color': 'black',
    'marker': None,
    'markersize': 0,
    'linestyle': 'solid',
    'linewidth': 1
}
all_line_style_default = {
    'color': 'gray',
    'marker': None,
    'markersize': 0,
    'linestyle': 'solid',
    'linewidth': 1
}
valid_line_style_default = {
    'color': 'green',
    'marker': None,
    'markersize': 0,
    'linestyle': 'solid',
    'linewidth': 1
}
invalid_line_style_default = {
    'color': 'red',
    'marker': None,
    'markersize': 0,
    'linestyle': 'solid',
    'linewidth': 1
}
selected_line_style_default = {
    'color': 'blue',
    'marker': None,
    'markersize': 0,
    'linestyle': 'solid',
    'linewidth': 1
}
central_line_style_default = {
    'color': 'black',
    'marker': 'o',
    'markersize': 5,
    'linestyle': 'solid',
    'linewidth': 1.5
}
dispersion_line_style_default = {
    'color': 'black',
    'marker': None,
    'markersize': 0,
    'linestyle': 'dashed',
    'linewidth': 0
}
fill_style_default ={
    'color': 'gray',
    'alpha': 0.5
}
insitu_line_style_default = {
    'color': 'red',
    'marker': '.',
    'markersize': 10,
    'linestyle': 'solid',
    'linewidth': 1
}
sat_line_style_default = {
    'color': 'blue',
    'marker': '.',
    'markersize': 10,
    'linestyle': 'solid',
    'linewidth': 1
}
insitu_central_style_default = {
    'color': 'red',
    'marker': 'o',
    'markersize': 5,
    'linestyle': 'solid',
    'linewidth': 1
}
sat_central_style_default = {
    'color': 'blue',
    'marker': 'o',
    'markersize': 5,
    'linestyle': 'solid',
    'linewidth': 1
}
insitu_dispersion_style_default = {
    'color': 'red',
    'marker': 'o',
    'markersize': 0,
    'linestyle': 'dashed',
    'linewidth': 0
}
sat_dispersion_style_default = {
    'color': 'red',
    'marker': 'o',
    'markersize': 0,
    'linestyle': 'dashed',
    'linewidth': 0
}
insitu_fill_style_default ={
    'color': 'red',
    'alpha': 0.5
}
sat_fill_style_default ={
    'color': 'blue',
    'alpha': 0.5
}
global_options = {
    'output_path': {
        'default': None,
        'type': 'directory'
    },
    'mu_valid_variable': {
        'default': 'mu_valid',
        'type': 'str',
        'values': ['mu_valid', 'mu_valid_common']
    },
    'fig_extension': {
        'default': 'tif',
        'type': 'str',
        'values': ['tif', 'jpg', 'png']
    },
    'fig_resolution': {
        'default': 300,
        'type': 'float'
    }
}

type_list = ['scatterplot','statswlplot','spectraplot','multipleplot','flagplot','histogram','timeseries','sequence','angleplot']

valid_stats = {
    'N':{
        'name': 'N',
        'desc': 'Number of data points',
        'format': 'i'
    },
    'NMATCH-UPS': {
        'name': 'NMU',
        'name_plot': 'NMu',
        'desc': 'Number of match-ups',
        'format': 'i'
    },
    'NGROUP':{
        'name': 'NGROUP',
        'name_plot': 'NGroup',
        'desc': 'Number of data points by group',
        'format': 'i'
    },
    'SLOPE_I':{
        'name': 'slope_I',
        'desc': 'Slope Regression Type I',
        'format': 'f2'
    },
    'SLOPE_II':{
        'name': 'slope_II',
        'desc': 'Slope Regression Type II',
        'format': 'f2'
    },
    'OFFSET_I':{
        'name': 'intercept_I',
        'desc': 'Intercept Regression Type I',
        'format': 'f2'
    },
    'OFFSET_II':{
        'name': 'intercept_II',
        'desc': 'Intercept Regression Type II',
        'format': 'f2'
    },
    'STD_ERR_I':{
        'name': 'std_err_I',
        'desc': 'Standard Error Regression Type I',
        'format': 'f2'
    },
    'STD_SLOPE_II':{
        'name': 'std_slope_II',
        'desc': 'Standard Error Slope Regression Type II',
        'format': 'f2'
    },
    'STD_OFFSET_II':{
        'name': 'std_intercept_II',
        'desc': 'Standard Error Intercept Regression Type II',
        'format': 'f2'
    },
    'R':{
        'name': 'PCC(r)',
        'name_plot': 'r',
        'desc': 'Pearson coefficient',
        'format' : 'f2'
    },
    'P_VALUE':{
        'name': 'p_value',
        'name_plot': 'p-value',
        'desc': 'Correlaton p-value',
        'format': 'f2'
    },
    'R2':{
        'name': 'DETER(r2)',
        'name_plot': f'R$^2$',
        'desc': 'Determination coefficient',
        'format': 'f2'
    },
    'RMSD':{
        'name': 'RMSD',
        'desc': 'Root Mean Square Deviation',
        'format': 'e1+units'
    },
    'BIAS':{
        'name': 'BIAS',
        'name_plot': 'bias',
        'desc': 'Bias value',
        'format': 'e1+units'
    },
    'APD':{
        'name': 'APD',
        'desc': 'Absolute percent difference',
        'format': 'i'
    },
    'RPD':{
        'name': 'RPD',
        'desc': 'Relative percent difference',
        'format': 'i'
    },
    'XAVG':{
        'name': 'XAVG',
        'desc': 'X average value',
        'format': 'f3+units'
    },
    'YAVG':{
        'name': 'YAVG',
        'desc': 'Y average value',
        'format': 'f3+units'
    },
    'CRMSE':{
        'name': 'CRMSE',
        'desc': 'Centered Root Mean Square Deviation',
        'format': 'f3+units'
    },
    'MAE':{
        'name': 'MAE',
        'desc': 'Mean absolute error',
        'format': 'f3+units'
    }
}

options_legend = {
    'legend':{
        'default': True,
        'type': 'boolean'
    },
    'legend_values':{
        'default': None,
        'type': 'strlist'
    }
}

options_stats = {
    'include_stats':{
        'default': False,
        'type':'boolean'
    },
    'stat_list':{
        'default': ['N'],
        'type': 'strlist',
        'values': [x.lower() for x in valid_stats.keys()]
    },
    'stats_xpos':{
        'default': 0.05,
        'type': 'float'
    },
    'stats_ypos':{
        'default': 0.70,
        'type': 'float'
    },
    'fontsizestats':{
        'default': fontsizestats_default,
        'type': 'float'
    }

}

options_title = {
    'title':{
        'default': None,
        'type': 'str'
    },
    'fontsizetitle':{
        'default': fontsizetitle_default,
        'type': 'float'
    }
}

options_size = {
    'fontsizeaxis':{
        'default': fontsizeaxis_default,
        'type': 'float'
    },
    'fontsizelabels':{
        'default': fontsizelabels_default,
        'type': 'float'
    }
}

options_multiple_plot = {
    'multiple_plot':{
        'default': None,
        'type': 'str'
    },
    'xfigsize':{
        'default': 7,
        'type': 'float'
    },
    'yfigsize':{
        'default': 7,
        'type': 'float'
    },
    'widthspace':{
        'default': 0.1,
        'type': 'float'
    },
    'heightspace':{
        'default': 0.1,
        'type': 'float'
    }

}

options_group = {
    'groupBy':{
        'default':None,
        'type': 'str'
    },
    'groupValues':{
        'default':None,
        'type': 'strlist'
    },
    'groupType':{
        'default':'flag',
        'type': 'str',
        'values': ['flag','float']
    }
}

options_select = {
    'selectBy':{
        'default':None,
        'type': 'str'
    },
    'selectValues':{
        'default':None,
        'type': 'strlist'
    },
    'selectType':{
        'default':'flag',
        'type': 'str',
        'values': ['flag','float']
    },
    'selectByWavelength':{
        'default': False,
        'type': 'boolean'
    },
    'wlvalues': {
        'default': None,
        'type': 'floatlist'
    }
}

options_time = {
    'start_date': {
        'default': None,
        'type': 'date'
    },
    'end_date': {
        'default': None,
        'type': 'date'
    },
    'start_time':{
        'default': None,
        'type': 'time'
    },
    'end_time':{
        'default': None,
        'type': 'time'
    }
}

options_scatterplots = {
    'type_scatterplot':{
        'default':'rrs',
        'type':'str',
        'values': ['rrs','chla','kd','general']
    },
    'xvar':{
        'default': 'mu_ins_rrs',
        'type': 'str'
    },
    'yvar':{
        'default': 'mu_sat_rrs',
        'type': 'str'
    },
    'scale_factor':{
        'default': 1000,
        'type': 'float'
    },
    'units':{
        'default': units_default['rrs'],
        'type': 'str'
    },
    'xlabel':{
        'default': None,
        'type': 'str'
    },
    'ylabel':{
        'default': None,
        'type': 'str'
    },
    'log_scale':{
        'default': False,
        'type': 'boolean'
    },
    'use_rhow':{
        'default': False,
        'type': 'boolean'
    },
    'min_xy':{
        'default': None,
        'type': 'float'
    },
    'max_xy':{
        'default': None,
        'type': 'float'
    },
    'ticks':{
        'default':None,
        'type': 'floatlist'
    },
    'regression_line':{
        'default': True,
        'type': 'boolean'
    },
    'regression_line_style':{
        'default': regression_line_style_default,
        'type': 'linestyle'
    },
    'identity_line':{
        'default': True,
        'type': 'boolean'
    },
    'identity_line_style':{
        'default': identity_line_style_default,
        'type': 'linestyle'
    },
    'apply_wavelength_color':{
        'default': True,
        'type': 'boolean'
    },
    'apply_density':{
        'default': True,
        'type': 'boolean'
    },
    'regression_line_groups':{
        'default': False,
        'type': 'boolean'
    },
    'type_regression':{
        'default': 'II',
        'type': 'str',
        'values': ['i','ii']
    },
    'individual_axis':{
        'default': False,
        'type': 'boolean'
    },
    'marker':{
        'default': [marker_default],
        'type': 'strlist'
    },
    'markersize':{
        'default': [marker_size_default],
        'type': 'intlist'
    },
    'color':{
        'default': [marker_color_default],
        'type': 'strlist'
    },
    'edgecolor':{
        'default': [edge_color_default],
        'type': 'strlist'
    },
    'linewidth':{
        'default': [line_width_default],
        'type': 'floatlist'
    },
    'apply_wavelength_color ':{
        'default': True,
        'type': 'boolean'
    }

}


options_spectraplots = {
    'type_rrs':{
        'default':'comparison_sat_insitu',
        'type':'str',
        'values': ['ins','sat','mu_ins','mu_comparison','comparison_sat_insitu','comparison_sat','user_defined']
    },
    'wl_variable':{
        'default': None,
        'type': 'str'
    },
    'y_variable':{
        'default': None,
        'type': 'str'
    },
    'wl_min':{
        'default': None,
        'type': 'float'
    },
    'wl_max':{
        'default': None,
        'type': 'float'
    },
    'y_min':{
        'default': None,
        'type': 'float'
    },
    'y_max':{
        'default': None,
        'type': 'float'
    },'scale_factor':{
        'default': 1000,
        'type': 'float'
    },'use_rhow':{
        'default': False,
        'type': 'boolean'
    },
    'xlabel':{
        'default': xlabel_wl_default,
        'type': 'str'
    },
    'ylabel':{
        'default': None,
        'type': 'str'
    },
    'title':{
        'default': None,
        'type': 'str'
    },
    'plot_spectra':{
        'default': ['valid'],
        'type': 'strlist',
        'values': ['none','all','valid','invalid','selected']
    },
    'plot_stats':{
        'default': True,
        'type': 'boolean'
    },
    'stat_plot_method':{
        'default': 'iqr',
        'type': 'str'
    },
    'all_line_style':{
        'default': all_line_style_default,
        'type': 'linestyle'
    },
    'linestyle':{
        'default': line_style_default,
        'type': 'linestyle'
    },
    'color':{
        'default': [line_style_default['color']],
        'type': 'strlist'
    },
    'marker':{
        'default': [line_style_default['marker']],
        'type': 'strlist'
    },
    'markersize':{
        'default': [line_style_default['markersize']],
        'type': 'floatlist'
    },
    'linestyle':{
        'default': [line_style_default['linestyle']],
        'type': 'strlist'
    },
    'linewidth': {
        'default': [line_style_default['linewidth']],
        'type': 'strlist'
    },
    'valid_line_style':{
        'default': valid_line_style_default,
        'type': 'linestyle'
    },
    'invalid_line_style':{
        'default': invalid_line_style_default,
        'type': 'linestyle'
    },
    'selected_line_style':{
        'default': selected_line_style_default,
        'type': 'linestyle'
    },
    'central_line_style':{
        'default': central_line_style_default,
        'type': 'linestyle'
    },
    'dispersion_line_style':{
        'default': dispersion_line_style_default,
        'type': 'linestyle'
    },
    'fill_style':{
        'default': fill_style_default,
        'type': 'fillstyle'
    },
    'mu_range': {
        'default': None,
        'type': 'intrange'
    },
    'mu_list': {
        'default': None,
        'type': 'intlist'
    },
    'insitu_line_style':{
        'default': insitu_line_style_default,
        'type': 'linestyle'
    },
    'sat_line_style':{
        'default': sat_line_style_default,
        'type': 'linestyle'
    },
    'insitu_central_style':{
        'default': insitu_central_style_default,
        'type': 'linestyle'
    },
    'sat_central_style':{
        'default': sat_central_style_default,
        'type': 'linestyle'
    },
    'insitu_dispersion_style':{
        'default': insitu_dispersion_style_default,
        'type': 'linestyle'
    },
    'sat_dispersion_style':{
        'default': sat_dispersion_style_default,
        'type': 'linestyle'
    },
    'insitu_fill_style':{
        'default': insitu_fill_style_default,
        'type': 'fillstyle'
    },
    'sat_fill_style':{
        'default': sat_fill_style_default,
        'type': 'fillstyle'
    }
}

options_histogram = {
    'hvar': {
        'default': None,
        'type': 'str'
    },
    'type_histo': {
        'default': 'int',
        'type': 'str',
        'values': ['int','float']
    },
    'int_values': {
        'default': None,
        'type': 'intlist'
    },
    'hticks':{
        'default': None,
        'type': 'strlist'
    },
    'xlabel':{
        'default': None,
        'type': 'str'
    },
    'ylabel':{
        'default': 'NValues',
        'type': 'str'
    }
}

options_timeseries = {
    'time_var':{
        'default': None,
        'type': 'str'
    },
    'avg_var':{
        'default': None,
        'type': 'strlist'
    },
    'dispersion_min_var':{
        'default': None,
        'type': 'strlist'
    },
    'dispersion_max_var':{
        'default': None,
        'type': 'strlist'
    },'xlabel':{
        'default': None,
        'type': 'str'
    },
    'ylabel':{
        'default': None,
        'type': 'str'
    },
    'y_min':{
        'default': None,
        'type': 'float'
    },
    'y_max':{
        'default': None,
        'type': 'float'
    },
    'log_scale':{
        'default': False,
        'type': 'boolean'
    },
    'type_time_axis':{
        'default': 'variable',
        'type': 'str',
        'values': ['variable','fix']
    },
    'method_fix_axis':{
        'default': 'all',
        'type': 'str',
        'values': ['all','nearest']
    },
    'xticks_range':{
        'default': 1,
        'type': 'int'
    },
    'xticks_labels_range':{
        'default': 1,
        'type': 'int'
    },
    'color':{
        'default': [marker_color_default],
        'type': 'strlist'
    },
    'wlref':{
        'default': -999.0,
        'type': 'float'
    }
}

options_sequences = {
    'start_time':{
        'default': None,
        'type': 'str'
    },
    'end_time':{
        'default': None,
        'type': 'str'
    },
    'start_date':{
        'default': None,
        'type': 'str'
    },
    'end_date':{
        'default': None,
        'type': 'str'
    },
    'frequency_units':{
        'default': 'minutes',
        'type': 'str',
        'values': ['minutes','hours','days','months','years']
    },
    'frequency':{
        'default': -999.0,
        'type': 'float',
    }
}

options_flagplot = {
    'type_flagplot':{
        'default': 'comparison',
        'type': 'str',
        'values': ['comparison']
    }
}

options_angleplot = {
    'angle_var':{
        'default': None,
        'type': 'str'
    },
    'avg_var':{
        'default': None,
        'type': 'str'
    },
    'color':{
        'default': [marker_color_default],
        'type': 'strlist'
    },
    'wlref':{
        'default': -999.0,
        'type': 'float'
    },
    'scale':{
        'default': 'log',
        'type': 'str',
        'values': ['log','linear','symlog','logit']
    },
    'rlim':{
        'default': (1e-5,1),
        'type': 'floattuple'
    },
    'theta_zero_location':{
        'default': 'N',
        'type': 'str',
        'values': ['N', 'NW', 'W', 'SW', 'S', 'SE', 'E', 'NE']
    },
    'theta_direction':{
        'default': -1,
        'type': 'int',
        'values': [-1,1]
    },
    'theta_min':{
        'default': 0,
        'type': 'float'
    },
    'theta_max':{
        'default': 360,
        'type': 'float'
    },
    'xticks':{
        'default': [0,45,90,135,180,225,270,315],
        'type': 'floatlist'
    },
    'rticks':{
        'default': [1e-4, 1e-3, 1e-2, 1e-1, 1],
        'type': 'floatlist'
    },
    'rlabel_position':{
        'default': None,
        'type': 'str'
    },
    'label_size':{
        'default': 10,
        'type': 'float'
    },
    'point_size':{
        'default': [15],
        'type': 'floatlist'
    },
    'point_marker':{
        'default': ['o'],
        'type': 'strlist'
    },
    'min_data':{
        'default': 1.0e-5,
        'type': 'float'
    }


}

def get_options_spectraplots():
    options = options_spectraplots
    for op in options_legend:
        options[op] = options_legend[op]
    for op in options_multiple_plot:
        options[op] = options_multiple_plot[op]
    for op in options_select:
        options[op] = options_select[op]
    for op in options_group:
        options[op] = options_group[op]
    for op in options_time:
        options[op] = options_time[op]
    return options

def get_options_scatterplots():
    options = options_scatterplots
    for op in options_legend:
        options[op] = options_legend[op]
    for op in options_stats:
        options[op] = options_stats[op]
    for op in options_title:
        options[op] = options_title[op]
    for op in options_size:
        options[op] = options_size[op]
    for op in options_multiple_plot:
        options[op] = options_multiple_plot[op]
    for op in options_group:
        options[op] = options_group[op]
    for op in options_select:
        options[op] = options_select[op]
    return options

def get_options_histogram():
    return options_histogram

def get_options_timeseries():
    options = options_timeseries

    for op in options_sequences:
        options[op] = options_sequences[op]
    for op in options_legend:
        options[op] = options_legend[op]
    for op in options_title:
        options[op] = options_title[op]
    for op in options_group:
        options[op] = options_group[op]

    return options

def get_options_angleplot():
    options = options_angleplot

    for op in options_legend:
        options[op] = options_legend[op]
    for op in options_title:
        options[op] = options_title[op]
    for op in options_group:
        options[op] = options_group[op]

    return options


def get_options_sequence():
    options =  options_sequences
    for op in options_title:
        options[op] = options_title[op]
    return options

def get_options_flag_plot():
    options =  options_flagplot
    for op in options_title:
        options[op] = options_title[op]
    return options
def get_scale_factor_str(scale_factor):
    import numpy as np
    scale_factor_str = ''
    n = int(np.log10(scale_factor))
    if n > 1:
        scale_factor_str = f'10$^-$$^{n}$'
    return scale_factor_str

#type_label: 'x' or 'y'
def get_label_scatterplot(type_label,type_rrs,use_rhow,scale_factor):
    if type_label=='x':
        ini = 'In situ'
    elif type_label=='y':
        ini = 'Satellite'
    if type_rrs=='rrs':
        if use_rhow:
            quantity = r'ρ$_w$'
            units = units_default['rhow']
        else:
            quantity = r'R$_r$$_s$'
            units = units_default['rrs']
    elif type_rrs=='chla':
        quantity = 'chl-a'
        units = units_default['chla']
    elif type_rrs=='kd':
        quantity = 'Kd'
        units = units_default['kd']

    scale_factor_str = get_scale_factor_str(scale_factor)

    if len(units)>0:
        if len(scale_factor_str)>0:
            label = f'{ini} {quantity} ({scale_factor_str} {units})'
        else:
            label = f'{ini} {quantity} ({units})'
    else:
        if len(scale_factor_str) > 0:
            label = f'{ini} {quantity} ({scale_factor_str}'
        else:
            label = f'{ini} {quantity}'

    return label





def get_color_default(value, min, max):
    cm = mpl.colormaps['jet']
    return cm((value - min) / (max - min))

def get_color_wavelength(wlvalue):
    dif_ref = 10000
    color_out = None
    for wlp in color_dict:
        wlpv = float(wlp)
        dif = abs(wlpv - wlvalue)
        if dif < dif_ref:
            dif_ref = dif
            color_out = color_dict[wlp]
    print(wlvalue,'->',color_out)
    return color_out

def get_color_flag(flagvalue):
    index = int(math.log2(flagvalue))
    return colors_default[index]

def get_color_list(n):
    if n <= len(colors_default):
        return colors_default[0:n]
    else:
        return colors_default[0]
