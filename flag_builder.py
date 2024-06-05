import os, math
import numpy as np
from netCDF4 import Dataset


class FlagBuilder:

    def __init__(self, path_mdb_file, options):
        self.nc_file = path_mdb_file
        self.VALID = os.path.isfile(path_mdb_file)
        self.options = options
        if not self.VALID:
            print(f'[ERROR] {path_mdb_file} is not a valid  file')

        self.flag_options = {
            'type': {'type_param': 'str', 'list_values': ['virtual_flag', 'spatial', 'temporal', 'ranges', 'csv']},
            'typevirtual': {'type_param': 'str', 'list_values': ['spatial', 'temporal', 'ranges', 'csv']},
            'var_limit_to_valid': {'type_param': 'str', 'default': None},
            'limit_to_central_pixel': {'type_param': 'boolean', 'default': False},
            'flag_spatial_index': {'type_group': ['spatial'], 'type_param': 'strlist'},
            'lat_variable': {'type_group': ['spatial'], 'type_param': 'str', 'default': 'insitu_latitude'},
            'lon_variable': {'type_group': ['spatial'], 'type_param': 'str', 'default': 'insitu_longitude'},
            'time_variable': {'type_group': ['temporal'], 'type_param': 'str', 'default': 'insitu_time'},
            'time_ini': {'type_group': ['temporal'], 'type_param': 'str', 'default': None},
            'time_fin': {'type_group': ['temporal'], 'type_param': 'str', 'default': None},
            'time_flag_type': {'type_group': ['temporal'], 'type_param': 'str', 'default': 'ranges',
                               'list_values': ['ranges', 'yearjday', 'yearmonth', 'jday', 'month', 'yearmonthday',
                                               'monthday', 'flag_satellite_yearmonthday']},
            'var_ranges': {'type_group': ['ranges'], 'type_param': 'str'},
            'flag_ranges_indexm': {'type_group': ['ranges'], 'type_param': 'strlist'},
            'path_csv': {'type_group': ['csv'], 'type_param': 'file', 'default': None},
            'col_date': {'type_group': ['csv'], 'type_param': 'str', 'default': None},
            'format_date': {'type_group': ['csv'], 'type_param': 'str', 'default': '%Y-%m-%d'},
            'col_time': {'type_group': ['csv'], 'type_param': 'str', 'default': None},
            'format_time': {'type_group': ['csv'], 'type_param': 'str', 'default': '%H:%M'},
            'col_flag': {'type_group': ['csv'], 'type_param': 'str'},
            'var_ref_time': {'type_group': ['csv'], 'type_param': 'str', 'default': 'insitu_time'},
            'flag_list': {'type_group': ['csv'], 'type_param': 'strlist'}

        }

    def get_virtual_flag_list(self):
        sfinal = []
        for flag in self.options.sections():
            if self.options.has_option(flag, 'type'):
                if self.options[flag]['type'] == 'virtual_flag': sfinal.append(flag)
        if len(sfinal) == 0: sfinal = None
        return sfinal

    def get_virtual_flags_options(self):
        vf_options = {}
        fl = self.get_virtual_flag_list()
        if fl is not None:
            for f in fl: vf_options[f] = self.get_options_dict(f)
        return vf_options

    def get_options_dict(self, flag_ref):
        options_dict = self.read_options_as_dict(flag_ref, self.flag_options)
        return options_dict

    def create_flag_array_ranges_v2(self, options_dict):

        default_value = 0
        flag_names = []
        flag_values = []
        fl_info = {}
        array = None
        nscans = 1
        dataset = Dataset(self.nc_file)
        for fl in options_dict:
            if fl.startswith('flag_ranges'):
                if options_dict[fl]['is_default']:
                    default_value = int(options_dict[fl]['flag_value'])
                value = options_dict[fl]['flag_value']
                if value not in flag_values:
                    flag_values.append(value)
                    flag_names.append(options_dict[fl]['flag_name'])
                value_s = str(options_dict[fl]['flag_value'])
                if value_s not in fl_info.keys():
                    fl_info[value_s] = [fl]
                else:
                    fl_info[value_s].append(fl)
                if array is None:
                    r_array = np.array(dataset.variables[options_dict[fl]['flag_var']][:])
                    array = r_array.copy().astype(np.int64)
                if array is not None and len(dataset.variables[options_dict[fl]['flag_var']].dimensions) == 2:
                    r_array = np.array(dataset.variables[options_dict[fl]['flag_var']][:])
                    array = r_array.copy().astype(np.int64)
                    nscans = array.shape[1]

        array[:] = default_value

        for svalue in fl_info:
            value = int(svalue)
            nranges = len(fl_info[svalue])
            if nranges == 1:
                fl = fl_info[svalue][0]
                r_array = np.array(dataset.variables[options_dict[fl]['flag_var']][:])
                if nscans > 1 and len(r_array.shape) == 1:
                    r_array = self.replicate_scans(r_array, nscans)
                v_min = options_dict[fl]['min_range']
                v_max = options_dict[fl]['max_range']
                if v_min is not None and v_max is not None:
                    array[np.logical_and(r_array >= v_min, r_array <= v_max)] = value
                elif v_min is None and v_max is not None:
                    array[r_array <= v_max] = value
                elif v_min is not None and v_max is None:
                    array[r_array >= v_min] = value
            else:
                array_check = np.zeros(array.shape)
                for fl in fl_info[svalue]:
                    r_array = np.array(dataset.variables[options_dict[fl]['flag_var']][:])
                    if nscans > 1 and len(r_array.shape) == 1:
                        r_array = self.replicate_scans(r_array, nscans)
                    v_min = options_dict[fl]['min_range']
                    v_max = options_dict[fl]['max_range']
                    if v_min is not None and v_max is not None:
                        array_check[np.logical_and(r_array >= v_min, r_array <= v_max)] = array_check[np.logical_and(
                            r_array >= v_min, r_array <= v_max)] + 1
                    elif v_min is None and v_max is not None:
                        array_check[r_array <= v_max] = array_check[r_array <= v_max] + 1
                    elif v_min is not None and v_max is None:
                        array_check[r_array >= v_min] = array_check[r_array >= v_min] + 1
                condition = options_dict[fl]['flag_condition']
                if condition == 'and':
                    array[array_check == nranges] = value
                elif condition == 'or':
                    array[array_check > 0] = value
        dataset.close()

        return array, flag_names, flag_values

    def replicate_scans(self, array, nscans):
        ndata = array.shape[0]
        new_array = np.zeros((ndata, nscans))
        for iscan in range(nscans):
            new_array[:, iscan] = array[:]
        return new_array

    def get_flag_value(self, vorig, use_pow2_vflags):
        if use_pow2_vflags:
            fvalue = int(math.pow(2, vorig))
        else:
            fvalue = int(vorig + 1)
        return fvalue

    ##when type option is selected, then the rest of options are only applied if type_group=type
    def read_options_as_dict(self, section, poptions):
        options_dict = {}
        type = None
        use_pow2_flags = False
        if self.options.has_option(section, 'type'):
            type = self.read_option(section, 'type', 'type', poptions, -1, use_pow2_flags)
        if (type == 'virtual_flag' or type is None) and self.options.has_option(section,
                                                                                'typevirtual'):  ##typevirtual overwrites type
            type = self.read_option(section, 'typevirtual', 'typevirtual', poptions, -1, use_pow2_flags)
        if type is None:
            print(
                f'[ERROR] type (or typevirtual) are not define. Potential types: {poptions["typevirtual"]["list_values"]}')
            return None
        if self.options.has_option(section, 'use_pow2_flags'):
            use_pow2_flags = self.options[section]['use_pow2_flags']
        options_dict = self.assign_options(options_dict, 'type', type)
        options_dict = self.assign_options(options_dict, 'use_pow2_flags', use_pow2_flags)
        for opt in poptions:
            if opt == 'type' or opt == 'typevirtual' or opt == 'use_pow2_flags':
                continue
            if type is not None and 'type_group' in poptions[opt]:
                if type not in poptions[opt]['type_group']:
                    continue
            if not opt.find('_index') >= 0:
                if self.options.has_option(section, opt):
                    value = self.read_option(section, opt, opt, poptions, -1, use_pow2_flags)
                    options_dict = self.assign_options(options_dict, opt, value)
                else:
                    if 'default' in poptions[opt].keys():
                        options_dict[opt] = poptions[opt]['default']
            else:
                if opt.find('_indexm') >= 0:
                    idx = 0
                    has_option = True
                    while has_option:
                        inner_idx = 0
                        has_inner = True
                        while has_inner:
                            opt_here = opt.replace('_indexm', f'_{idx}_{inner_idx}')
                            if self.options.has_option(section, opt_here):
                                value = self.read_option(section, opt, opt_here, poptions, idx, use_pow2_flags)
                                options_dict = self.assign_options(options_dict, opt_here, value)
                                inner_idx = inner_idx + 1
                            else:
                                has_inner = False
                            if inner_idx == 0:
                                has_option = False
                        idx = idx + 1
                else:
                    idx = 0
                    has_option = True
                    while has_option:
                        opt_here = opt.replace('_index', f'_{idx}')
                        if self.options.has_option(section, opt_here):
                            value = self.read_option(section, opt, opt_here, poptions, idx, use_pow2_flags)
                            options_dict = self.assign_options(options_dict, opt_here, value)
                        else:
                            has_option = False
                        idx = idx + 1

        return options_dict

    def assign_options(self, options_dict, key, value):
        keys = key.split('.')
        if len(keys) == 1:
            options_dict[key] = value
        elif len(keys) == 2:
            if keys[0] in options_dict.keys():
                options_dict[keys[0]][keys[1]] = value
            else:
                options_dict[keys[0]] = {keys[1]: value}
        return options_dict

    def read_option(self, section, opt, opt_here, poptions, idx, use_pow2_flags):
        default = None
        if 'default' in poptions[opt].keys():
            default = poptions[opt]['default']
        value = self.get_value_param(section, opt_here, default, poptions[opt]['type_param'])

        if 'list_values' in poptions[opt]:
            if value not in poptions[opt]['list_values']:
                value = None

        return value

    def get_value(self, section, key):
        value = None
        if self.options.has_option(section, key):
            value = self.options[section][key]
        return value

    def get_value_param(self, section, key, default, type):
        value = self.get_value(section, key)

        if value is None:
            return default
        if type == 'str':
            return value
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
        if type == 'intlist':
            list_str = value.split(',')
            list = []
            for vals in list_str:
                vals = vals.strip()
                list.append(int(vals))
            return list