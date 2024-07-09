import os, math
import numpy as np
from netCDF4 import Dataset
from flag_manager import Flags
from flag_manager import OptionsManager
from datetime import datetime as dt


class FlagBuilder:

    def __init__(self, path_mdb_file, options):
        self.nc_file = path_mdb_file
        self.VALID = os.path.isfile(path_mdb_file)

        if not self.VALID:
            print(f'[ERROR] {path_mdb_file} is not a valid  file')

        self.omanager = OptionsManager(options)

        self.flag_list = self.omanager.get_section_list(None)
        if self.flag_list is None:
            self.VALID = False
            print(f'[ERROR] Flag configuration could not be started')

        self.flag_options = {
            'type': {'type_param': 'str', 'list_values': ['virtual_flag', 'spatial', 'temporal', 'ranges', 'csv']},
            'typevirtual': {'type_param': 'str', 'list_values': ['spatial', 'temporal', 'flag', 'ranges', 'csv']},
            'level': {'type_param': 'str', 'list_values': ['l1', 'l2']},
            'var_limit_to_valid': {'type_param': 'str', 'default': None},
            'limit_to_central_pixel': {'type_param': 'boolean', 'default': False},
            'output_dim': {'type_param': 'str', 'default': 'satellite_id'},
            'default_flag': {'type_param': 'str', 'default': 'DEFAULT'},
            'default_value': {'type_param': 'int', 'default': 0},
            'default_masked': {'type_param': 'boolean', 'default': False},
            'flag_spatial_index': {'type_group': ['spatial'], 'type_param': 'strlist'},
            'lat_variable': {'type_group': ['spatial'], 'type_param': 'str', 'default': 'insitu_latitude'},
            'lon_variable': {'type_group': ['spatial'], 'type_param': 'str', 'default': 'insitu_longitude'},
            'time_variable': {'type_group': ['temporal'], 'type_param': 'str', 'default': 'insitu_time'},
            'time_ini_abs': {'type_group': ['temporal'], 'type_param': 'str', 'default': None},
            'time_fin_abs': {'type_group': ['temporal'], 'type_param': 'str', 'default': None},
            'time_flag_type': {'type_group': ['temporal'], 'type_param': 'str', 'default': 'ranges',
                               'list_values': ['ranges', 'instants']},
            'instant_format': {'type_group': ['temporal'], 'type_param': 'str', 'default': '%Y'},
            'flag_ranges_indexm': {'type_group': ['ranges'], 'type_param': 'strlist'},
            'flag_indexm': {'type_group': ['flag'], 'type_param': 'strlist'},
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
        return self.omanager.get_virtual_flag_list()

    def get_virtual_flags_options(self):
        vf_options = {}
        fl = self.get_virtual_flag_list()
        if fl is not None:
            for f in fl:
                vf_options[f] = self.get_options_dict(f)
        return vf_options

    def get_options_dict(self, flag_ref):
        options_dict = self.omanager.read_options_as_dict(flag_ref, self.flag_options)
        return options_dict

    def create_flag_array(self, flag_ref, create_copy):
        options_dict = self.omanager.read_options_as_dict(flag_ref, self.flag_options)
        if options_dict is None:
            return None, None, None

        type = options_dict['type']
        if type is None:
            type = options_dict['typevirtual']
        if type is None:
            return
        if type == 'spatial':
            array, dims, flag_names, flag_values = self.create_flag_array_spatial(options_dict)
            if create_copy:
                self.create_copy_with_flag_band(flag_ref, array, flag_names, flag_names, dims)

        if type == 'temporal':
            array, dims, flag_names, flag_values = self.create_flag_array_temporal_v2(options_dict)
            if create_copy:
                self.create_copy_with_flag_band(flag_ref, array, flag_names, flag_names, dims)

        if type == 'ranges':
            array, dims, flag_names, flag_values = self.create_flag_array_ranges(options_dict)
            if create_copy:
                self.create_copy_with_flag_band(flag_ref, array, flag_names, flag_names, dims)

        if type == 'csv':
            array, dims, flag_names, flag_values = self.create_flag_array_csv(options_dict)
            if create_copy:
                self.create_copy_with_flag_band(flag_ref, array, flag_names, flag_values, dims)

        if type == 'flag':
            array, dims, flag_names, flag_values = self.create_flag_array_flag(options_dict)
            if create_copy:
                self.create_copy_with_flag_band(flag_ref, array, flag_names, flag_values, dims)

        return flag_values, flag_names, array

    def check_level(self, dataset, level, required_arrays):
        # getting nseries and nscan
        nseries = -1
        nscan = -1
        for variable in dataset.variables:
            if nseries == -1 and variable.startswith('l1') and len(dataset.variables[variable].shape) == 2:
                nseries = dataset.variables[variable].shape[0]
                nscan = dataset.variables[variable].shape[1]
        # start output array
        if level == 'l1':
            array = np.ma.zeros((nseries, nscan))
        if level == 'l2':
            array = np.ma.zeros((nseries,))

        ##checking required arrays
        for name_var in required_arrays:
            array_check = required_arrays[name_var]['array']
            if level == 'l1' and len(array_check.shape) == 1 and array_check.shape[0] == nseries:
                ##convert l2 into level1
                array_check = self.replicate_scans(array_check, nscan)
            if level == 'l2' and len(array_check.shape) == 2 and array_check.shape[0] == nseries and array_check.shape[
                1] == nscan:
                ##convert l1 into level2
                array_check = np.squeeze(array_check[:, 0])
            required_arrays[name_var]['array'] = array_check

        return array, required_arrays

    def create_flag_array_temporal_v2(self, options_dict):

        dataset = Dataset(self.nc_file)

        default_value = options_dict['default_value']
        default_flag = options_dict['default_flag']
        var_time = options_dict['time_variable']

        use_pow2_flags = options_dict['use_pow2_flags']
        time_array = dataset.variables[var_time][:]
        if var_time=='insitu_time' and  options_dict['output_dim']=='satellite_id':
            pass

        array = time_array.copy().astype(np.int32)
        array[:] = default_value
        if default_flag is None or default_flag=='DEFAULT':
            flag_names = []
            flag_values = []
        else:
            flag_names = [default_flag]
            flag_values = [default_value]
        flag_type = options_dict['time_flag_type']

        if flag_type=='instants':
            ins_format = options_dict['instant_format']
            index_flag = 0
            for idx,t in enumerate(time_array):
                print(idx,t)
                str_instant = dt.utcfromtimestamp(float(t)).strftime(ins_format)
                if str_instant not in flag_names:
                    flag_names.append(str_instant)
                    flag_value = self.get_flag_value(index_flag,use_pow2_flags)
                    flag_values.append(flag_value)
                    index_flag = index_flag + 1
                else:
                    index_here = flag_names.index(str_instant)
                    flag_value = self.get_flag_value(index_here,use_pow2_flags)
                array[idx] = np.int64(flag_value)

        dataset.close()

        return array, None, flag_names, flag_values

    def create_flag_array_ranges_v2(self, options_dict):
        dataset = Dataset(self.nc_file)
        # print('****************************************************************************************************')
        # print(options_dict)
        # print('======================================================================================================')
        default_value = options_dict['default_value']
        default_flag = options_dict['default_flag']
        indices_ranges = {}
        required_arrays = {}
        for fl in options_dict:
            if fl.startswith('flag_ranges'):
                index_range = fl.split('_')[2]
                index_and = fl.split('_')[3]
                if not index_range in indices_ranges.keys():
                    indices_ranges[index_range] = {
                        'flag_name': options_dict[fl]['flag_name'],
                        'flag_value': options_dict[fl]['flag_value'],
                        index_and: options_dict[fl]['condition_list'],
                        'n_and': 1
                    }
                else:
                    indices_ranges[index_range][index_and] = options_dict[fl]['condition_list']
                    indices_ranges[index_range]['n_and'] = indices_ranges[index_range]['n_and'] + 1

                for condition in options_dict[fl]['condition_list']:
                    name_var = condition['name_var']
                    is_flag = True if condition['flag_or_range'] == 'flag' else False
                    if name_var not in required_arrays.keys():
                        if name_var in dataset.variables:
                            if is_flag:
                                flag_info = self.get_flag_info_from_variable(dataset.variables[name_var])
                                required_arrays[name_var] = flag_info
                            else:
                                required_arrays[name_var] = {
                                    'array': dataset.variables[name_var][:]
                                }
                        else:
                            if name_var in self.omanager.get_virtual_flag_list():
                                flag_values, flag_names, array = self.create_flag_array(name_var, False)
                                options_dict_f = self.omanager.read_options_as_dict(name_var, self.flag_options)
                                required_arrays[name_var] = {
                                    'array': array,
                                    'flag_meanings': ' '.join(flag_names),
                                    'flag_values': flag_values,
                                    'user_pow2_flags': options_dict_f['use_pow2_flags']
                                }

        ##checking level
        level = options_dict['level']
        if level == 'l1' or level == 'l2':
            array, required_arrays = self.check_level(dataset, level, required_arrays)

        # starting array
        array[:] = default_value
        if default_flag != 'DEFAULT':
            flag_names = [default_flag]
            flag_values = [default_value]
        else:
            flag_names = []
            flag_values = []

        # creating array
        for irange in indices_ranges:
            flag_value = indices_ranges[irange]['flag_value']
            flag_values.append(flag_value)
            flag_names.append(indices_ranges[irange]['flag_name'])
            array_here = np.zeros(array.shape)
            n_and = indices_ranges[irange]['n_and']
            for iand in range(n_and):
                iand_s = f'{iand:.0f}'
                array_or = None
                for orcondition in indices_ranges[irange][iand_s]:
                    array_info = required_arrays[orcondition['name_var']]
                    r_array = array_info['array']
                    if orcondition['flag_or_range'] == 'flag':
                        name_flag = orcondition['name_flag']
                        or_flag_list = array_info['flag_meanings'].split(' ')
                        if name_flag in or_flag_list:
                            index_name_flag = or_flag_list.index(name_flag)
                            v_flag_value = array_info['flag_values'][index_name_flag]
                            v_min = v_flag_value
                            v_max = v_flag_value
                        else:
                            print(f'[WARNING] {name_flag} is not defined in flag band {orcondition["name_var"]}')
                            v_min = -1
                            v_min = -1
                    else:
                        v_min = orcondition['min_val']
                        v_max = orcondition['max_val']
                    if array_or is None:
                        array_or = np.zeros(r_array.shape)
                    if v_min is not None and v_max is not None:
                        array_or[np.logical_and(r_array >= v_min, r_array <= v_max)] = 1
                    elif v_min is None and v_max is not None:
                        array_or[r_array <= v_max] = 1
                    elif v_min is not None and v_max is None:
                        array_or[r_array >= v_min] = 1

                array_here = array_here + array_or

            array[array_here == n_and] = array[array_here == n_and] + flag_value

        dataset.close()

        return array, flag_names, flag_values

    def replicate_scans(self, array, nscan):
        nseries = array.shape[0]
        array = np.repeat(array.reshape(1, nseries), nscan, axis=1).reshape((nseries, nscan))
        return array


    def create_flag_array_flag(self, options_dict):
        dataset = Dataset(self.nc_file)
        default_value = options_dict['default_value']
        default_flag = options_dict['default_flag']
        indices_ranges = {}
        required_arrays = {}
        ##check options
        for fl in options_dict:
            if fl.startswith('flag_'):
                index_range = fl.split('_')[1]
                index_and = fl.split('_')[2]
                if not index_range in indices_ranges.keys():
                    indices_ranges[index_range] = {
                        'flag_name': options_dict[fl]['flag_name'],
                        'flag_value': options_dict[fl]['flag_value'],
                        index_and: options_dict[fl]['condition_list'],
                        'n_and': 1
                    }
                else:
                    indices_ranges[index_range][index_and] = options_dict[fl]['condition_list'],
                    indices_ranges[index_range]['n_and'] = indices_ranges[index_range]['n_and'] + 1

                for condition in options_dict[fl]['condition_list']:
                    name_var = condition['flag_var']
                    if name_var not in required_arrays.keys():
                        if name_var in dataset.variables:
                            info = self.get_flag_info_from_variable(dataset.variables[name_var])
                            required_arrays[name_var] = info
                        else:
                            print('NO IMPLEMENTED YET')
                            # if name_var in self.omanager.get_virtual_flag_list():
                            #     print('aqui no deberiamos estar-->', name_var)
                            #     self.create_flag_array(name_var, False)
                            #     # print(options_vf)

        ##checking level
        level = options_dict['level']
        if level == 'l1' or level == 'l2':
            array, required_arrays = self.check_level(dataset, level, required_arrays)

        # starting array
        array[:] = default_value
        if default_flag != 'DEFAULT':
            flag_names = [default_flag]
            flag_values = [default_value]
        else:
            flag_names = []
            flag_values = []

        # creating array
        for irange in indices_ranges:
            flag_value = indices_ranges[irange]['flag_value']
            flag_values.append(flag_value)
            flag_names.append(indices_ranges[irange]['flag_name'])
            array_here = np.zeros(array.shape)
            n_and = indices_ranges[irange]['n_and']
            for iand in range(n_and):
                iand_s = f'{iand:.0f}'
                array_or = None
                for orcondition in indices_ranges[irange][iand_s]:
                    info_flag = required_arrays[orcondition['flag_var']]
                    if array_or is None:
                        array_or = self.get_array_or(info_flag, orcondition['flag_list'])
                    else:
                        array_or_here = self.get_array_or(info_flag, orcondition['flag_list'])
                        array_or[array_or_here == 1] = 1

                array_here = array_here + array_or

            array[array_here == n_and] = array[array_here == n_and] + flag_value

        dims = None
        if level == 'l1':
            dims = ('series', 'scan')
        if level == 'l2':
            dims = ('series',)
        dataset.close()

        return array, dims, flag_names, flag_values

    # def create_flag_array_ranges(self, options_dict):
    #     var_ranges = options_dict['var_ranges']
    #     r_array = self.mfile.get_full_array(var_ranges)
    #     array = r_array.copy().astype(np.int32)
    #     default_value = 0
    #     flag_names = []
    #     flag_values = []
    #     for fl in options_dict:
    #         if fl.startswith('flag_ranges'):
    #             if options_dict[fl]['is_default']:
    #                 default_value = int(options_dict[fl]['flag_value'])
    #             flag_values.append(options_dict[fl]['flag_value'])
    #             flag_names.append(options_dict[fl]['flag_name'])
    #     ##Default value
    #     array[:] = default_value
    #     ##Flag limites by lat-long values
    #     for fl in options_dict:
    #         if fl.startswith('flag_ranges') and not options_dict[fl]['is_default']:
    #             v_min = options_dict[fl]['min_range']
    #             v_max = options_dict[fl]['max_range']
    #             value = int(options_dict[fl]['flag_value'])
    #             array[np.logical_and(r_array >= v_min, r_array < v_max)] = value
    #     dims = self.mfile.get_dims_variable(var_ranges)
    #
    #     if default_value != 0:
    #         fillValue = self.mfile.get_fill_value(var_ranges)
    #         if fillValue is not None:
    #             array[r_array == fillValue] = 0
    #
    #     # for idx in range(1800):
    #     #     print(idx,':',r_array[0][idx],'-->',array[0][idx],' with default: ',default_value)
    #
    #     return array, dims, flag_names, flag_values
    #
    # def create_flag_array_temporal(self, options_dict):
    #
    #     var_time = options_dict['time_variable']
    #
    #     flag_type = options_dict['time_flag_type']
    #     use_pow2_vflags = options_dict['use_pow2_flags']
    #
    #     time_array = self.mfile.get_full_array(var_time)
    #     array = time_array.copy().astype(np.int32)
    #     flag_names = []
    #     flag_values = []
    #
    #     # default_value = 0
    #
    #     orig_shape = array.shape
    #     if len(orig_shape) > 1:
    #         array1D = array.flatten()
    #     else:
    #         array1D = array
    #
    #     ##limit to central pixels
    #     array1Dcentral = None
    #     if options_dict['limit_to_central_pixel']:
    #         array_central = self.mfile.get_full_array('insitu_spatial_index')
    #         if array_central.shape == orig_shape:
    #             if len(orig_shape) > 1:
    #                 array1Dcentral = array_central.flatten()
    #             else:
    #                 array1Dcentral = array_central
    #
    #     ##limit to valid
    #     array1Dvalid = None
    #     var_limit_to_valid = options_dict['var_limit_to_valid']
    #     if var_limit_to_valid is not None:
    #         array_valid = self.mfile.get_full_array(var_limit_to_valid)
    #         if array_valid.shape == orig_shape:
    #             if len(orig_shape) > 1:
    #                 array1Dvalid = array_valid.flatten()
    #             else:
    #                 array1Dvalid = array_valid
    #
    #     # print('°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°')
    #     # array1Dvalid[6990] = 1
    #
    #     ##time ini + time_fin
    #     time_ini = None
    #     time_fin = None
    #     if options_dict['time_ini'] is not None and options_dict['time_fin'] is not None:
    #         time_ini = dt.strptime(options_dict['time_ini'], '%Y-%m-%d %H:%M').replace(tzinfo=pytz.utc)
    #         time_fin = dt.strptime(options_dict['time_fin'], '%Y-%m-%d %H:%M').replace(tzinfo=pytz.utc)
    #     fvalue = 0
    #     for idx in range(len(array1D)):
    #         val = array1D[idx]
    #
    #         if val < 0:
    #             array1D[idx] = 0
    #             continue
    #
    #         if array1Dvalid is not None:
    #             if array1Dvalid[idx] != 1:
    #                 continue
    #
    #         if array1Dcentral is not None:
    #             if array1Dcentral[idx] != 0:
    #                 continue
    #         time_here = dt.utcfromtimestamp(val).replace(tzinfo=pytz.utc)
    #
    #         if time_ini is not None and time_fin is not None:
    #             if time_here < time_ini or time_here > time_fin:
    #                 continue
    #         # if time_ini is not None and time_fin is not None:
    #         #     if time_ini < time_here < time_fin:
    #         #         continue
    #
    #         if flag_type == 'yearjday':
    #             time_here_str = time_here.strftime('%Y%j')
    #         elif flag_type == 'yearmonth':
    #             time_here_str = time_here.strftime('%Y%m')
    #         elif flag_type == 'jday':
    #             time_here_str = time_here.strftime('%j')
    #         elif flag_type == 'month':
    #             time_here_str = time_here.strftime('%m')
    #         elif flag_type.endswith('yearmonthday'):
    #             time_here_str = time_here.strftime('%Y-%m-%d')
    #         elif flag_type == 'monthday':
    #             time_here_str = time_here.strftime('%m-%d')
    #
    #         if flag_type.startswith('flag_satellite'):
    #             array_flag_satellite = self.mfile.get_full_array('flag_satellite')
    #             flag_satellite_here = array_flag_satellite[idx]
    #             if flag_satellite_here == 1:
    #                 time_here_str = f'S3A_{time_here_str}'
    #             if flag_satellite_here == 2:
    #                 time_here_str = f'S3B_{time_here_str}'
    #
    #         if time_here_str not in flag_names:
    #             flag_names.append(time_here_str)
    #             vflag = self.get_flag_value(fvalue, use_pow2_vflags)
    #             flag_values.append(vflag)
    #             print(f'[INFO] Adding temporal flag: {time_here_str} with value: {vflag}----->{val}')
    #             fvalue = fvalue + 1
    #         else:
    #             index_flag = flag_names.index(time_here_str)
    #             vflag = flag_values[index_flag]
    #
    #         array1D[idx] = vflag
    #
    #     if len(orig_shape) > 1:
    #         array = array1D.reshape(orig_shape)
    #     else:
    #         array = array1D
    #
    #     dims = self.mfile.get_dims_variable(var_time)
    #     # if default_value != 0:
    #     #     array[lat_array == -999.0] = 0
    #
    #     return array, dims, flag_names, flag_values
    #
    # def create_flag_array_csv(self, options_dict):
    #     if not self.VALID:
    #         print(f'[ERROR] Flag builder instance is not valid')
    #         return None
    #
    #     path_csv = options_dict['path_csv']
    #     if path_csv is None:
    #         print(f'[ERROR] path_csv is not available or does not exist. Review configuration file.')
    #         return None
    #     col_date = options_dict['col_date']
    #     if col_date is None:
    #         print(f'[ERROR] col_date option is required. Review configuration file.')
    #         return None
    #     col_flag = options_dict['col_flag']
    #     if col_flag is None:
    #         print(f'[ERROR] col_flag option is required. Review configuration file.')
    #         return None
    #     import pandas as pd
    #     try:
    #         dset = pd.read_csv(path_csv, sep=';')
    #     except:
    #         print(f'[ERROR] Problems reading CSV path: {path_csv}')
    #         return None
    #     col_names = dset.columns.tolist()
    #     if col_date not in col_names:
    #         print(f'[ERROR] col_date {col_date} is not avaialable in the CSV file')
    #         return None
    #     if col_flag not in col_names:
    #         print(f'[ERROR] col_date {col_flag} is not avaialable in the CSV file')
    #         return None
    #     col_time = options_dict['col_time']
    #     if col_time is not None and col_time not in col_names:
    #         print(f'[ERROR] col_time {col_time} is not avaialable in the CSV file')
    #         return None
    #
    #     date_csv = dset[col_date]
    #     time_csv = dset[col_time] if col_time is not None else None
    #     flag_csv = dset[col_flag]
    #     all_flag_list = np.unique(flag_csv).tolist()
    #
    #     if 'flag_list' in options_dict and options_dict['flag_list'] is not None:
    #         flag_list = options_dict['flag_list']
    #         for f in flag_list:
    #             if f not in all_flag_list:
    #                 print(f'[ERROR] Flag: {f} is not in the flag list: {all_flag_list} in col_flag {col_flag}')
    #                 return None
    #     else:
    #         flag_list = all_flag_list
    #
    #     if self.mfile is None:
    #         dataset = Dataset(self.path_mdb_file)
    #     else:
    #         dataset = self.mfile.nc
    #     var_time = options_dict['var_ref_time']
    #     if var_time not in dataset.variables:
    #         print(f'[ERROR] {var_time} is not included in the input dataset.')
    #         if self.mfile is None:
    #             dataset.close()
    #         return None
    #     time_array = dataset.variables[var_time][:]
    #     array = np.zeros(time_array.shape)
    #     array[time_array.mask] = -999
    #
    #     format_date = options_dict['format_date']
    #     if col_time is not None:
    #         format_date = f'{format_date}T{options_dict["format_time"]}'
    #     print(f'[INFO] Format date: {format_date}')
    #     list_times = {}
    #     if len(time_array.shape) == 2:
    #         for isatellite in range(time_array.shape[0]):
    #             for iinsitu in range(time_array.shape[1]):
    #                 if not np.ma.is_masked(time_array[isatellite, iinsitu]):
    #                     ts = float(time_array[isatellite, iinsitu])
    #                     ts_s = dt.utcfromtimestamp(ts).strftime(format_date)
    #                     list_times[ts_s] = {
    #                         'satellite_id': isatellite,
    #                         'insitu_id': iinsitu,
    #                         'assigned': False
    #                     }
    #     elif len(time_array.shape) == 1:
    #         for isatellite in range(time_array.shape[0]):
    #             if not np.ma.is_masked(time_array[isatellite]):
    #                 ts = float(time_array[isatellite])
    #                 ts_s = dt.utcfromtimestamp(ts).strftime(format_date)
    #                 list_times[ts_s] = {
    #                     'satellite_id': isatellite,
    #                     'assigned': False
    #                 }
    #
    #     print(f'[INFO] Retrieved {len(list_times)} time stamps')
    #     use_pow2_flags = options_dict['use_pow2_flags']
    #     if use_pow2_flags:
    #         flag_values = [2 ** x for x in range(len(flag_list))]
    #     else:
    #         flag_values = [x + 1 for x in range(len(flag_list))]
    #     nassigned = 0
    #     for index, date_c in enumerate(date_csv):
    #         if col_time is not None:
    #             date_c = f'{date_c}T{time_csv[index]}'
    #         flag_here = flag_csv[index]
    #         index_flag_here = flag_list.index(flag_here)
    #         if date_c in list_times:
    #             nassigned = nassigned + 1
    #             if len(time_array.shape) == 2:
    #                 satellite_id = list_times[date_c]['satellite_id']
    #                 insitu_id = list_times[date_c]['insitu_id']
    #                 array[satellite_id, insitu_id] = np.int32(flag_values[index_flag_here])
    #                 list_times[date_c]['assigned'] = True
    #             elif len(time_array.shape) == 1:
    #                 satellite_id = list_times[date_c]['satellite_id']
    #                 array[satellite_id] = np.int32(flag_values[index_flag_here])
    #                 list_times[date_c]['assigned'] = True
    #     print(f'[INFO] Number of assigned values: {nassigned}')
    #
    #     dims = dataset.variables[var_time].dimensions
    #     if self.mfile is None:
    #         dataset.close()
    #
    #     return array, dims, flag_list, flag_values

    def get_flag_info_from_variable(self, variable):
        array = variable[:]
        flag_meanings = None
        if 'flag_meanings' in variable.ncattrs():
            if isinstance(flag_meanings, list):
                flag_meanings = ' '.join(flag_meanings)
            else:
                flag_meanings = variable.flag_meanings
        if 'flag_values' in variable.ncattrs():
            flag_values = np.array([int(x) for x in variable.flag_values.split(',')])
        elif 'flag_masks' in variable.ncattrs():
            flag_values = np.array([int(x) for x in variable.flag_masks.split(',')])
        else:
            flag_values = np.unique(array).tolist()
        # print(flag_values, type(flag_values))
        flag_values = np.sort(flag_values).tolist()
        use_pow2_flags = False
        if len(flag_values) >= 3 and flag_values[2] == 4:
            use_pow2_flags = True

        return {'array': array, 'flag_meanings': flag_meanings, 'flag_values': flag_values,
                'user_pow2_flags': use_pow2_flags}

    def get_array_or(self, flag_info, flag_list):
        flag_list_real = []
        compute_all = False
        is_all = False
        for flag in flag_list:
            if flag == '[NONE]' or flag == '[ALL]':
                compute_all = True
                if flag == '[ALL]': is_all = True
            else:
                flag_list_real.append(flag)

        flag_array = flag_info['array']
        fflags = Flags(flag_info['flag_values'], flag_info['flag_meanings'])
        if compute_all:
            array_all = fflags.Mask(flag_array, flag_info['flag_meanings'].split(' '))
            if is_all:  ##or condition, if ALL selected, or arry is array_all
                array_all[array_all > 0] = 1
                return array_all

        if len(flag_list_real) == 0 and compute_all:
            array_here = np.zeros(array_all.shape)
            array_here[array_all == 0] = 1
            array_here[array_all > 0] = 0
            return array_here

        if len(flag_list_real) > 0:
            array_here = fflags.Mask(flag_array, flag_list_real)
            array_here[array_here > 0] = 1
            if compute_all:
                array_here[array_all == 0] = 1
            return array_here

    def get_flag_value(self, vorig, use_pow2_vflags):
        if use_pow2_vflags:
            fvalue = int(math.pow(2, vorig))
        else:
            fvalue = int(vorig + 1)
        return fvalue

    # def create_flag_array_spatial(self, options_dict):
    #     var_latitude = options_dict['lat_variable']
    #     var_longitude = options_dict['lon_variable']
    #     lat_array = self.mfile.get_full_array(var_latitude)
    #     lon_array = self.mfile.get_full_array(var_longitude)
    #     array = lat_array.copy().astype(np.int32)
    #     default_value = 0
    #     flag_names = []
    #     flag_values = []
    #     for fl in options_dict:
    #         if fl.startswith('flag_spatial'):
    #             if options_dict[fl]['is_default']:
    #                 default_value = int(options_dict[fl]['flag_value'])
    #             flag_values.append(options_dict[fl]['flag_value'])
    #             flag_names.append(options_dict[fl]['flag_name'])
    #     ##Default value
    #     array[:] = default_value
    #     ##Flag limites by lat-long values
    #     for fl in options_dict:
    #         if fl.startswith('flag_spatial') and not options_dict[fl]['is_default']:
    #             lat_min = options_dict[fl]['lat_min']
    #             lat_max = options_dict[fl]['lat_max']
    #             lon_min = options_dict[fl]['lon_min']
    #             lon_max = options_dict[fl]['lon_max']
    #             value = int(options_dict[fl]['flag_value'])
    #             array[np.logical_and(np.logical_and(lat_array >= lat_min, lat_array <= lat_max),
    #                                  np.logical_and(lon_array >= lon_min, lon_array <= lon_max))] = value
    #     dims = self.mfile.get_dims_variable(var_longitude)
    #     if default_value != 0:
    #         array[lat_array == -999.0] = 0
    #
    #     return array, dims, flag_names, flag_values
    #
    # def create_copy_with_flag_band(self, name_flag, array, flag_list, flag_values, dims):
    #     folder = os.path.dirname(self.path_mdb_file)
    #     time_now = dt.now().strftime('%Y%m%d%H%M%S')
    #     temp_file = os.path.join(folder, f'Temp_{time_now}.nc')
    #     ncout = Dataset(temp_file, 'w', format='NETCDF4')
    #     file_nc = Dataset(self.path_mdb_file, 'r')
    #
    #     # copy global attributes all at once via dictionary
    #     ncout.setncatts(file_nc.__dict__)
    #
    #     # copy dimensions
    #     for name, dimension in file_nc.dimensions.items():
    #         ncout.createDimension(
    #             name, (len(dimension) if not dimension.isunlimited() else None))
    #
    #     # copy variables
    #     for name, variable in file_nc.variables.items():
    #         fill_value = None
    #         if '_FillValue' in list(file_nc.ncattrs()):
    #             fill_value = variable._FillValue
    #
    #         ncout.createVariable(name, variable.datatype, variable.dimensions, fill_value=fill_value, zlib=True,
    #                              complevel=6)
    #         # copy variable attributes all at once via dictionary
    #         ncout[name].setncatts(file_nc[name].__dict__)
    #
    #         ncout[name][:] = file_nc[name][:]
    #
    #     # flag_variable
    #     if name_flag in file_nc.variables:
    #         print(f'[INFO] Flag name {name_flag} is already in the original file')
    #     else:
    #         var = ncout.createVariable(name_flag, 'i4', dims, zlib=True, complevel=6, fill_value=-999)
    #
    #     flag_meanings = ' '.join(flag_list)
    #     var.flag_meanings = flag_meanings
    #     var.flag_values = flag_values
    #     var[:] = array[:]
    #
    #     file_nc.close()
    #     ncout.close()
    #
    #     os.rename(temp_file, self.path_mdb_file)
