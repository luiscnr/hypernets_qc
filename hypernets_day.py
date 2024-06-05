from datetime import datetime as dt
import os, subprocess, pytz
from netCDF4 import Dataset
from plot_multiple import PlotMultiple
from hypernets_day_file import HYPERNETS_DAY_FILE

class HYPERNETS_DAY:

    def __init__(self, path_data, path_output):
        if path_data is not None:
            self.path_data = path_data
        else:
            self.path_data = '/mnt/c/DATA_LUIS/INSITU_HYPSTAR'
        if path_output is not None:
            self.path_output = path_output
        else:
            self.path_output = self.path_data
        print(f'[INFO] Started HYPERNETS_DAY with:')
        print(f'[INFO] ->Input path data: {self.path_data}')
        print(f'[INFO] ->Output path data: {self.path_output}')

        rsync_user = 'hypstar'
        self.url_base = f'{rsync_user}@enhydra.naturalsciences.be'
        self.url_base_npl = f'cnr@hypernetssvr1.npl.co.uk'
        self.base_folder = '/waterhypernet/hypstar/processed_v2/'
        self.base_folder_l2_npl = '/home/cnr/'
        self.base_folder_l2_rbins = '/home/hypstar/'

        self.ssh_base = 'ssh -X -Y -p 9022'
        self.ssh_base_npl = 'ssh'
        self.rsync_base = f'rsync -a -e \'ssh -p 9022\' {self.url_base}:{self.base_folder}'
        self.rsync_url = f'rsync -a -e \'ssh -p 9022\' {self.url_base}'
        self.rsync_url_npl = f'rsync -a -e \'ssh\' {self.url_base_npl}'

        self.n_series = 50
        self.files_dates = {}
        self.dataset_w = None

    def get_files_date_nodownload(self, site, date_here):
        self.files_dates = {}
        date_folder = self.get_folder_date(site, date_here)
        if date_folder is None:
            return
        list_sequences = self.get_sequences_date_from_file_list(site, date_here)
        if len(list_sequences) == 0:
            print(f'[WARNING] No sequences found for date: {date_here}')
            return
        list_seq_refs = [x[3:-2] for x in list_sequences]

        for name in os.listdir(date_folder):
            if name.find('L2A_REF') > 0 and name.endswith('nc'):
                sequence_ref = name.split('_')[5]
                try:
                    list_seq_refs.index(sequence_ref)
                    if sequence_ref not in self.files_dates.keys():
                        self.files_dates[sequence_ref] = {
                            'file_l2': os.path.join(date_folder, name),
                            'file_l1': None,
                            'file_images': None,
                            'valid': True
                        }
                    else:
                        self.files_dates[sequence_ref]['file_l2'] = os.path.join(date_folder, name)
                except:
                    pass
            if name.find('L1C_ALL') > 0 and name.endswith('nc'):
                sequence_ref = name.split('_')[5]
                try:
                    list_seq_refs.index(sequence_ref)
                    if sequence_ref not in self.files_dates.keys():
                        self.files_dates[sequence_ref] = {
                            'file_l2': None,
                            'file_l1': os.path.join(date_folder, name),
                            'file_images': None,
                            'valid': True
                        }
                    else:
                        self.files_dates[sequence_ref]['file_l1'] = os.path.join(date_folder, name)
                except:
                    pass
            if name.find('IMG') > 0 and name.endswith('jpg'):
                sequence_ref = name.split('_')[4]
                try:
                    list_seq_refs.index(sequence_ref)
                    if sequence_ref not in self.files_dates.keys():
                        self.files_dates[sequence_ref] = {
                            'file_l2': None,
                            'file_l1': None,
                            'file_images': [os.path.join(date_folder, name)],
                            'valid': True
                        }
                    else:
                        file_images = self.files_dates[sequence_ref]['file_images']
                        if file_images is None:
                            file_images = [os.path.join(date_folder, name)]
                        else:
                            file_images.append(os.path.join(date_folder, name))
                        self.files_dates[sequence_ref]['file_images'] = file_images
                except:
                    pass

    def get_sun_images_date(self, site, date_here,ndw):
        sun_images = {}

        date_folder = self.get_folder_date(site, date_here)
        if date_folder is not None:
            nimages = 0
            for name in os.listdir(date_folder):
                if name.find('_0_0') > 0 and name.endswith('.jpg'):
                    seq = f'{name.split("_")[4]}00'
                    sun_images[seq] = os.path.join(date_folder, name)
                    nimages = nimages + 1

        if ndw:
            return sun_images

        if date_folder is None:
            path_site = os.path.join(self.path_data, site)
            path_year = os.path.join(path_site, date_here.strftime('%Y'))
            path_month = os.path.join(path_year, date_here.strftime('%m'))
            path_day = os.path.join(path_month, date_here.strftime('%d'))
            self.create_if_not_exists(path_site)
            self.create_if_not_exists(path_year)
            self.create_if_not_exists(path_month)
            self.create_if_not_exists(path_day)
            date_folder = path_day

        if site == 'JSIT':
            list_sequences = self.get_sequences_date_l2(site, date_here)
        else:
            list_sequences = self.get_sequences_date(site, date_here)
            if len(list_sequences) == 0:
                list_sequences = self.get_sequences_date_l2(site, date_here)

        if len(list_sequences) == 0:
            print(f'[WARNING] No sequences found for date: {date_here}')
            return

        for seq in list_sequences:
            print(f'[INFO] Checking sun picture in sequence: {seq}')
            # seq_time =  dt.strptime(seq[3:], '%Y%m%dT%H%M%S').replace(tzinfo=pytz.UTC)
            files_sun = self.get_sun_image_sequence(site, seq)

            if len(files_sun) == 1:
                name_file = files_sun[0].split('/')[-1]
                file_sun = os.path.join(date_folder, name_file)
                if not os.path.exists(file_sun):
                    if site=='JSIT':
                        self.transfer_file_ssh_npl(files_sun[0],file_sun)
                    else:
                        self.transfer_file_ssh(files_sun[0], file_sun)
                if os.path.exists(file_sun):
                    name_new = self.get_name_new_file_sun(file_sun,site,seq)
                    file_sun_new = os.path.join(date_folder,name_new)
                    os.rename(file_sun,file_sun_new)
                    sun_images[seq[3:]] = file_sun_new
                else:
                    sun_images[seq[3:]] = None

        return sun_images

    def get_name_new_file_sun(self,file_img,site,seq):
        type = 'W'
        picture = '016'
        if site=='JSIT':
            type = 'L'
            picture = '090'
        date_img =  dt.fromtimestamp(os.path.getmtime(file_img))
        date_img_str = date_img.strftime('%Y%m%dT%H%M')
        name_new = f'HYPTERNETS_{type}_{site}_IMG_{seq[3:-2]}_{date_img_str}_{picture}_0_0_v2.0.jpg'
        return name_new

    def get_files_date(self, site, date_here):
        self.files_dates = {}
        date_folder = self.get_folder_date(site, date_here)
        if date_folder is None:
            return

        list_sequences = self.get_sequences_date(site, date_here)
        if len(list_sequences) == 0:
            print(f'[WARNING] No sequences found for date: {date_here}')
            return

        for name in os.listdir(date_folder):
            if name.find('L2A_REF') > 0:
                sequence_ref = name.split('_')[5]
                prefix = f'SEQ{sequence_ref}'
                sequence_folder = self.get_sequence_folder_from_list(list_sequences, prefix)
                ssh_path = self.get_ssh_path(site, date_here, sequence_folder)
                path_image = f'{ssh_path}/image'
                files_images = self.get_images_sequence(site, date_here, sequence_folder)
                print(f'[INFO] Checking sequence: {sequence_folder}')
                file_l1 = os.path.join(date_folder, name.replace('L2A_REF', 'L1C_ALL'))
                files_images_folder = [os.path.join(date_folder, x) for x in files_images]
                self.files_dates[sequence_ref] = {
                    'file_l2': os.path.join(date_folder, name),
                    'file_l1': file_l1,
                    'file_images': files_images_folder,
                    'valid': True
                }
                for img in files_images_folder:
                    if not os.path.exists(img):
                        name_img = img.split('/')[-1]
                        print(f'[INFO] --> File {name_img} was not found. Launching download...')
                        self.transfer_file_ssh(f'{path_image}/{name_img}', img)
                if not os.path.exists(file_l1):
                    name_l1 = file_l1.split('/')[-1]
                    print(f'[INFO] --> File {name_l1} was not found. Launching download...')
                    self.transfer_file_ssh(f'{ssh_path}/{name_l1}', file_l1)

    def save_sun_images(self, output_file, sun_images,time_list):
        pm = PlotMultiple()
        nrow = len(sun_images)
        ncol = len(sun_images[0])
        xfigsize = 2.2 * ncol
        yfigsize = 2.7 * nrow



        pm.start_multiple_plot_advanced(nrow, ncol, xfigsize, yfigsize, 0, 0.15, True)
        for irow in range(nrow):
            for icol in range(ncol):
                file_img = sun_images[irow][icol]
                if file_img is not None:
                    name = file_img.split('/')[-1]
                    title = name.split('_')[4]
                else:
                    title = time_list[irow][icol]
                if file_img is not None:
                    # pm.plot_image_title(file_img,index_row,index_col,title)
                    pm.plot_image_hypernets(file_img, irow, icol, title)
                else:
                    pm.plot_blank_with_title(irow, icol, title)

        print(f'[INFO] Saving sun plot to: {output_file}')
        pm.save_fig(output_file)
        pm.close_plot()

    def get_file_date_complete(self, site, date_here):
        folder_date = self.get_output_folder_date(site, date_here)
        if folder_date is None:
            return None
        date_here_str = date_here.strftime('%Y%m%d')
        file_date = os.path.join(folder_date, f'HYPERNETES_W_DAY_{date_here_str}.nc')

        return file_date

    def get_hypernets_day_file(self, site, date_here):


        file_date = self.get_file_date_complete(site, date_here)
        if file_date is None:
            return None
        if os.path.exists(file_date):
            return HYPERNETS_DAY_FILE(file_date, self.path_data)
        else:
            return None

    def start_file_date_complete(self, site, date_here, overwrite):


        file_date = self.get_file_date_complete(site, date_here)
        if file_date is None:
            print(f'[WARNING] Date folder for {site} and {date_here} is not avaiable. Skipping...')
            return -999
        if os.path.exists(file_date) and not overwrite:
            print(f'[WARNING] File: {file_date} already exists. Skipping')
            return -100
        if len(self.files_dates) == 0:
            print(f'[WARNING] Dates are not avaiable for this date. Skipping...')
            return -1
        seq_list = list(self.files_dates.keys())
        seq_list.sort()

        index_seq_ref = -1
        dims = None
        nseq_valid = 0
        for iseq in range(len(seq_list)):
            seq = seq_list[iseq]
            dims_here = self.check_dimensions(seq)
            if dims_here is not None:
                dims = dims_here
                index_seq_ref = iseq
                nseq_valid = nseq_valid + 1
            else:
                self.files_dates[seq]['valid'] = False

        if dims is None:
            print(f'[WARNING] L1 files are not available for this date')
            return -2

        try:
            self.dataset_w = Dataset(file_date, 'w')
        except:
            print(f'[WARNING] File {file_date} can not be created. Please check folder permissions')
            return -3
        ##dimensions
        print(f'[INFO] Creating dimensions...')
        self.dataset_w.createDimension('series')
        self.dataset_w.createDimension('scan', dims['scan'])
        self.dataset_w.createDimension('wavelength', dims['wavelength'])

        ##rgb variables
        print(f'[INFO] Creating image variables...')
        rgb_variables = {
            'pictures_sky_irr_1': {
                'ref': '003',
                'oza': 180,
                'oaa': 90,
                'prefix': 'HYPERNETS_W_VEIT_IMG',
                'suffix': '003_180_90_v2.0.jpg'
            },
            'pictures_sky_rad_1': {
                'ref': '006',
                'oza': 140,
                'oaa': 90,
                'prefix': 'HYPERNETS_W_VEIT_IMG',
                'suffix': '006_140_90_v2.0.jpg'
            },
            'pictures_water_rad': {
                'ref': '009',
                'oza': 40,
                'oaa': 90,
                'prefix': 'HYPERNETS_W_VEIT_IMG',
                'suffix': '009_40_90_v2.0.jpg'
            },
            'pictures_sky_rad_2': {
                'ref': '012',
                'oza': 140,
                'oaa': 90,
                'prefix': 'HYPERNETS_W_VEIT_IMG',
                'suffix': '012_140_90_v2.0.jpg'
            },
            'pictures_sky_irr_2': {
                'ref': '015',
                'oza': 180,
                'oaa': 90,
                'prefix': 'HYPERNETS_W_VEIT_IMG',
                'suffix': '015_180_90_v2.0.jpg'
            },
            'pictures_sun': {
                'ref': '016',
                'oza': 0,
                'oaa': 0,
                'prefix': 'HYPERNETS_W_VEIT_IMG',
                'suffix': '016_0_0_v2.0.jpg'

            }
        }
        for rgb_var in rgb_variables:
            var = self.dataset_w.createVariable(rgb_var, 'f8', ('series',), zlib=True, complevel=6)
            for at in rgb_variables[rgb_var]:
                var.setncattr(at, rgb_variables[rgb_var][at])

        ##level1 and level 2 variables
        print(f'[INFO] Creating level 1 variables...')
        self.create_variables(1, seq_list[index_seq_ref])
        print(f'[INFO] Creating level 2 variables...')
        self.create_variables(2, seq_list[index_seq_ref])

        ##var sequence time
        self.dataset_w.createVariable('sequence_ref', 'f8', ('series',), zlib=True, complevel=6)

        return nseq_valid

    def set_data(self, site, date_here):
        self.set_netcdf_data(1)
        self.set_netcdf_data(2)
        self.set_rgb_images_data()
        self.set_sequence_data()
        self.set_global_attributtes(site, date_here)

    def set_sequence_data(self):
        print(f'[INFO] Set sequence reference data...')

        seq_list = list(self.files_dates.keys())
        seq_list.sort()
        for idx in range(len(seq_list)):
            seq = seq_list[idx]
            if not self.files_dates[seq]['valid']:
                continue
            seq_time = dt.strptime(seq, '%Y%m%dT%H%M').replace(tzinfo=pytz.UTC)
            seq_time_stamp = float(seq_time.timestamp())
            self.dataset_w.variables['sequence_ref'][idx] = seq_time_stamp

    def set_global_attributtes(self, site, date_here):
        print(f'[INFO] Set global attributes...')
        self.dataset_w.n_series = len(list(self.files_dates.keys()))
        self.dataset_w.creation_time = dt.utcnow().strftime('%Y-%m-%d %H:%M:%S')
        if date_here is not None:
            self.dataset_w.date = date_here.strftime('%Y-%m-%d')
        if site is not None:
            self.dataset_w.site = site

    def set_netcdf_data(self, level):

        seq_list = list(self.files_dates.keys())
        seq_list.sort()
        for idx in range(len(seq_list)):
            seq = seq_list[idx]
            if not self.files_dates[seq]['valid']:
                continue
            if level == 1:
                file = self.files_dates[seq]['file_l1']
                prename = 'l1'
            if level == 2:
                file = self.files_dates[seq]['file_l2']
                prename = 'l2'
            if file is None:
                continue
            print(f'[INFO] Set level{level} data for sequence {seq}')
            dataset = Dataset(file)
            for var_name in dataset.variables:
                if var_name.startswith('u_rel') or var_name.startswith('err'):
                    continue
                if var_name == 'wavelength' or var_name == 'bandwidth':
                    continue
                var_name_new = f'{prename}_{var_name}'
                dimensions = self.dataset_w.variables[var_name_new].dimensions
                ndim = len(dimensions)
                if level == 1:
                    if ndim == 2:
                        self.dataset_w.variables[var_name_new][idx, :] = dataset.variables[var_name][:]
                    elif ndim == 3:
                        self.dataset_w.variables[var_name_new][idx, :, :] = dataset.variables[var_name][:, :]
                if level == 2:
                    if ndim == 1 and dimensions[0] == 'series':
                        self.dataset_w.variables[var_name_new][idx] = dataset.variables[var_name][0]
                    if ndim == 2 and dimensions[0] == 'series':
                        self.dataset_w.variables[var_name_new][idx, :] = dataset.variables[var_name][:, 0]

            dataset.close()

    def set_rgb_images_data(self):
        seq_list = list(self.files_dates.keys())
        seq_list.sort()
        rgb_variables = {
            '003': {'name_var': 'pictures_sky_irr_1', 'check_at': False},
            '006': {'name_var': 'pictures_sky_rad_1', 'check_at': False},
            '009': {'name_var': 'pictures_water_rad', 'check_at': False},
            '012': {'name_var': 'pictures_sky_rad_2', 'check_at': False},
            '015': {'name_var': 'pictures_sky_irr_2', 'check_at': False},
            '016': {'name_var': 'pictures_sun', 'check_at': False}
        }
        for idx in range(len(seq_list)):
            seq = seq_list[idx]
            if not self.files_dates[seq]['valid']:
                continue
            if self.files_dates[seq]['file_images'] is None:
                continue
            print(f'[INFO] Saving RGB images for sequence: {seq}')
            for file_img in self.files_dates[seq]['file_images']:
                name = file_img.split('/')[-1]
                name_s = name.split('_')
                ref = name_s[6]
                if ref not in rgb_variables.keys():
                    continue
                var_name = rgb_variables[ref]['name_var']
                variable = self.dataset_w.variables[var_name]

                time_here = dt.strptime(name_s[5], '%Y%m%dT%H%M').replace(tzinfo=pytz.UTC)
                time_stamp = float(time_here.timestamp())
                variable[idx] = time_stamp
                if not rgb_variables[ref]['check_at']:  ##check attributes
                    variable.oza = int(name_s[7])
                    variable.oaa = int(name_s[8])
                    variable.prefix = '_'.join(name_s[0:4])
                    variable.suffix = '_'.join(name_s[6:])
                    rgb_variables[ref]['check_at'] = True

    def close_datafile_complete(self):
        self.dataset_w.close()
        print(f'[INFO] Completed')

    def check_dimensions(self, seq):
        file_l1 = self.files_dates[seq]['file_l1']
        if file_l1 is None:
            return None
        dataset = Dataset(file_l1)
        dim_out = {
            'scan': dataset.dimensions['scan'].size,
            'wavelength': dataset.dimensions['wavelength'].size
        }
        dataset.close()

        return dim_out

    def create_variables(self, level, seq):
        if level == 1:
            file = self.files_dates[seq]['file_l1']
            prename = 'l1'
        if level == 2:
            file = self.files_dates[seq]['file_l2']
            prename = 'l2'
        dataset = Dataset(file)
        for var_name in dataset.variables:
            if var_name.startswith('u_rel') or var_name.startswith('err'):
                continue
            if var_name == 'wavelength' or var_name == 'bandwidth':
                if level == 1:
                    continue
                elif level == 2:
                    var_name_new = var_name
            else:
                var_name_new = f'{prename}_{var_name}'
            dimensions = dataset.variables[var_name].dimensions
            if level == 1:
                dimensions = tuple(['series'] + list(dimensions))
            if level == 2 and len(dimensions) == 2 and dimensions[0] == 'wavelength' and dimensions[1] == 'series':
                dimensions = ('series', 'wavelength')
            fillValue = None
            if '_FillValue' in dataset.variables[var_name].ncattrs():
                fillValue = dataset.variables[var_name]._FillValue
            var_new = self.dataset_w.createVariable(var_name_new, dataset.variables[var_name].datatype, dimensions,
                                                    zlib=True, complevel=6, fill_value=fillValue)
            for at in dataset.variables[var_name].ncattrs():
                if at == '_FillValue':
                    continue
                var_new.setncattr(at, dataset.variables[var_name].getncattr(at))
            if (var_name == 'wavelength' or var_name == 'bandwidth') and level == 2:
                var_new[:] = dataset.variables[var_name][:]
        dataset.close()

    def transfer_file_ssh(self, input_path, output_path):
        cmd = f'{self.rsync_url}:{input_path} {output_path}'
        prog = subprocess.Popen(cmd, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
        out, err = prog.communicate()
        if err:
            print(f'[ERROR]{err}')

    def transfer_file_ssh_npl(self, input_path, output_path):
        cmd = f'{self.rsync_url_npl}:{input_path} {output_path}'
        prog = subprocess.Popen(cmd, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
        out, err = prog.communicate()
        if err:
            print(f'[ERROR]{err}')

    def get_sequence_folder_from_list(self, list_sequences, prefix):
        sequence = None
        for seq in list_sequences:
            if seq.startswith(prefix):
                sequence = seq
        return sequence

    def get_folder_date(self, site, date_here):
        folder_date = os.path.join(self.path_data, site, date_here.strftime('%Y'), date_here.strftime('%m'),
                                   date_here.strftime('%d'))
        if os.path.exists(folder_date):
            return folder_date
        else:
            return None

    def get_output_folder_date(self, site, date_here):
        folder_site = os.path.join(self.path_output, site)
        folder_year = os.path.join(folder_site, date_here.strftime('%Y'))
        folder_month = os.path.join(folder_year, date_here.strftime('%m'))
        folder_date = os.path.join(folder_month, date_here.strftime('%d'))
        self.create_if_not_exists(folder_site)
        self.create_if_not_exists(folder_year)
        self.create_if_not_exists(folder_month)
        self.create_if_not_exists(folder_date)

        if os.path.exists(folder_date):
            return folder_date
        else:
            return None

    def create_if_not_exists(self, folder):
        if not os.path.exists(folder):
            try:
                os.mkdir(folder)
                os.chmod(folder,0o775)
            except:
                print(f'[ERROR] Error creating folder: {folder}')


    def get_ssh_path(self, site, date_here, sequence_folder):
        year_str = date_here.strftime('%Y')
        month_str = date_here.strftime('%m')
        day_str = date_here.strftime('%d')
        path = f'{self.base_folder}{site}/{year_str}/{month_str}/{day_str}'
        if sequence_folder is not None:
            path = f'{path}/{sequence_folder}'
        return path

    def get_ssh_path_l2(self, site, sequence_folder):
        base_folder_l2 = self.base_folder_l2_rbins
        if site == 'JSIT':
            base_folder_l2 = self.base_folder_l2_npl
        path = f'{base_folder_l2}{site}/DATA/{sequence_folder}'
        return path

    def get_sequences_date(self, site, date_here):
        ssh_path = self.get_ssh_path(site, date_here, None)
        cmd = f'{self.ssh_base} {self.url_base} ls {ssh_path}'
        list_sequences = self.get_list_files_from_ls_cmd(cmd)
        return list_sequences

    def get_sequences_date_l2(self, site, date_here):
        base_folder_l2 = self.base_folder_l2_rbins
        url_base = self.url_base
        ssh_base = self.ssh_base
        if site == 'JSIT':
            base_folder_l2 = self.base_folder_l2_npl
            url_base = self.url_base_npl
            ssh_base = self.ssh_base_npl
        path_search = f'{base_folder_l2}{site}/DATA/SEQ{date_here.strftime("%Y%m%d")}*'

        cmd = f'{ssh_base} {url_base} ls -d {path_search}'

        list_sequences = self.get_list_files_from_ls_cmd(cmd)
        list_sequences = [x.split('/')[-1] for x in list_sequences]

        return list_sequences

    def get_sequences_date_from_file_list(self, site, date_here):
        folder_date = self.get_folder_date(site, date_here)
        file_list = os.path.join(folder_date, 'sequence_list.txt')
        list_sequences = []
        if os.path.exists(file_list):
            f1 = open(file_list, 'r')
            for line in f1:
                if len(line) > 0:
                    list_sequences.append(line.strip())
            f1.close()
            list_sequences.sort()
        return list_sequences

    def get_images_sequence(self, site, date_here, sequence_folder):
        ssh_path = self.get_ssh_path(site, date_here, sequence_folder)
        path_image = f'{ssh_path}/image'
        cmd = f'{self.ssh_base} {self.url_base} ls {path_image}'
        list_images = self.get_list_files_from_ls_cmd(cmd)
        return list_images

    def get_sun_image_sequence(self, site, sequence_folder):
        url_base = self.url_base
        ssh_base = self.ssh_base
        if site == 'JSIT':
            url_base = self.url_base_npl
            ssh_base = self.ssh_base_npl
        ssh_path = self.get_ssh_path_l2(site, sequence_folder)

        path_image = f'{ssh_path}/RADIOMETER/*_0000_0_0000.jpg'

        cmd = f'{ssh_base} {url_base} ls {path_image}'

        list_images = self.get_list_files_from_ls_cmd(cmd)
        return list_images

    def get_list_files_from_ls_cmd(self, cmd):
        prog = subprocess.Popen(cmd, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
        out, err = prog.communicate()
        list = out.decode('utf-8').split('\n')
        listd = []
        for l in list:
            if l == '':
                continue
            try:
                listd.append(l)
            except:
                pass
        return listd
