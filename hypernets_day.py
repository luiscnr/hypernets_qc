from datetime import datetime as dt
from datetime import timedelta
import os, subprocess, pytz, configparser
from netCDF4 import Dataset
from plot_multiple import PlotMultiple
from hypernets_day_file import HYPERNETS_DAY_FILE
from hypernets_day_file_land import HYPERNETS_DAY_FILE_LAND


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
        self.format_img = '.png'

        print(f'[INFO] Started HYPERNETS_DAY with:')
        print(f'[INFO] ->Input path data: {self.path_data}')
        print(f'[INFO] ->Output path data: {self.path_output}')

        self.url_base = f'hypstar@enhydra.naturalsciences.be'
        self.url_base_npl = f'cnr@hypernetssvr1.npl.co.uk'
        self.base_folder = '/waterhypernet/hypstar/processed_v2/'
        self.base_folder_raw_npl = '/home/cnr/'
        self.base_folder_raw_rbins = '/waterhypernet/HYPSTAR/Raw'  # KGR

        self.ssh_base = 'ssh -X -Y -p 9022'
        self.ssh_base_npl = 'ssh'
        self.rsync_base = f'rsync -a -e \'ssh -p 9022\' {self.url_base}:{self.base_folder}'
        self.rsync_url = f'rsync -a -e \'ssh -p 9022\' {self.url_base}'
        self.rsync_url_npl = f'rsync -a -e \'ssh\' {self.url_base_npl}'

        self.n_series = 50
        self.files_dates = {}
        self.dataset_w = None

        self.rgb_refs = ['003', '006', '009', '012', '015', '016']
        self.rgb_oza = [180, 140, 40, 140, 180, 0]
        self.rgb_oaa = [90, 90, 90, 90, 90, 0]
        self.rgb_pictures_names = ['sky_irr_1', 'sky_rad_1', 'water_rad', 'sky_rad_2', 'sky_irr_2', 'sun']

        self.raw_folder_organization = 'YYYYMMDD'

    def get_files_date_nodownload(self, site, date_here):

        self.files_dates = {}
        date_folder = self.get_folder_date(site, date_here)
        if date_folder is None:
            return

        ##list_sequences is obtained from sequence_list.txt (CNR) or sequence folders (RBINS)
        use_seq_folders = self.check_use_seq_folders(site, date_here)
        list_sequences = self.get_sequences_date_from_file_list(site, date_here)

        if len(list_sequences) == 0:
            print(f'[WARNING] No sequences found for date: {date_here}')
            return
        list_seq_refs = [x[3:-2] for x in list_sequences]

        if use_seq_folders:  ##rbins, files organized in sequence folders
            folders_to_check = [os.path.join(date_folder, seq) for seq in list_sequences]
            folders_to_check_images = [os.path.join(x, 'image') for x in folders_to_check]
        else:  ##cnr, files in the same folder
            folders_to_check = [date_folder]
            folders_to_check_images = [date_folder]

        ##check nc files
        for folder_to_check in folders_to_check:
            if not os.path.exists(folder_to_check):
                print(f'[WARNING] Folder {folder_to_check} does not exist. Skipping...')
                continue
            for name in os.listdir(folder_to_check):
                if name.find('L2A_REF') > 0 and name.endswith('nc'):
                    sequence_ref = name.split('_')[5]
                    try:
                        list_seq_refs.index(sequence_ref)
                        if sequence_ref not in self.files_dates.keys():
                            self.files_dates[sequence_ref] = {
                                'file_l2': os.path.join(folder_to_check, name),
                                'file_l1': None,
                                'file_images': None,
                                'valid': True
                            }

                        else:
                            self.files_dates[sequence_ref]['file_l2'] = os.path.join(folder_to_check, name)
                    except:
                        pass
                if name.find('L1C_ALL') > 0 and name.endswith('nc'):
                    sequence_ref = name.split('_')[5]
                    try:
                        list_seq_refs.index(sequence_ref)
                        if sequence_ref not in self.files_dates.keys():
                            self.files_dates[sequence_ref] = {
                                'file_l2': None,
                                'file_l1': os.path.join(folder_to_check, name),
                                'file_images': None,
                                'valid': True
                            }
                        else:
                            self.files_dates[sequence_ref]['file_l1'] = os.path.join(folder_to_check, name)
                    except:
                        pass

        ##check pictures
        for folder_to_check in folders_to_check_images:
            if not os.path.exists(folder_to_check):
                print(f'[WARNING] Folder {folder_to_check} does not exist. Skipping...')
                continue
            for name in os.listdir(folder_to_check):
                if name.find('IMG') > 0 and name.endswith('jpg'):
                    sequence_ref = name.split('_')[4]
                    try:
                        list_seq_refs.index(sequence_ref)
                        if sequence_ref not in self.files_dates.keys():
                            self.files_dates[sequence_ref] = {
                                'file_l2': None,
                                'file_l1': None,
                                'file_images': [os.path.join(folder_to_check, name)],
                                'valid': True
                            }
                        else:
                            file_images = self.files_dates[sequence_ref]['file_images']
                            if file_images is None:
                                file_images = [os.path.join(folder_to_check, name)]
                            else:
                                file_images.append(os.path.join(folder_to_check, name))
                            self.files_dates[sequence_ref]['file_images'] = file_images
                    except:
                        pass

    def get_disk_usage_log_file(self, site, ndw):
        file_log = os.path.join(self.path_data, site, f'disk-usage_{site}.log')  ##CNR Implementation
        if not os.path.exists(file_log):  ##RBINS Implementation
            if site == 'JSIT':
                file_log = f'{self.base_folder_raw_npl}/{site}/LOGS/disk-usage.log'
            else:
                file_log = f'{self.base_folder_raw_rbins}/{site}/LOGS/disk-usage.log'  ##KGR
        print("file_log=", file_log)
        if ndw:
            if os.path.exists(file_log):
                return file_log
            else:
                return None
        ##ndw (no download) option is only false for testing#############################################
        if site == 'JSIT':
            path_log = f'{self.base_folder_raw_npl}/{site}/LOGS/disk-usage.log'
            self.transfer_file_ssh_npl(path_log, file_log)
        else:
            path_log = f'{self.base_folder_raw_rbins}/{site}/LOGS/disk-usage.log'
            self.transfer_file_ssh(path_log, file_log)
        if os.path.exists(file_log):
            return file_log
        else:
            return None

    def get_last_available_log(self, site, type_log, ndw):
        file_log = os.path.join(self.path_data, site, f'last_{type_log}_{site}.log')  ##CNR IMPLEMENTATION
        if not os.path.exists(file_log):  ##RBINS IMPLEMENTATION
            list_files = self.get_list_log_files_date(dt.now(), site, type_log)  # check same day
            if len(list_files) == 0:  ##check day before
                list_files = self.get_list_log_files_date(dt.now() - timedelta(hours=24), site, type_log)
            if len(list_files) > 0:
                file_log = self.get_last_log_file_from_list(list_files)

        if ndw:  ##OPERATIONAL IMPLEMENENTATION
            if os.path.exists(file_log):
                return file_log
            else:
                return None

        ##Testing uisng ndw==False, to download file logs using ssh (only CNR IMPLEMENTATION)
        list_files = self.get_list_log_files_date(dt.now(), site, type_log)
        if len(list_files) == 0: list_files = self.get_list_log_files_date(dt.now() - timedelta(hours=24), site,
                                                                           type_log)
        if len(list_files) == 0:
            print(f'[WARNING]{type_log} log for {site} was not found in the last 2 days')
            return None
        file_server_last = self.get_last_log_file_from_list(list_files)
        if file_server_last is not None:
            if site == 'JSIT':
                self.transfer_file_ssh_npl(file_server_last, file_log)
            else:
                self.transfer_file_ssh(file_server_last, file_log)
        if os.path.exists(file_log):
            return file_log
        else:
            return None

    def get_last_log_file_from_list(self, list_files):
        date_ref = dt(2020, 1, 1)
        file_server_last = None
        for name in list_files:
            try:
                date_file = dt.strptime(name.split('/')[-1][:15], '%Y-%m-%d-%H%M')
                if date_file > date_ref:
                    file_server_last = name
                    date_ref = date_file
            except:
                pass
        return file_server_last

    # type_log: hello, access, webcam, sequence
    def get_list_log_files_date(self, work_date, site, type_log):
        url_base = self.url_base
        ssh_base = self.ssh_base
        base_folder = self.base_folder_raw_rbins
        if site == 'JSIT':
            url_base = self.url_base_npl
            ssh_base = self.ssh_base_npl
            base_folder = self.base_folder_raw_npl
        work_date_str = work_date.strftime('%Y-%m-%d')
        log_folder = os.path.join(base_folder, site, 'LOGS')
        path_log = f'{base_folder}/{site}/LOGS/{work_date_str}*{type_log}.log'
        if os.path.isdir(log_folder):
            cmd = f'ls {path_log}'
        else:
            cmd = f'{ssh_base} {url_base} ls {path_log}'
        list_files = self.get_list_files_from_ls_cmd(cmd)
        return sorted(list_files)

    def get_sun_images_date(self, site, date_here, ndw):
        sun_images = {}

        date_folder = self.get_folder_date(site, date_here)
        if date_folder is not None:
            for name in os.listdir(date_folder):
                if name.find('_0_0') > 0 and name.endswith('.jpg'):
                    seq = f'{name.split("_")[4]}00'
                    sun_images[seq] = os.path.join(date_folder, name)
                ##checking sun images also in is SEQ*FOLDERS
                if name.startswith('SEQ'):
                    folder_seq = os.path.join(date_folder, name, 'image')
                    if os.path.isdir(folder_seq):
                        for name_s in os.listdir(folder_seq):
                            if name_s.find('_0_0') > 0 and name_s.endswith('.jpg'):
                                seq = f'{name_s.split("_")[4]}00'
                                sun_images[seq] = os.path.join(folder_seq, name)
        if ndw:
            return sun_images

        ##download via ssh, only relevant for CNR implementation

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
            list_sequences = self.get_sequences_date_raw(site, date_here)
        else:
            list_sequences = self.get_sequences_date(site, date_here)
            if len(list_sequences) == 0:
                list_sequences = self.get_sequences_date_raw(site, date_here)

        if len(list_sequences) == 0:
            print(f'[WARNING] No sequences found for date: {date_here}')
            return

        date_folder_raw = None
        if self.raw_folder_organization == 'YYYYMMDD':
            date_folder_raw = f'{date_here.strftime("%Y")}/{date_here.strftime("%m")}/{date_here.strftime("%d")}'

        for seq in list_sequences:
            print(f'[INFO] Checking sun picture in sequence: {seq}')
            # seq_time =  dt.strptime(seq[3:], '%Y%m%dT%H%M%S').replace(tzinfo=pytz.UTC)
            files_sun = self.get_sun_image_sequence(site, date_folder_raw, seq)

            if len(files_sun) == 1:
                name_file = files_sun[0].split('/')[-1]
                file_sun = os.path.join(date_folder, name_file)
                if not os.path.exists(file_sun):
                    if site == 'JSIT':
                        self.transfer_file_ssh_npl(files_sun[0], file_sun)
                    else:
                        self.transfer_file_ssh(files_sun[0], file_sun)
                if os.path.exists(file_sun):
                    name_new = self.get_name_new_file_sun(file_sun, site, seq)
                    file_sun_new = os.path.join(date_folder, name_new)
                    os.rename(file_sun, file_sun_new)
                    sun_images[seq[3:]] = file_sun_new
                else:
                    sun_images[seq[3:]] = None

        return sun_images

    def get_name_new_file_sun(self, file_img, site, seq):
        type = 'W'
        picture = '016'
        if site == 'JSIT':
            type = 'L'
            picture = '090'
        date_img = dt.fromtimestamp(os.path.getmtime(file_img))
        date_img_str = date_img.strftime('%Y%m%dT%H%M')
        name_new = f'HYPTERNETS_{type}_{site}_IMG_{seq[3:-2]}_{date_img_str}_{picture}_0_0_v2.0.jpg'
        return name_new

    def get_sequences_info(self, site, date_here, sequences_with_data, sequences_abs_range):

        all_sequences = self.get_sequences_date_from_file_list(site, date_here)

        all_sequences = [x[:-2] for x in all_sequences]

        if sequences_abs_range is not None:
            all_sequences_orig = all_sequences.copy()
            all_sequences = []
            for seq in all_sequences_orig:
                date_here = dt.strptime(seq[3:], '%Y%m%dT%H%M').replace(tzinfo=pytz.utc).timestamp()
                if sequences_abs_range[0] <= date_here <= sequences_abs_range[1]:
                    all_sequences.append(seq)

        sequences_with_data = [f'SEQ{x}' for x in sequences_with_data]

        sequences_without_data = []
        all_sequences_info = {}
        for seq in all_sequences:
            if seq not in sequences_with_data:
                sequences_without_data.append(seq)
                all_sequences_info[seq] = -1
            else:
                all_sequences_info[seq] = sequences_with_data.index(seq)

        return sequences_without_data, all_sequences_info

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

    def save_sun_images(self, output_file, sun_images, time_list):
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
        pm.save_fig_with_resolution(output_file, 150)
        pm.close_plot()

    def get_file_date_complete(self, site, date_here):
        folder_date = self.get_output_folder_date(site, date_here)
        if folder_date is None:
            return None
        date_here_str = date_here.strftime('%Y%m%d')
        file_date = os.path.join(folder_date, f'HYPERNETES_W_DAY_{date_here_str}.nc')

        return file_date

    def get_file_date_land_complete(self, site, date_here):
        folder_date = self.get_output_folder_date(site, date_here)
        if folder_date is None:
            return None
        date_here_str = date_here.strftime('%Y%m%d')
        file_date = os.path.join(folder_date, f'HYPERNETES_L_DAY_{date_here_str}.nc')

        return file_date

    def get_files_img_for_sequences_no_data(self, site, date_here, seq, use_seq_folders):
        files_img = {}
        for name_img, ref in zip(self.rgb_pictures_names, self.rgb_refs):
            files_img[ref] = {
                'name_img': name_img,
                'file_img': None
            }
        date_folder = self.get_folder_date(site, date_here)
        if date_folder is None:
            return

        if use_seq_folders:
            for name in os.listdir(date_folder):
                if name.startswith(seq):
                    date_folder = os.path.join(date_folder, name, 'image')
                    break

        seq_ref = seq.replace('SEQ', 'IMG_')
        for name in os.listdir(date_folder):
            if name.endswith('.jpg') and name.find(seq_ref) > 0:
                name_s = name.split('_')
                ref = name_s[6]
                if ref not in files_img.keys():
                    continue
                files_img[ref]['file_img'] = os.path.join(date_folder, name)
        return files_img

    def get_sequence_range(self, date_here, config_file_summary, absolute):

        # from datetime import datetime as dt
        options = configparser.ConfigParser()
        options.read(config_file_summary)

        if options.has_option('sequence_info', 'start_time') and options.has_option('sequence_info', 'end_time'):
            start_time_str = options['sequence_info']['start_time']
            end_time_str = options['sequence_info']['end_time']
            frequency = 30
            if options.has_option('sequence_info', 'frequency'):
                try:
                    frequency = int(options['sequence_info']['frequency'])
                except:
                    pass

            try:
                start_time = dt.strptime(f'{date_here.strftime("%Y-%m-%d")}T{start_time_str}', '%Y-%m-%dT%H:%M')
                end_time = dt.strptime(f'{date_here.strftime("%Y-%m-%d")}T{end_time_str}', '%Y-%m-%dT%H:%M')
                frequency_seconds = frequency * 60
                nsequences = ((end_time.timestamp() - start_time.timestamp()) / frequency_seconds) + 1
                if absolute:
                    end_time = end_time + timedelta(minutes=frequency)
                range = [start_time.replace(tzinfo=pytz.UTC).timestamp(), end_time.replace(tzinfo=pytz.UTC).timestamp(),
                         int(nsequences)]
                return range
            except:
                return None

        return None

    def get_hypernets_day_file(self, site, date_here):

        if site == 'JSIT':
            file_date = self.get_file_date_land_complete(site, date_here)
        else:
            file_date = self.get_file_date_complete(site, date_here)
        if file_date is None:
            return None
        if os.path.exists(file_date):
            if site == 'JSIT':
                return HYPERNETS_DAY_FILE_LAND(file_date, self.path_data)
            else:
                return HYPERNETS_DAY_FILE(file_date, self.path_data)
        else:
            print(f'[WARNING] Expected HYPERNETS day file {file_date} does not exist')
            return None

    def start_file_date_land_complete(self, site, date_here, overwrite):

        file_date = self.get_file_date_land_complete(site, date_here)
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

            dims_here = self.check_dimensions_land(seq)
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
        self.dataset_w.createDimension('sequence')
        self.dataset_w.createDimension('series', dims['series'])
        self.dataset_w.createDimension('wavelength', dims['wavelength'])

        rgb_variables = {}
        for idx in range(6):
            rgb_variables[f'pictures_{self.rgb_pictures_names[idx]}'] = {
                'ref': self.rgb_refs[idx],
                'oza': self.rgb_oza[idx],
                'oaa': self.rgb_oaa[idx],
                'prefix': f'HYPERNETS_L_{site}_IMG',
                'suffix': f'{self.rgb_refs[idx]}_{self.rgb_oza[idx]}_{self.rgb_oaa[idx]}_v2_0.jpg'
            }
        for rgb_var in rgb_variables:
            var = self.dataset_w.createVariable(rgb_var, 'f8', ('sequence',), zlib=True, complevel=6)
            for at in rgb_variables[rgb_var]:
                var.setncattr(at, rgb_variables[rgb_var][at])

        ##level1 and level 2 variables
        print(f'[INFO] Creating level 1 variables...')
        self.create_variables_land(1, seq_list[index_seq_ref])
        print(f'[INFO] Creating level 2 variables...')
        self.create_variables_land(2, seq_list[index_seq_ref])

        ##var sequence time
        self.dataset_w.createVariable('sequence_ref', 'f8', ('sequence',), zlib=True, complevel=6)

        return nseq_valid

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
        # self.rgb_refs = ['003', '006', '009', '012', '015', '016']
        rgb_variables = {}
        for idx in range(6):
            rgb_variables[f'pictures_{self.rgb_pictures_names[idx]}'] = {
                'ref': self.rgb_refs[idx],
                'oza': self.rgb_oza[idx],
                'oaa': self.rgb_oaa[idx],
                'prefix': f'HYPERNETS_W_{site}_IMG',
                'suffix': f'{self.rgb_refs[idx]}_{self.rgb_oza[idx]}_{self.rgb_oaa[idx]}_v2_0.jpg'
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

    def set_data_land(self, site, date_here):
        self.set_netcdf_data_land(1)
        self.set_netcdf_data_land(2)
        self.set_rgb_images_data()
        self.set_sequence_data()
        self.set_global_attributtes(site, date_here)

    def set_sequence_data(self):
        print(f'[INFO] Set sequence reference data...')

        seq_list = list(self.files_dates.keys())
        seq_list.sort()
        index_add = -1
        for idx in range(len(seq_list)):
            seq = seq_list[idx]
            if not self.files_dates[seq]['valid']:
                continue
            index_add = index_add + 1
            seq_time = dt.strptime(seq, '%Y%m%dT%H%M').replace(tzinfo=pytz.UTC)
            seq_time_stamp = float(seq_time.timestamp())
            self.dataset_w.variables['sequence_ref'][index_add] = seq_time_stamp

    def set_global_attributtes(self, site, date_here):
        print(f'[INFO] Set global attributes...')
        self.dataset_w.n_sequences = len(list(self.files_dates.keys()))
        self.dataset_w.creation_time = dt.utcnow().strftime('%Y-%m-%d %H:%M:%S')
        if date_here is not None:
            self.dataset_w.date = date_here.strftime('%Y-%m-%d')
        if site is not None:
            self.dataset_w.site = site

    def set_netcdf_data(self, level):

        seq_list = list(self.files_dates.keys())
        seq_list.sort()
        index_add = -1
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
            index_add = index_add + 1
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
                        self.dataset_w.variables[var_name_new][index_add, :] = dataset.variables[var_name][:]
                    elif ndim == 3:
                        self.dataset_w.variables[var_name_new][index_add, :, :] = dataset.variables[var_name][:, :]
                if level == 2:
                    if ndim == 1 and dimensions[0] == 'series':
                        self.dataset_w.variables[var_name_new][index_add] = dataset.variables[var_name][0]
                    if ndim == 2 and dimensions[0] == 'series':
                        self.dataset_w.variables[var_name_new][index_add, :] = dataset.variables[var_name][:, 0]

            dataset.close()

    def set_netcdf_data_land(self, level):

        seq_list = list(self.files_dates.keys())
        seq_list.sort()
        index_add = -1
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

            index_add = index_add + 1
            print(f'[INFO] Set level{level} data for sequence {seq} [{index_add}]')
            dataset = Dataset(file)
            for var_name in dataset.variables:
                if var_name.startswith('u_rel') or var_name.startswith('err'):
                    continue
                if var_name == 'wavelength' or var_name == 'bandwidth':
                    continue
                var_name_new = f'{prename}_{var_name}'
                dimensions = self.dataset_w.variables[var_name_new].dimensions
                ndim = len(dimensions)

                if ndim == 2:
                    self.dataset_w.variables[var_name_new][index_add, :] = dataset.variables[var_name][:]
                elif ndim == 3:
                    self.dataset_w.variables[var_name_new][index_add, :, :] = dataset.variables[var_name][:].transpose()

            dataset.close()

    def set_raw_data_organization(self, config_file_summary):
        options = configparser.ConfigParser()
        options.read(config_file_summary)
        if options.has_option('GLOBAL_OPTIONS', 'raw_data_organization'):
            values = ['NONE', 'YYYYMMDD']
            new_value = options['GLOBAL_OPTIONS']['raw_data_organization'].strip().upper()
            if new_value in values:
                self.raw_folder_organization = new_value
            else:
                print(
                    f'[WARNING] Option {new_value} in not valid for raw_folder_organization. Plese use of the following values {values}')
                print(f'[WARNING] Default value {self.raw_folder_organization} will be used')

    def set_rgb_refs(self, config_file_summary):
        options = configparser.ConfigParser()
        options.read(config_file_summary)
        if options.has_option('GLOBAL_OPTIONS', 'rgb_refs'):
            value = options['GLOBAL_OPTIONS']['rgb_refs'].strip()
            if len(value.split(',')) == 6:
                for idx, x in enumerate(value.split(',')):
                    xs = x.split('_')
                    if len(xs) == 1:
                        self.rgb_refs[idx] = xs[0].strip()
                    elif len(xs) == 3:
                        self.rgb_refs[idx] = xs[0].strip()
                        self.rgb_oza[idx] = int(xs[1].strip())
                        self.rgb_oaa[idx] = int(xs[2].strip())
                # self.rgb_refs = [x.strip() for x in value.split(',')]
                print(f'[INFO] Camera image refs set to: {self.rgb_refs}')
                print(f'[INFO]   Zenith angles:{self.rgb_oza}')
                print(f'[INFO]   Azimuth angles:{self.rgb_oaa}')

        if options.has_option('GLOBAL_OPTIONS', 'rgb_names'):
            value = options['GLOBAL_OPTIONS']['rgb_names'].strip()
            if len(value.split(',')) == 6:
                self.rgb_pictures_names = [x.strip() for x in value.split(',')]
                print(f'[INFO] Camera images names set to: {self.rgb_pictures_names}')

    def set_rgb_images_data(self):
        seq_list = list(self.files_dates.keys())
        seq_list.sort()
        # self.rgb_refs = ['003', '006', '009', '012', '015', '016']
        rgb_variables = {
            self.rgb_refs[0]: {'name_var': f'pictures_{self.rgb_pictures_names[0]}', 'check_at': False},
            self.rgb_refs[1]: {'name_var': f'pictures_{self.rgb_pictures_names[1]}', 'check_at': False},
            self.rgb_refs[2]: {'name_var': f'pictures_{self.rgb_pictures_names[2]}', 'check_at': False},
            self.rgb_refs[3]: {'name_var': f'pictures_{self.rgb_pictures_names[3]}', 'check_at': False},
            self.rgb_refs[4]: {'name_var': f'pictures_{self.rgb_pictures_names[4]}', 'check_at': False},
            self.rgb_refs[5]: {'name_var': f'pictures_{self.rgb_pictures_names[5]}', 'check_at': False}
        }
        index_add = -1
        for idx in range(len(seq_list)):
            seq = seq_list[idx]
            if not self.files_dates[seq]['valid']:
                continue
            index_add = index_add + 1
            if self.files_dates[seq]['file_images'] is None:
                continue
            print(f'[INFO] Saving RGB images for sequence: {seq} [{index_add}]')
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
                variable[index_add] = time_stamp
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

    def check_dimensions_land(self, seq):
        file_l1 = self.files_dates[seq]['file_l1']
        if file_l1 is None or self.files_dates[seq]['file_l2'] is None:
            return None
        dataset = Dataset(file_l1)
        dim_out = {
            'wavelength': dataset.dimensions['wavelength'].size,
            'series': dataset.dimensions['series'].size
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

    def create_variables_land(self, level, seq):
        from netCDF4 import Dataset
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
            if len(dimensions) == 1 and dimensions[0] == 'series':
                dimensions = tuple(['sequence'] + list(dimensions))
            if len(dimensions) == 2 and dimensions[0] == 'wavelength' and dimensions[1] == 'series':
                dimensions = ('sequence', 'series', 'wavelength')

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
                os.chmod(folder, 0o775)
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

    def get_ssh_path_raw(self, site, date_folder, sequence_folder):
        raw_folder = self.base_folder_raw_rbins
        if site == 'JSIT':
            raw_folder = self.base_folder_raw_npl
        if date_folder is None:
            path = f'{raw_folder}{site}/DATA/{sequence_folder}'
        else:
            path = f'{raw_folder}{site}/DATA/{date_folder}/{sequence_folder}'
        return path

    def get_sequences_date(self, site, date_here):
        ssh_path = self.get_ssh_path(site, date_here, None)
        cmd = f'{self.ssh_base} {self.url_base} ls {ssh_path}'
        list_sequences = self.get_list_files_from_ls_cmd(cmd)
        return sorted(list_sequences)

    def get_sequences_date_raw(self, site, date_here):
        raw_folder = self.base_folder_raw_rbins
        url_base = self.url_base
        ssh_base = self.ssh_base
        if site == 'JSIT':
            raw_folder = self.base_folder_raw_npl
            url_base = self.url_base_npl
            ssh_base = self.ssh_base_npl
        path_search = f'{raw_folder}{site}/DATA/SEQ{date_here.strftime("%Y%m%d")}*'

        cmd = f'{ssh_base} {url_base} ls -d {path_search}'

        list_sequences = self.get_list_files_from_ls_cmd(cmd)
        list_sequences = [x.split('/')[-1] for x in list_sequences]

        return sorted(list_sequences)

    def check_use_seq_folders(self, site, date_here):
        folder_date = self.get_folder_date(site, date_here)
        if folder_date is None:
            return None
        file_list = os.path.join(folder_date, 'sequence_list.txt')
        use_seq_folder = False if os.path.exists(file_list) else True
        return use_seq_folder

    def get_sequences_date_from_file_list(self, site, date_here):
        list_sequences = []
        folder_date = self.get_folder_date(site, date_here)
        if folder_date is None:
            return sorted(list_sequences)
        file_list = os.path.join(folder_date, 'sequence_list.txt')
        if os.path.exists(file_list):  ##CNR
            f1 = open(file_list, 'r')
            for line in f1:
                if len(line) > 0:
                    list_sequences.append(line.strip())
            f1.close()
            list_sequences.sort()
        else:  # RBINS
            for name in os.listdir(folder_date):
                if name.startswith('SEQ') and os.path.isdir(os.path.join(folder_date, name)):
                    list_sequences.append(name)
            list_sequences.sort()

        return sorted(list_sequences)

    def get_images_sequence(self, site, date_here, sequence_folder):
        ssh_path = self.get_ssh_path(site, date_here, sequence_folder)
        path_image = f'{ssh_path}/image'
        cmd = f'{self.ssh_base} {self.url_base} ls {path_image}'
        list_images = self.get_list_files_from_ls_cmd(cmd)
        return list_images

    def get_sun_image_sequence(self, site, date_folder, sequence_folder):
        url_base = self.url_base
        ssh_base = self.ssh_base
        if site == 'JSIT':
            url_base = self.url_base_npl
            ssh_base = self.ssh_base_npl
        ssh_path = self.get_ssh_path_raw(site, date_folder, sequence_folder)

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
        return sorted(listd)

    def save_report_image_only_pictures(self, site, delete_images, overwrite, seq, files_img, input_path_report):
        print(f'[INFO] Sequence {seq} (No Level-2 data available)')
        seq_time_str = seq[3:]
        seq_time = dt.strptime(seq_time_str, '%Y%m%dT%H%M')
        file_out = os.path.join(input_path_report, f'{site}_{seq_time_str}_Report{self.format_img}')
        if os.path.exists(file_out) and not overwrite:
            return
        names_img = ['sky_irr_1', 'sky_rad_1', 'water_rad', 'sky_rad_2', 'sky_irr_2', 'sun']
        names_img_files = {}
        for ref in files_img:
            name_img = files_img[ref]['name_img']
            if name_img in names_img:
                names_img_files[name_img] = files_img[ref]['file_img']

        dir_img = os.path.join(input_path_report, 'IMG')
        if not os.path.exists(dir_img):
            os.mkdir(dir_img)

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
        pm.save_fig(file_out)
        pm.close_plot()

        if delete_images:
            for name in os.listdir(dir_img):
                file_here = os.path.join(dir_img, name)
                os.remove(file_here)
            os.rmdir(dir_img)
