import os,configparser
from datetime import datetime as dt
from datetime import timezone
from datetime import timedelta
from matplotlib import pyplot as plt
from plot_options import PlotOptions
from flag_builder import FlagBuilder
from hypernets_day import HYPERNETS_DAY_BASE
from matplotlib.backends.backend_pdf import PdfPages

class ConfigOptions:
    def __init__(self,config_file):
        try:
            self.options = configparser.ConfigParser()
            self.options.read(config_file)
        except Exception as ex:
            print(f'[ERROR] {config_file} could not be read: {ex}')
            self.options = None

    def get_rgb_refs(self):
        if self.options is None:
            return None

        rgb_refs = []
        rgb_oza = []
        rgb_oaa = []
        rgb_pictures_names = []

        if self.options.has_option('GLOBAL_OPTIONS', 'rgb_refs'):
            value = self.options['GLOBAL_OPTIONS']['rgb_refs'].strip()
            if len(value.split(',')) == 6:
                for idx, x in enumerate(value.split(',')):
                    xs = x.split('_')
                    if len(xs) == 1:
                        rgb_refs.append(xs[0].strip())
                    elif len(xs) == 3:
                        rgb_refs.append(xs[0].strip())
                        rgb_oza.append(xs[1].strip())
                        rgb_oaa.append(xs[2].strip())
                print(f'[INFO] Camera image refs set to: {rgb_refs}')
                print(f'[INFO]   Zenith angles:{rgb_oza}')
                print(f'[INFO]   Azimuth angles:{rgb_oaa}')

        if self.options.has_option('GLOBAL_OPTIONS', 'rgb_names'):
            value = self.options['GLOBAL_OPTIONS']['rgb_names'].strip()
            if len(value.split(',')) == 6:
                rgb_pictures_names = [x.strip() for x in value.split(',')]
            print(f'[INFO] Camera images names set to: {rgb_pictures_names}')

        result = {
            'rgb_refs': rgb_refs if len(rgb_refs) == 6 else None,
            'rgb_oza': rgb_oza if len(rgb_oza) == 6 else None,
            'rgb_oaa': rgb_oaa if len(rgb_oaa) == 6 else None,
            'rgb_pictures_names': rgb_pictures_names if len(rgb_pictures_names) == 6 else None
        }
        return result

    def get_sequence_times_and_frequency(self,date_here):
        frequency = 30
        start_time = None
        end_time = None
        if self.options.has_option('sequence_info', 'start_time') and self.options.has_option('sequence_info', 'end_time'):
            start_time_str = self.options['sequence_info']['start_time']
            end_time_str = self.options['sequence_info']['end_time']
            try:
                start_time = dt.strptime(f'{date_here.strftime("%Y-%m-%d")}T{start_time_str}', '%Y-%m-%dT%H:%M').replace(tzinfo=timezone.utc)
                end_time = dt.strptime(f'{date_here.strftime("%Y-%m-%d")}T{end_time_str}', '%Y-%m-%dT%H:%M').replace(tzinfo=timezone.utc)
            except Exception as ex:
                print(f'[ERROR] Start ({self.options["sequence_info"]["start_time"]}) and/or end time ({self.options["sequence_info"]["end_time"]}) are not in the correct format HH:MM: {ex}')

        if self.options.has_option('sequence_info', 'frequency'):
            try:
                frequency = int(self.options['sequence_info']['frequency'])
            except Exception as ex:
                print(
                    f'[WARNING] Option sequence_info/frequency {self.options["sequence_info"]["frequency"]} is not a valid integer value: {ex}')
                print(f'[WARNING] Using default value: {frequency}')
                pass


        return start_time,end_time,frequency

    def get_sequence_range(self, date_here, absolute):

        start_time, end_time, frequency = self.get_sequence_times_and_frequency(date_here)
        frequency_seconds = frequency * 60
        if start_time is not None and end_time is not None:
            n_sequences = ((end_time.timestamp() - start_time.timestamp()) / frequency_seconds) + 1
            if absolute:
                end_time = end_time + timedelta(minutes=frequency)
            range_here = [start_time.timestamp(), end_time.timestamp(),int(n_sequences)]
            return range_here
        else:
            return None

class DailyMail:

    def __init__(self,file_qc_mail):
        self.fout = open(file_qc_mail, 'w')

    def add_new_line(self, str):
        self.fout.write('\n')
        self.fout.write(str)

    def create_file(self,site,start_date,daily_sequences_summary,extra_info):
        self.fout.write(f'QUALITY CONTROL - {site} - {start_date.strftime("%Y-%m-%d")}')
        self.add_new_line('===================================')
        self.add_new_line('')
        if daily_sequences_summary is not None:
            if daily_sequences_summary['NTotal'] == 0:
                self.add_new_line(f'WARNING: No folder sequences found for {site} on {start_date.strftime("%Y-%m-%d")}')
                self.add_new_line(f'Please review if the system is working and folders are available in the server')
                self.add_new_line('')
            self.add_new_line('SEQUENCES SUMMARY')
            self.add_new_line('=================')
            self.add_new_line(f'Start time:  {daily_sequences_summary["start_time"]}')
            self.add_new_line(f'End time: {daily_sequences_summary["end_time"]}')
            self.add_new_line(f'Expected sequences: {daily_sequences_summary["expected_sequences"]}')
            self.add_new_line(f'Available sequences: {daily_sequences_summary["NTotal"]}')
            self.add_new_line(f'Sequences processed to L2: {daily_sequences_summary["NAvailable"]}')
            if 'VALID' in daily_sequences_summary.keys():
                self.add_new_line(f'Valid sequences after quality control: {daily_sequences_summary["VALID"]}')
            self.add_new_line('')

        self.add_new_line('SYSTEM STATUS')
        self.add_new_line('=============')

        ##disk usage
        file_log_disk_usage = extra_info['file_log_disk_usage']
        if file_log_disk_usage is not None:
            lines_disk_usage = self.get_lines_disk_usage(file_log_disk_usage)
            for line in lines_disk_usage:
                self.add_new_line(line)

        ##last log sequence
        file_log_last_sequence = extra_info['file_log_last_sequence']  ##None, not implemented
        if file_log_last_sequence is not None:
            ##no implemented
            pass

        self.add_new_line('')
        self.add_new_line('DAILY CHECKING FILES')
        self.add_new_line('====================')
        self.add_new_line(f'Output folder: {extra_info["folder_day"]}')
        file_summary = os.path.join(extra_info["folder_day"], extra_info["name_summary"])
        self.add_new_line(f'Summary file: {file_summary if os.path.exists(file_summary) else "Not. Av."}')
        self.add_new_line('')
        self.add_new_line(f'PDF file: {extra_info["file_pdf"] if os.path.exists(extra_info["file_pdf"]) else "Not. Av."}')
        self.add_new_line('')

        self.fout.close()

    def get_lines_disk_usage(self,file_log):
        lines = ['']
        if not os.path.exists(file_log):
            return lines
        df = pd.read_csv(file_log, sep=' ')
        lines.append('DISK USAGE')
        lines.append('----------')
        last_line = df.iloc[-1]
        used = float(last_line.iloc[1]) / (1024 * 1024)
        av = float(last_line.iloc[2]) / (1024 * 1024)
        lines.append(
            f' Last measurement: {last_line.iloc[0]} Used: {used:.2f} Gb. Available: {av:.2f} Gb. %Use: {last_line.iloc[4]}')

        porc_ref = float(str(last_line.iloc[4])[:-1])

        nlines = len(df.index)
        last_five_dates = {}
        date_ref = dt.strptime(last_line.iloc[0][:15], '%Y-%m-%d-%H%M').replace(hour=12, minute=0, second=12)
        date_ref_str_loop = date_ref.strftime('%Y-%m-%d')
        for i in range(5):
            date_ref = date_ref - timedelta(hours=24)
            date_ref_str = date_ref.strftime('%Y-%m-%d')
            last_five_dates[date_ref_str] = None

        first_date_here_str = None
        last_date_here_str = None
        used_array = []
        total_array = []
        porc_use_array = []

        for idx in range(nlines - 1, 0, -1):
            line_here = df.loc[idx]
            if pd.isna(line_here.iloc[0]): continue
            date_here_str = str(line_here.iloc[0])[:10]

            if date_here_str != date_ref_str_loop:
                date_ref_str_loop = date_here_str
                date_here_str_basic = dt.strptime(date_here_str, '%Y-%m-%d').strftime('%Y-%m-%d')
                if date_here_str_basic in last_five_dates.keys():
                    last_five_dates[date_here_str_basic] = line_here
                if last_date_here_str is None:
                    last_date_here_str = str(line_here.iloc[0])

                porc_here = float(str(line_here.iloc[4])[:-1])
                if abs(porc_ref - porc_here) < 2:
                    porc_ref = porc_here
                    used_array.append(float(line_here.iloc[1]))
                    total_array.append(float(line_here.iloc[1]) + float(line_here.iloc[2]))
                    porc_use_array.append(line_here.iloc[4])
                    first_date_here_str = str(line_here.iloc[0])
                else:
                    break

        lines.append(f' Overall period:')
        start_used = used_array[-1] / (1024 * 1024)
        end_used = used_array[0] / (1024 * 1024)
        used_increase = []
        porc_used_increase = []
        for idx in range(1, len(used_array)):
            used_increase.append(used_array[idx - 1] - used_array[idx])
            porc_used_increase.append(float(str(porc_use_array[idx - 1])[:-1]) - float(str(porc_use_array[idx])[:-1]))

        avg_increase = np.mean(np.array(used_increase))
        avg_increase_mb = avg_increase / 1024
        avg_increase_porc = np.mean(np.array(porc_used_increase))

        total_99 = np.min(total_array) * 0.99

        remaining = total_99 - used_array[-1]
        ndays = np.floor(remaining / avg_increase)
        last_day = dt.strptime(last_date_here_str[:10], '%Y-%m-%d')
        day_fill = last_day + timedelta(days=ndays)

        lines.append(f'  Start: {first_date_here_str} Used: {start_used:.2f} Gg. %Used: {porc_use_array[-1]}')
        lines.append(f'  End: {last_date_here_str} Used: {end_used:.2f} Gg. %Used: {porc_use_array[0]}')
        lines.append(f'  Average daily increase: {avg_increase_mb:.2f} Mb. ({avg_increase_porc:.3f}%).')
        lines.append(f'  Number of days until 99%: {ndays:.0f} ({day_fill.strftime("%Y-%m-%d")}).')

        lines.append(' Last five days: ')
        for last_five_date in last_five_dates:
            line_here = last_five_dates[last_five_date]
            if line_here is None:
                lines.append(f'  {last_five_date}: Data are not available')
            else:
                used = float(line_here.iloc[1]) / (1024 * 1024)
                av = float(line_here.iloc[2]) / (1024 * 1024)
                lines.append(
                    f'  {line_here.iloc[0]} Used: {used:.2f} Gb. Available: {av:.2f} Gb. %Use: {line_here.iloc[4]}')
        lines.append('')

        return lines

def plot_from_options(hfile, config_file, output_path_images, sequences_no_data, only_sequences_summary,verbose=True):
    if not os.path.exists(config_file):
        print(f'[ERROR] Plot configuration file: {config_file} does not exist. ')
        return None

    if verbose:
        print(f'[INFO] Started plotting from file: {hfile.file_nc}')

    options = configparser.ConfigParser()
    options.read(config_file)

    poptions = PlotOptions(options, None)
    poptions.set_global_options()

    if output_path_images is not None:
        if not os.path.exists(output_path_images):
            try:
                os.mkdir(output_path_images)
            except:
                pass
        if os.path.exists(output_path_images):
            poptions.global_options['output_path'] = output_path_images

    if poptions.global_options['output_path'] is None:
        poptions.global_options['output_path'] = os.path.dirname(hfile.file_nc)

    fbuilder = FlagBuilder(hfile.file_nc, options)
    hfile.flag_builder = fbuilder

    list_figures = poptions.get_list_figures()

    daily_sequences_summary = None
    for figure in list_figures:
        print('------------------------------------------------------------------------------------------')
        print(f'[INFO] Starting figure: {figure}')
        options_figure = poptions.get_options(figure)
        if options_figure is None:
            continue

        if options_figure['type'] == 'sequence':
            start_time_str = options_figure['start_time']
            end_time_str = options_figure['end_time']
            sequences_no_data_real = []
            for seq in sequences_no_data:
                seq_time = dt.strptime(seq[3:], '%Y%m%dT%H%M')
                seq_min = dt.strptime(f'{seq_time.strftime("%Y%m%d")}T{start_time_str}', '%Y%m%dT%H:%M')
                seq_max = dt.strptime(f'{seq_time.strftime("%Y%m%d")}T{end_time_str}', '%Y%m%dT%H:%M')
                if seq_min <= seq_time <= seq_max:
                    # print(seq)
                    sequences_no_data_real.append(seq)
            hfile.sequences_no_data = sequences_no_data_real

        if options_figure['apply'] and options_figure['type'] == 'sequence':
            daily_sequences_summary = hfile.plot_from_options_impl(options_figure)
        else:
            if not only_sequences_summary:
                hfile.plot_from_options_impl(options_figure)

    return daily_sequences_summary

def get_start_and_end_dates(args):
    start_date = dt.now()
    if args.start_date:
        try:
            start_date = dt.strptime(args.start_date, '%Y-%m-%d')
        except Exception as ex:
            print(f'[ERROR] Start date is not in valid format YYYY-mm-dd. Exception: {ex}')
            return None, None
    if args.end_date:
        try:
            end_date = dt.strptime(args.end_date, '%Y-%m-%d')
        except Exception as ex:
            print(f'[ERROR] End date is not in valid format YYYY-mm-DD. Exception: {ex}')
            return None, None
    else:
        end_date = start_date

    return start_date, end_date

def get_start_and_end_times(args):
    start_time = '00:00'
    end_time = '23:59'
    if args.start_time:
        try:
            start_time_p = dt.strptime(args.start_time, '%H:%M')
            start_time = start_time_p.strftime('%H:%M')
        except Exception as ex:
            print(f'[WARNING] Start time is not in valid format HH:MM. Exception: {ex}]')
            pass
    if args.end_time:
        try:
            end_time_p = dt.strptime(args.end_time, '%H:%M')
            end_time = end_time_p.strftime('%H:%M')
        except Exception as ex:
            print(f'[WARNING] End time is not in valid format HH:MM. Exception: {ex}]')
            pass

    return start_time, end_time

def check_input_path(args):
    if args.input_path:
        if os.path.isdir(args.input_path):
            return args.input_path
        else:
            print(f'[ERROR] Input path {args.input_path} is not a valid directory')
            return None
    else:
        print(f'[ERROR] Argument input path is required')
        return None

def check_output_path(output_path):
    if not os.path.isdir(output_path):
        try:
            os.makedirs(output_path, exist_ok=False)
        except OSError as ex:
            print(f'[ERROR] Output path {output_path} is not a valid directory and could not be created: {ex}')
            return None
    return output_path

def get_path_date(path_base,site,date_here,create_dir=False):
    path_date = os.path.join(path_base, site, date_here.strftime('%Y'), date_here.strftime('%m'),date_here.strftime('%d'))
    if not os.path.isdir(path_date) and create_dir:
        try:
            os.makedirs(path_date)
        except OSError as ex:
            print(f'[ERROR] {path_date} could not be created: {ex}')
            return None
    if os.path.isdir(path_date):
        return path_date
    else:
        return None

def create_empty_image(file_img, site, date_here):
    plt.figure(figsize=(6, 0.75))
    plt.title(f'L2 data were not available for {site} on {date_here.strftime("%Y-%m-%d")}')
    plt.xticks([])
    plt.yticks([])
    plt.axis('off')
    plt.tight_layout()
    plt.savefig(file_img, dpi=300)

def create_daily_pdf_report(input_path, output_path, site, date_here, file_summary, sequences,overwrite=False):
    hday = HYPERNETS_DAY_BASE(input_path, output_path)
    folder_day = hday.get_output_folder_date(site, date_here)
    date_here_str = date_here.strftime('%Y%m%d')
    file_pdf = os.path.join(folder_day, f'Report_{site}_{date_here_str}.pdf')
    print(f'[INFO] PDF report: {file_pdf}')
    if os.path.exists(file_pdf) and not overwrite:
        print(f'[WARNING] PDF report already exists')
        return
    pdf = PdfPages(file_pdf)
    if file_summary is not None and os.path.exists(file_summary):
        plt.close()
        fig = plt.figure(figsize=(10, 18))
        plt.imshow(plt.imread(file_summary))
        plt.axis('off')
        fig.tight_layout()
        pdf.savefig(dpi=300, bbox_inches='tight')
    for sequence in sequences:
        print(f'[INFO] Adding sequence to PDF file: {sequence}')
        if sequence is not None:
            file_img = os.path.join(folder_day, f'{site}_{sequence[3:]}_Report.png')
            if os.path.exists(file_img):
                plt.close()
                fig = plt.figure(figsize=(10, 18))
                plt.imshow(plt.imread(file_img))
                plt.axis('off')
                fig.tight_layout()
                pdf.savefig(dpi=300)
    pdf.close()

def get_color_list(n):
    colors_default = ['Blue', 'Red', 'Green', 'm', 'Cyan', 'Orange', 'Yellow']
    if n <= len(colors_default):
        return colors_default[0:n]
    else:
        return colors_default[0]