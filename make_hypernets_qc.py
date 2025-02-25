import os, shutil, argparse, configparser
from datetime import datetime as dt
from datetime import timedelta
import numpy as np
import pandas as pd
#import owncloud
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from PIL import Image

from hypernets_day_file import HYPERNETS_DAY_FILE
from hypernets_day import HYPERNETS_DAY
from plot_options import PlotOptions
from flag_builder import FlagBuilder

print("starting to read arguments")
parser = argparse.ArgumentParser(description="Creation of insitu nc files")
parser.add_argument('-m', "--mode",
                    choices=['GETFILES', 'CREATEDAYFILES', 'REPORTDAYFILES', 'SUMMARYFILES', 'NCFROMCSV', 'PLOT',
                             'SUNDOWNLOAD', 'SUNPLOTS', 'SUNMAIL', 'CORRECTANGLES', 'COPYFROMCSV', 'SINGLESUN'],
                    required=True)
parser.add_argument('-sd', "--start_date", help="Start date. Optional with --listdates (YYYY-mm-dd)")
parser.add_argument('-ed', "--end_date", help="End date. Optional with --listdates (YYYY-mm-dd)")
parser.add_argument('-st', "--start_time", help="Start time. (HH:MM)")
parser.add_argument('-et', "--end_time", help="End time. (HH:MM)")
parser.add_argument('-i', "--input_path", help="Input path")
parser.add_argument('-o', "--output_path", help="Output path")
parser.add_argument('-c', "--config_path", help="Configuration file path")
parser.add_argument('-site', "--site_name", help="Site name")
parser.add_argument('-sopt', "--summary_options", help="Summary options,separated by '_': csv,nc,copy")
parser.add_argument('-ndays', "--ndays_interval", help="Interval days between start date and end date")
parser.add_argument('-ndel', "--nodelfiles", help="Do not delete temp files.", action="store_true")
parser.add_argument("-ndw", "--nodownload", help="No download (for launching without connection with RBINS).",
                    action="store_true")
parser.add_argument("-ow", "--overwrite", help="Overwrite output file(s).", action="store_true")
parser.add_argument("-v", "--verbose", help="Verbose mode.", action="store_true")

print(parser.parse_args())  # print the arguments

args = parser.parse_args()


def make_report_files(input_path, output_path, site, start_date, end_date):
    if args.verbose:
        print(f'[INFO] Started making report')
    work_date = start_date.replace(hour=0, minute=0, second=0, microsecond=0)
    hday = HYPERNETS_DAY(input_path, output_path)
    interval = 24
    if args.ndays_interval:
        interval = 24 * int(args.ndays_interval)

    if args.config_path:
        config_file_summary = args.config_path
    else:
        config_file_summary = os.path.join(output_path, site, 'ConfigPlotSummary.ini')
        if not os.path.exists(config_file_summary):
            config_file_summary = os.path.join(output_path, 'ConfigPlotSummary.ini')
    print('[INFO] Config file summary: ', config_file_summary)

    daily_sequences_summary = None
    while work_date <= end_date:
        if args.verbose:
            print(f'--------------------------------------------------------------------------------------------------')
            print(f'[INFO] Date: {work_date}')

        sequence_abs_range = None
        if os.path.exists(config_file_summary):
            sequence_abs_range = hday.get_sequence_range(work_date, config_file_summary, True)
        if sequence_abs_range is not None:
            print(
                f'[INFO] Absolute sequence range: {dt.utcfromtimestamp(sequence_abs_range[0]).strftime("%Y-%m-%d %H:%M")}-{dt.utcfromtimestamp(sequence_abs_range[1]).strftime("%Y-%m-%d %H:%M")}')

        hdayfile = hday.get_hypernets_day_file(site, work_date)
        use_seq_folders = hday.check_use_seq_folders(site, work_date)
        if use_seq_folders is None:##no data folder
            sequences_all = []
            sequences_no_data = sequences_all
        else:
            if hdayfile is None:
                sequences_all = hday.get_sequences_date_from_file_list(site, work_date)
                sequences_all = [x[:-2] if len(x) == 18 else x for x in sequences_all]
                sequences_no_data = sequences_all
            else:
                sequences_no_data, sequences_all = hday.get_sequences_info(site, work_date, hdayfile.get_sequences(),
                                                                       sequence_abs_range)

        print(f'[INFO] Total number of sequences with folders: {len(sequences_all)}')
        print(f'[INFO] Sequences without L2 data: {len(sequences_no_data)}')

        output_folder_date = hday.get_output_folder_date(site, work_date)
        if output_folder_date is None:
            print(f'[ERROR] Path image date could not be created in {output_path}. Please review permissions')
            work_date = work_date + timedelta(hours=interval)
            continue
        if hdayfile is None:
            print(
                f'[WARNING] HYPERNETS day file for date {work_date} is not available. Creating empty daily summary plot and skipping...')
            file_img = os.path.join(output_folder_date, f'{site}_{work_date.strftime("%Y%m%d")}_DailySummary.png')
            create_empty_image(file_img, site, work_date)
            delete = False if args.nodelfiles else True
            for seq in sequences_all:
                if os.path.exists(config_file_summary):
                    hday.set_rgb_refs(config_file_summary)
                files_img = hday.get_files_img_for_sequences_no_data(site, work_date, seq, use_seq_folders)
                hday.save_report_image_only_pictures(site, delete, args.overwrite, seq, files_img, output_folder_date)
            sequence_range = hday.get_sequence_range(work_date, config_file_summary, False)
            daily_sequences_summary = {
                'NTotal': 0,
                'NAvailable': 0,
                'start_time': dt.utcfromtimestamp(sequence_range[0]).strftime('%Y-%m-%d %H:%M'),
                'end_time': dt.utcfromtimestamp(sequence_range[1]).strftime('%Y-%m-%d %H:%M'),
                'expected_sequences': sequence_range[2]
            }
            work_date = work_date + timedelta(hours=interval)
            continue
        hdayfile.set_path_images_date(site, work_date)

        file_summary = None
        if os.path.exists(config_file_summary):
            dir_img_summary = os.path.join(os.path.dirname(hdayfile.file_nc), 'SUMMARY')
            file_summary = os.path.join(os.path.dirname(hdayfile.file_nc),
                                        f'{site}_{work_date.strftime("%Y%m%d")}_DailySummary{hdayfile.format_img}')
            if os.path.exists(file_summary) and not args.overwrite:
                print(f'[WARNING] Summary file: {output_path} alreaday exist. Skipping...')
                if start_date == end_date:  ##required daily sequences summary:
                    print(f'[INFO] Retrieving daily sequences summary...')
                    daily_sequences_summary = plot_from_options(hdayfile, config_file_summary, dir_img_summary,
                                                                sequences_no_data, True)

                    if os.path.isdir(dir_img_summary):
                        for name in os.listdir(dir_img_summary):
                            os.remove(os.path.join(dir_img_summary,name))
                        os.rmdir(dir_img_summary)
            else:
                daily_sequences_summary = plot_from_options(hdayfile, config_file_summary, dir_img_summary,
                                                            sequences_no_data, False)
                hdayfile.save_report_summary_image(site, work_date, dir_img_summary, daily_sequences_summary, input_path)

        delete = False if args.nodelfiles else True
        for seq in sequences_all:
            isequence = sequences_all[seq]
            if isequence >= 0:
                hdayfile.isequence = isequence
                hdayfile.save_report_image(site, delete, args.overwrite)
            else:
                hday.set_rgb_refs(config_file_summary)
                files_img = hday.get_files_img_for_sequences_no_data(site, work_date, seq,use_seq_folders)
                hdayfile.save_report_image_only_pictures(site, delete, args.overwrite, seq, files_img)

        create_daily_pdf_report(input_path, output_path, site, work_date, file_summary, sequences_all)

        work_date = work_date + timedelta(hours=interval)

    if start_date == end_date:
        hday = HYPERNETS_DAY(input_path, output_path)
        folder_day = hday.get_output_folder_date(site, start_date)
        date_str = start_date.strftime("%Y%m%d")
        name_summary = f'{site}_{date_str}_DailySummary.png'
        name_pdf = f'Report_{site}_{date_str}.pdf'
        file_pdf = os.path.join(folder_day, name_pdf)
        file_qc_mail = os.path.join(output_path, site, 'QCMail.mail')
        public_link = ''
        owncloud_info = None
        if os.path.exists(config_file_summary):
            options = configparser.ConfigParser()
            options.read(config_file_summary)
            if options.has_option('GLOBAL_OPTIONS', f'public_link'):
                public_link = options['GLOBAL_OPTIONS'][f'public_link'].strip()
            owncloud_options = ['owncloud_client', 'owncloud_user', 'owncloud_password']
            owncloud_info = {}
            for owc in owncloud_options:
                if options.has_option('GLOBAL_OPTIONS', owc):
                    owncloud_info[owc] = options['GLOBAL_OPTIONS'][owc]
                else:
                    owncloud_info = None
                    break

        # public_link = 'https://file.sic.rm.cnr.it/index.php/s/rBeO2UMtdJ4F3Gx'
        print(f'[INFO] Creating e-mail file: {file_qc_mail}')
        extra_info = {
            'folder_day': folder_day,
            'name_summary': name_summary,
            'file_pdf': file_pdf,
            'public_link': public_link,
            'file_log_disk_usage': hday.get_disk_usage_log_file(site, True),
            'file_log_last_sequence': None#hday.get_last_available_log(site, 'sequence', True)
        }
        print("extra_info: ", extra_info)
        
        create_daily_mail_file(file_qc_mail, site, start_date, daily_sequences_summary, extra_info)
        print(f'[INFO] PDF file: {file_pdf} --> {os.path.exists(file_pdf)}')
        print(f'[INFO] Owncloud info: {owncloud_info}')
        if os.path.exists(file_pdf) and owncloud_info is not None:
            import owncloud
            session = owncloud.Client(owncloud_info['owncloud_client'])
            session.login(owncloud_info['owncloud_user'], owncloud_info['owncloud_password'])
            session.put_file(f'/ESA-HYP-POP/LastQC_Reports/{site}_LastQC.pdf', file_pdf)


def create_daily_mail_file(file_qc_mail, site, start_date, daily_sequences_summary, extra_info):
    fout = open(file_qc_mail, 'w')
    fout.write(f'QUALITY CONTROL - {site} - {start_date.strftime("%Y-%m-%d")}')
    add_new_line(fout, '===================================')
    add_new_line(fout, '')
    if daily_sequences_summary is not None:
        if daily_sequences_summary['NTotal'] == 0:
            add_new_line(fout, f'WARNING: No folder sequences found for {site} on {start_date.strftime("%Y-%m-%d")}')
            add_new_line(fout, f'Please review if the system is working and folders are available in the server')
            add_new_line(fout, '')

        add_new_line(fout, 'SEQUENCES SUMMARY')
        add_new_line(fout, '=================')
        add_new_line(fout, f'Start time:  {daily_sequences_summary["start_time"]}')
        add_new_line(fout, f'End time: {daily_sequences_summary["end_time"]}')
        add_new_line(fout, f'Expected sequences: {daily_sequences_summary["expected_sequences"]}')
        add_new_line(fout, f'Available sequences: {daily_sequences_summary["NTotal"]}')
        add_new_line(fout, f'Sequences processed to L2: {daily_sequences_summary["NAvailable"]}')
        if 'VALID' in daily_sequences_summary.keys():
            add_new_line(fout, f'Valid sequences after quality control: {daily_sequences_summary["VALID"]}')
        add_new_line(fout, '')




    add_new_line(fout, 'SYSTEM STATUS')
    add_new_line(fout, '=============')

    ##disk usage
    file_log_disk_usage = extra_info['file_log_disk_usage']
    if file_log_disk_usage is not None:
        lines_disk_usage = get_lines_disk_usage(file_log_disk_usage)
        for line in lines_disk_usage:
             add_new_line(fout, line)

    ##last log sequence
    file_log_last_sequence = extra_info['file_log_last_sequence']  ##None, not implemented
    if file_log_last_sequence is not None:
        ##no implemented
        pass

    add_new_line(fout,'')
    add_new_line(fout, 'DAILY CHECKING FILES')
    add_new_line(fout, '====================')
    add_new_line(fout, f'Ouput folder: {extra_info["folder_day"]}')
    file_summary = os.path.join(extra_info["folder_day"], extra_info["name_summary"])
    add_new_line(fout, f'Summary file: {file_summary if os.path.exists(file_summary) else "Not. Av."}')
    add_new_line(fout, '')
    add_new_line(fout, f'PDF file: {extra_info["file_pdf"] if os.path.exists(extra_info["file_pdf"]) else "Not. Av."}')
    if os.path.exists(extra_info["file_pdf"]):
        add_new_line(fout, f'Link to PDF file: {extra_info["public_link"]}')
    add_new_line(fout, '')

    fout.close()


def get_lines_disk_usage(file_log):
    lines = ['']
    if not os.path.exists(file_log): return lines
    df = pd.read_csv(file_log, sep=' ')
    lines.append('DISK USAGE')
    lines.append('----------')
    last_line = df.iloc[-1]
    used = float(last_line.iloc[1]) / (1024 * 1024)
    av = float(last_line.iloc[2]) / (1024 * 1024)
    lines.append(f' Last measurement: {last_line.iloc[0]} Used: {used:.2f} Gb. Available: {av:.2f} Gb. %Use: {last_line.iloc[4]}')

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
        if pd.isna(line_here.iloc[0]):continue
        date_here_str = str(line_here.iloc[0])[:10]


        if date_here_str != date_ref_str_loop:
            date_ref_str_loop = date_here_str
            date_here_str_basic = dt.strptime(date_here_str, '%Y-%m-%d').replace(hour=12, minute=0, second=12).strftime(
                '%Y-%m-%d')
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
            lines.append(f'  {line_here.iloc[0]} Used: {used:.2f} Gb. Available: {av:.2f} Gb. %Use: {line_here.iloc[4]}')
    lines.append('')

    return lines


def add_new_line(fout, str):
    fout.write('\n')
    fout.write(str)


def create_empty_image(file_img, site, date_here):
    plt.figure(figsize=(6, 0.75))
    plt.title(f'L2 data were not available for {site} on {date_here.strftime("%Y-%m-%d")}')
    plt.xticks([])
    plt.yticks([])
    plt.axis('off')
    plt.tight_layout()
    plt.savefig(file_img, dpi=300)


def create_daily_pdf_report(input_path, output_path, site, date_here, file_summary, sequences):
    hday = HYPERNETS_DAY(input_path, output_path)
    folder_day = hday.get_output_folder_date(site, date_here)
    date_here_str = date_here.strftime('%Y%m%d')
    file_pdf = os.path.join(folder_day, f'Report_{site}_{date_here_str}.pdf')
    print(f'[INFO] PDF report: {file_pdf}')
    if os.path.exists(file_pdf) and not args.overwrite:
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


def make_get_files(input_path, site, start_date, end_date):
    if args.verbose:
        print(f'[INFO] Started getting files')
    work_date = start_date.replace(hour=0, minute=0, second=0, microsecond=0)

    hday = HYPERNETS_DAY(input_path, None)
    while work_date <= end_date:
        if args.verbose:
            print(f'--------------------------------------------------------------------------------------------------')
            print(f'[INFO] Date: {work_date}')
        hday.get_files_date(site, work_date)
        work_date = work_date + timedelta(hours=24)


def make_create_dayfiles(input_path, output_path, site, start_date, end_date):
    if args.verbose:
        print(f'[INFO] Started creating files')

    if args.config_path:
        config_file_summary = args.config_path
    else:
        config_file_summary = os.path.join(output_path, site, 'ConfigPlotSummary.ini')
        if not os.path.exists(config_file_summary):
            config_file_summary = os.path.join(output_path, 'ConfigPlotSummary.ini')
    print('[INFO] Config file summary: ', config_file_summary)

    work_date = start_date.replace(hour=0, minute=0, second=0, microsecond=0)
    interval = 24
    if args.ndays_interval:
        interval = 24 * int(args.ndays_interval)
    hday = HYPERNETS_DAY(input_path, output_path)
    if os.path.exists(config_file_summary):
        hday.set_rgb_refs(config_file_summary)

    while work_date <= end_date:
        if args.verbose:
            print(f'--------------------------------------------------------------------------------------------------')
            print(f'[INFO] Date: {work_date}')
        if args.nodownload:
            hday.get_files_date_nodownload(site, work_date)
        else:
            hday.get_files_date(site, work_date)

        if len(hday.files_dates) == 0:
            print(f'[WARNING] No data files found for the sequences on date: {work_date}. Skipping...')
            work_date = work_date + timedelta(hours=interval)
            continue

        if args.verbose:
            print(f'[INFO] Number of sequences: {len(hday.files_dates)}')
        if site == 'JSIT':  ##LAND
            nseq = hday.start_file_date_land_complete(site, work_date, args.overwrite)
        else:
            nseq = hday.start_file_date_complete(site, work_date, args.overwrite)
        if nseq < 0:
            work_date = work_date + timedelta(hours=interval)
            continue
        else:
            if args.verbose:
                print(f'[INFO] Sequences with L2 data: : {nseq}')
        if site == 'JSIT':  ##LAND
            hday.set_data_land(site, work_date)
        else:
            hday.set_data(site, work_date)
        hday.close_datafile_complete()
        work_date = work_date + timedelta(hours=interval)


def plot_from_options(hfile, config_file, output_path_images, sequences_no_data, only_sequences_summary):
    if not os.path.exists(config_file):
        print(f'[ERROR] Plot configuration file: {config_file} does not exist. ')
        return None

    if args.verbose:
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


def make_flagged_nc_from_csv(config_file):
    if args.verbose:
        print(f'[INFO] Started concatenaded nc file from csv using configuration file: {config_file}')

    options = configparser.ConfigParser()
    options.read(config_file)
    if not options.has_section('basic'):
        print(f'[ERROR] Basic section should be included in the configuration file')
        return
    input_csv = None
    if options.has_option('basic', 'input_csv'):
        input_csv = options['basic']['input_csv'].strip()
    if input_csv is None:
        print(f'[ERROR] input_csv should be included in the basic section of the configuration file')
        return
    if not os.path.exists(input_csv):
        print(f'[ERROR] File {input_csv} does not exist')
        return

    if options.has_option('basic', 'output_path'):
        output_path = options['basic']['output_path'].strip()
    else:
        output_path = os.path.dirname(input_csv)
    if not os.path.isdir(output_path):
        try:
            os.mkdir(output_path)
        except:
            print(f'[ERROR] Output path {output_path} does not exist')
            return

    col_sequence = 'sequence_ref'
    col_file_nc = 'file_nc'
    col_isequence = 'isequence'
    flags = get_flags(options)
    name_csv = input_csv.split('/')[-1]
    name_csv = name_csv[:-4]
    file_nc = os.path.join(output_path, f'CSV_NC_{name_csv}.nc')

    df = pd.read_csv(input_csv, sep=';')

    file_nc_prev = ''
    dataset_w = None
    hdayfile = None
    sindices = []
    pindices = []
    index_w = 0
    for index, row in df.iterrows():
        file_nc_here = str(row[col_file_nc])
        if not os.path.exists(file_nc_here):
            continue
        isequence = int(row[col_isequence])

        if file_nc_here != file_nc_prev:
            if dataset_w is not None:  ##set data
                dataset_w = hdayfile.set_data_dataset_w(dataset_w, sindices, index_w)
                for flag in flags:
                    print(pindices)
                    print(flags)
                    print(flag)
                    flag_array = get_flag_array(pindices, df, flags, flag)
                    print(flag_array)
                    dataset_w = hdayfile.set_data_flag(dataset_w, index_w, flag, flag_array)
                index_w = index_w + len(sindices)
                hdayfile = HYPERNETS_DAY_FILE(file_nc_here, None)
                sindices = [isequence]
                pindices = [int(index)]
            else:
                hdayfile = HYPERNETS_DAY_FILE(file_nc_here, None)
                dataset_w = hdayfile.start_dataset_w(file_nc)
                dataset_w = hdayfile.add_flag_variables(dataset_w, flags)
                sindices = [isequence]
                pindices = [int(index)]
            file_nc_prev = file_nc_here
        else:
            sindices.append(isequence)
            pindices.append(int(index))

    if dataset_w is not None:
        dataset_w.close()


def make_copy_from_csv(site, input_path, output_directory):
    import pandas as pd
    from datetime import datetime as dt
    base_folder = '/store3/HYPERNETS/INPUT_HYPSTARv2.0_QC/'
    df = pd.read_csv(input_path, sep=';')
    for index, row in df.iterrows():
        sref = row['sequence_ref']
        date_ref = dt.strptime(sref, '%Y%m%dT%H%M')
        date_folder = os.path.join(base_folder, site, date_ref.strftime('%Y'), date_ref.strftime('%m'),
                                   date_ref.strftime('%d'))
        report_name = f'{site}_{sref}_Report.png'
        input_file = os.path.join(date_folder, report_name)
        output_file = os.path.join(output_directory, report_name)
        if os.path.exists(input_file):
            shutil.copy(input_file, output_file)
            print(f'[INFO] Input file: {input_file} coppied to {output_directory}')
        else:
            print(f'[WARNING] Input file: {input_file} is not available. Skipping...')


def get_flags(options):
    flags = {}
    if options.has_section('flags'):
        index = 0
        while options.has_option('flags', f'flag_{index}.name'):
            name = options['flags'][f'flag_{index}.name'].strip()
            flag_values = []
            flag_meanings = []
            if options.has_option('flags', f'flag_{index}.values'):
                flag_values_str = options['flags'][f'flag_{index}.values'].strip()
                flag_values = [int(x) for x in flag_values_str.split(',')]
            if options.has_option('flags', f'flag_{index}.values'):
                flag_meanings_str = options['flags'][f'flag_{index}.meanings'].strip()
                flag_meanings = [str(x) for x in flag_meanings_str.split(',')]
            flags[name] = {
                'values': flag_values,
                'meanings': flag_meanings
            }
            index = index + 1
    return flags


def get_flag_array(indices, df, flags, flag):
    print('Columna: ', df.columns)
    print('Indices: ', indices)
    print(type(indices))
    print('Flag', flag)
    data_flag = df.loc[indices, flag]
    try:
        data_flag_array = np.array(data_flag).astype(np.int64)
    except:
        for idx in range(len(flags[flag]['values'])):
            val = flags[flag]['values'][idx]
            meaning = flags[flag]['meanings'][idx].strip()
            data_flag[data_flag == meaning] = val
        data_flag_array = np.array(data_flag).astype(np.int64)

    return data_flag_array


def make_summary_files(input_path, output_path, site, start_date, end_date, start_time, end_time, options):
    if args.verbose:
        print(f'[INFO] Started summary files with options: {options}')

    import pandas as pd
    work_date = start_date.replace(hour=0, minute=0, second=0, microsecond=0)

    hday = HYPERNETS_DAY(None, input_path)
    interval = 24
    if args.ndays_interval:
        interval = 24 * int(args.ndays_interval)

    if not os.path.isdir(output_path):
        try:
            os.mkdir(output_path)
        except:
            print(f'[ERROR] Output path {output_path} is not a valid directory and could not be created')
            return

    start_date_str = start_date.strftime('%Y%m%d')
    end_date_str = end_date.strftime('%Y%m%d')
    dataset_w = None
    index_w = 0
    df_csv = None
    col_names_csv = None

    while work_date <= end_date:
        if args.verbose:
            print(f'--------------------------------------------------------------------------------------------------')
            print(f'[INFO] Date: {work_date}')
        hdayfile = hday.get_hypernets_day_file(site, work_date)
        if hdayfile is None:
            print(f'[WARNING] HYPERNETS day file for date {work_date} is not available. Skipping...')
            work_date = work_date + timedelta(hours=interval)
            continue
        work_time_str = work_date.strftime('%Y-%m-%d')
        start_work_time = dt.strptime(f'{work_time_str}T{start_time}', '%Y-%m-%dT%H:%M')
        end_work_time = dt.strptime(f'{work_time_str}T{end_time}', '%Y-%m-%dT%H:%M')
        sequences, sequence_indices = hdayfile.get_sequences_interval(start_work_time, end_work_time)
        if 'copy' in options:
            report_files = hdayfile.get_report_files_interval(sequences, site, start_time, end_time)
            if args.verbose:
                print(f'[INFO] -> Copy {len(report_files)} report files')
            if len(report_files) > 0:
                for rinfile in report_files:
                    routfile = os.path.join(output_path, os.path.basename(rinfile))
                    shutil.copy(rinfile, routfile)
        if 'nc' in options and sequence_indices is not None:
            if dataset_w is None:
                file_nc = os.path.join(output_path, f'Comparison_{site}_{start_date_str}_{end_date_str}.nc')
                dataset_w = hdayfile.start_dataset_w(file_nc)
            print('Setting data->', len(sequence_indices))
            dataset_w = hdayfile.set_data_dataset_w(dataset_w, sequence_indices, index_w)
            index_w = index_w + len(sequence_indices)

        if 'csv' in options and sequence_indices is not None:
            if col_names_csv is None:
                col_names_csv = hdayfile.get_csv_col_names()
                df_csv = hdayfile.get_dataframe_lines(sequence_indices, col_names_csv)
            else:
                df_here = hdayfile.get_dataframe_lines(sequence_indices, col_names_csv)
                df_csv = pd.concat([df_csv, df_here])

        work_date = work_date + timedelta(hours=interval)

    if df_csv is not None:
        file_csv = os.path.join(output_path, f'Comparison_{site}_{start_date_str}_{end_date_str}.csv')
        df_csv.to_csv(file_csv, sep=';')
    if dataset_w is not None:
        dataset_w.close()


def make_sun_download(input_path, site, start_date, end_date):
    if args.verbose:
        print(f'[INFO] Downloading sun images')
    interval = 24
    if args.ndays_interval:
        interval = 24 * int(args.ndays_interval)
    work_date = start_date.replace(hour=0, minute=0, second=0, microsecond=0)
    hday = HYPERNETS_DAY(input_path, input_path)

    while work_date <= end_date:
        if args.verbose:
            print(f'--------------------------------------------------------------------------------------------------')
            print(f'[INFO] Date: {work_date}')
            hday.get_sun_images_date(site, work_date, False)

        work_date = work_date + timedelta(hours=interval)


def make_sun_plots(input_path, output_path, site, start_date, end_date, ndw):
    if args.verbose:
        print(f'[INFO] Started creating files')
        print(f'[INFO] No download set to: {ndw}')
    interval = 24
    if args.ndays_interval:
        interval = 24 * int(args.ndays_interval)
    work_date = start_date.replace(hour=0, minute=0, second=0, microsecond=0)

    if args.config_path:
        config_file_summary = args.config_path
    else:
        config_file_summary = os.path.join(output_path, site, 'ConfigPlotSummary.ini')
        if not os.path.exists(config_file_summary):
            config_file_summary = os.path.join(output_path, 'ConfigPlotSummary.ini')
    print('[INFO] Config file summary: ', config_file_summary)

    hday = HYPERNETS_DAY(input_path, output_path)
    ndays = 5
    hours = ['11:00', '12:00', '13:00', '14:00', '15:00']
    max_time_diff = 30 * 60  ##30 minutes
    if os.path.exists(config_file_summary):
        import configparser
        options = configparser.ConfigParser()
        options.read(config_file_summary)
        if options.has_option('GLOBAL_OPTIONS', 'sun_times'):
            value = options['GLOBAL_OPTIONS']['sun_times'].strip()
            value_s = value.split(',')
            sun_hours = []
            for x in value_s:
                try:
                    dt.strptime(x.strip(), '%H:%M')
                    sun_hours.append(x.strip())
                except:
                    print(f'[WARNING] {x.strip()} is not a valid sun time in format HH:MM. Using default values')
            if len(sun_hours) == 6:
                print(f'[INFO] Sun time hours set to: {sun_hours}')
                hours = sun_hours

        if options.has_option('GLOBAL_OPTIONS', 'sun_max_time_diff'):
            value = options['GLOBAL_OPTIONS']['sun_max_time_diff'].strip()
            try:
                max_time_diff = float(value) * 60
                print(f'[INFO] Maximum time difference to select sun images set to: {max_time_diff} minutes')
            except:
                print(
                    f'[WARNING] {value} is not a valid max. time difference for sun images. Use a numeric value (minutes)')

    nhours = len(hours)
    sun_images_list = [[''] * ndays] * nhours
    sun_images_time_list = [[''] * ndays] * nhours

    while work_date <= end_date:
        if args.verbose:
            print(f'--------------------------------------------------------------------------------------------------')
            print(f'[INFO] Date: {work_date}')

        for iday in range(ndays):
            work_date_here = work_date - timedelta(days=iday)
            print(f'[INFO] --> {work_date_here}')
            times_tmp = []
            for ihour in range(nhours):
                times_tmp.append(f'{work_date_here.strftime("%Y%m%d")}T{hours[ihour].replace(":", "")}')
            sun_images_time_list[iday] = times_tmp
            sun_images = hday.get_sun_images_date(site, work_date_here, ndw)
            sun_images_list_day = [None] * nhours
            sun_images_time_diffs = [max_time_diff] * nhours
            sun_images_time_refs = [dt.strptime(f'{work_date_here.strftime("%Y-%m-%d")} {x}', '%Y-%m-%d %H:%M') for x in
                                    hours]
            for sun_image_time_str in sun_images:
                sun_image_time = dt.strptime(sun_image_time_str, '%Y%m%dT%H%M%S')
                for idx in range(nhours):
                    diff_here = abs((sun_image_time - sun_images_time_refs[idx]).total_seconds())
                    if diff_here <= sun_images_time_diffs[idx]:
                        sun_images_time_diffs[idx] = diff_here
                        sun_images_list_day[idx] = sun_images[sun_image_time_str]

            sun_images_list[iday] = sun_images_list_day

        dir_out = hday.get_output_folder_date(site, work_date)

        file_out = os.path.join(dir_out, f'{site}_SunImages_{work_date.strftime("%Y%m%d")}.png')
        hday.save_sun_images(file_out, sun_images_list, sun_images_time_list)
        work_date = work_date + timedelta(hours=interval)


def get_start_and_end_dates():
    start_date = dt.now()
    if args.start_date:
        try:
            start_date = dt.strptime(args.start_date, '%Y-%m-%d')
        except:
            print(f'[ERROR] Start date is not in valid format YYYY-mm-DD')
            return None, None
    if args.end_date:
        try:
            end_date = dt.strptime(args.end_date, '%Y-%m-%d')
        except:
            print(f'[ERROR] End date is not in valid format YYYY-mm-DD')
            return None, None
    else:
        end_date = start_date

    return start_date, end_date


def get_start_and_end_times():
    start_time = '00:00'
    end_time = '23:59'
    if args.start_time:
        try:
            start_time_p = dt.strptime(args.start_time, '%H:%M')
            start_time = start_time_p.strftime('%H:%M')
        except:
            pass
    if args.end_time:
        try:
            end_time_p = dt.strptime(args.end_time, '%H:%M')
            end_time = end_time_p.strftime('%H:%M')
        except:
            pass

    return start_time, end_time


def prepare_sun_plot_email(input_path, output_path, site, start_date):
    output_path_site = os.path.join(output_path, site)
    file_mail_text = os.path.join(output_path_site, 'SunMail.mail')
    lines = [f'SUN PICTURES CHECK {start_date.strftime("%Y-%m-%d")}']
    ini_date = start_date - timedelta(days=5)
    lines.append(f'Start date for checking: {ini_date.strftime("%Y-%m-%d")}')
    lines.append(f'End date for checking: {start_date.strftime("%Y-%m-%d")}')
    path_sun = os.path.join(output_path, site, start_date.strftime('%Y'), start_date.strftime('%m'),
                            start_date.strftime('%d'))
    path_file = os.path.join(path_sun, f'{site}_SunImages_{start_date.strftime("%Y%m%d")}.png')
    if os.path.exists(path_file):
        lines.append(f'Image file path: {path_file}')
    else:
        lines.append(f'Error: No image file was found for date: {start_date.strftime("%Y-%m-%d")}')
    fout = open(file_mail_text, 'w')
    fout.write(lines[0])
    for idx in range(1, len(lines)):
        fout.write('\n')
        fout.write(lines[idx])
    fout.close()
    if os.path.exists(path_file):
        print(path_file)
    else:
        print('NONE')


def correct_angles(input_path, output_path, site, start_date, end_date):
    print('[INFO] Correctiong angles')
    interval = 24
    if args.ndays_interval:
        interval = 24 * int(args.ndays_interval)
    work_date = start_date.replace(hour=0, minute=0, second=0, microsecond=0)
    hd = HYPERNETS_DAY(input_path, output_path)
    while work_date <= end_date:
        if args.verbose:
            print(f'[INFO] Date: {work_date.strftime("%Y-%m-%d")}')
        hday = hd.get_hypernets_day_file(site, work_date)
        input_path_date = hd.get_folder_date(site, work_date)
        if hday is not None:
            list_metadata = hday.get_metadata_files(input_path_date)
            array_angles = hday.get_angles_from_metadata_files(list_metadata)
            hday.creating_copy_with_new_angles(array_angles)
        else:
            print(f'[WARNING] QC file for date {work_date.strftime("%Y-%m-%d")} is not available. Skipping...')
        work_date = work_date + timedelta(hours=interval)


def make_single_sun(site, sequence, output_path):
    if not os.path.exists(output_path):
        try:
            os.mkdir(output_path)
        except:
            print(f'[ERROR] {output_path} does not exist and could not be created')
            return
    import subprocess
    base_download = f'rsync -a -e \'ssh -p 9022\' hypstar@enhydra.naturalsciences.be:/home/hypstar/'
    suffix = 'RADIOMETER/01_016_0000_0_0000.jpg'
    cmd = f'{base_download}{site}/DATA/{sequence}/{suffix} {output_path}'
    prog = subprocess.Popen(cmd, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
    out, err = prog.communicate()
    if err:
        print(err)

    file_sun = os.path.join(output_path, '01_016_0000_0_0000.jpg')
    if not os.path.exists(file_sun):
        print(f'[ERROR] {file_sun} could not be downloaded')
        return
    file_out = os.path.join(output_path, f'{site}_SunImage_{sequence}.jpg')

    image = Image.open(file_sun)
    rimage = image.rotate(270, expand=True)
    fig = plt.figure()
    axhere = fig.gca()
    axhere.imshow(rimage)

    w, h = rimage.size

    ##central point
    axhere.axvline(w / 2, 0.48, 0.52, color='red', linewidth=0.25)
    axhere.axhline(h / 2, 0.48, 0.52, color='red', linewidth=0.25)
    ##grid
    incremx = int(w / 4)
    incremy = int(h / 4)
    for x in range(0, w, incremx):
        axhere.axvline(x, color='red', linewidth=0.5)
    for y in range(0, h, incremy):
        axhere.axhline(y, color='red', linewidth=0.5)

    axhere.set_xticks([])
    axhere.set_yticks([])
    axhere.set_title(sequence)
    plt.tight_layout()
    plt.savefig(file_out, bbox_inches='tight')
    os.remove(file_sun)
    print(f'Completed. File saved: {file_out}')


def do_test():
    ##for launching test
    return True


def main():
    # if do_test():
    #     return

    if args.verbose:
        print('STARTED')

    start_date, end_date = get_start_and_end_dates()
    if start_date is None:
        return
    start_time, end_time = get_start_and_end_times()

    site = 'VEIT'
    if args.site_name:
        site = args.site_name
    input_path = None
    output_path = None

    if args.input_path:
        input_path = args.input_path
        if args.verbose:
            print(f'[INFO] Input path set to: {input_path}')
    if args.output_path:
        output_path = args.output_path
        if args.verbose:
            print(f'[INFO] Output path set to: {output_path}')

    summary_options = ['copy']
    if args.summary_options:
        summary_options = args.summary_options.split('_')

    if args.mode == 'GETFILES':
        make_get_files(input_path, site, start_date, end_date)

    if args.mode == 'CREATEDAYFILES':
        make_create_dayfiles(input_path, output_path, site, start_date, end_date)

    if args.mode == 'REPORTDAYFILES':
        if args.config_path and not os.path.exists(args.config_path):
            print(f'[ERROR] Configuration file: {args.config_path} does not exist')
            return
        make_report_files(input_path, output_path, site, start_date, end_date)

    if args.mode == 'SUMMARYFILES':
        make_summary_files(input_path, output_path, site, start_date, end_date, start_time, end_time, summary_options)

    if args.mode == 'NCFROMCSV':
        make_flagged_nc_from_csv(args.config_path)

    if args.mode == 'COPYFROMCSV':
        make_copy_from_csv(site, input_path, output_path)

    if args.mode == 'PLOT':
        plot_from_options(input_path, args.config_path, None)

    if args.mode == 'SUNDOWNLOAD':
        make_sun_download(input_path, site, start_date, end_date)

    if args.mode == 'SUNPLOTS':
        make_sun_plots(input_path, output_path, site, start_date, end_date, args.nodownload)

    if args.mode == 'SUNMAIL':
        prepare_sun_plot_email(input_path, output_path, site, start_date)

    if args.mode == 'CORRECTANGLES':
        correct_angles(input_path, output_path, site, start_date, end_date)

    if args.mode == 'SINGLESUN':
        sequence = args.input_path
        make_single_sun(site, sequence, output_path)


# %%
if __name__ == '__main__':
    main()
