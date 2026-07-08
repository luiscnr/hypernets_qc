
import os, shutil, argparse, configparser
import common_functions as cf
from datetime import timedelta
from datetime import datetime as dt
from datetime import timezone
from hypernets_day import HYPERNETS_DAY_LAND

parser = argparse.ArgumentParser(description="Creation of insitu nc files")
parser.add_argument('-m', "--mode",choices=['CREATEDAYFILES', 'REPORTDAYFILES'],required=True)
parser.add_argument('-sd', "--start_date", help="Start date. Optional with --listdates (YYYY-mm-dd)")
parser.add_argument('-ed', "--end_date", help="End date. Optional with --listdates (YYYY-mm-dd)")
parser.add_argument('-st', "--start_time", help="Start time. (HH:MM)")
parser.add_argument('-et', "--end_time", help="End time. (HH:MM)")
parser.add_argument('-i', "--input_path", help="Input path",required=True)
parser.add_argument('-o', "--output_path", help="Output path")
parser.add_argument('-c', "--config_path", help="Configuration file path")
parser.add_argument('-site', "--site_name", help="Site name")
parser.add_argument("-ow", "--overwrite", help="Overwrite output file(s).", action="store_true")
parser.add_argument('-ndel', "--nodelfiles", help="Do not delete temp files.", action="store_true")
parser.add_argument("-v", "--verbose", help="Verbose mode.", action="store_true")
#print(parser.parse_args())  # print the arguments

args = parser.parse_args()

#log should be WARNING or ERROR
def get_config_file(output_path,site,log='WARNING'):
    if args.config_path:
        config_file_summary = args.config_path
    else:
        config_file_summary = os.path.join(output_path, site, 'ConfigPlotSummary.ini')
        if not os.path.exists(config_file_summary):
            config_file_summary = os.path.join(output_path, 'ConfigPlotSummary.ini')

    if os.path.isfile(config_file_summary):
        print('[INFO] Configuration file: ', config_file_summary)
        return config_file_summary
    else:
        print(f'[{log}] Configuration file {config_file_summary} is not available')
        return None

def  make_report_files(input_path, output_path, site, start_date, end_date):
    if args.verbose:
        print(f'[INFO] Started creating reports')
    config_file_summary = get_config_file(output_path,site,log='ERROR')
    if config_file_summary is None:
        return
    options_c = cf.ConfigOptions(config_file_summary)
    work_date = start_date.replace(hour=0, minute=0, second=0, microsecond=0)
    hday = HYPERNETS_DAY_LAND(input_path, output_path)
    interval = 24

    daily_sequences_summary = None
    while work_date <= end_date:
        if args.verbose:
            print(f'--------------------------------------------------------------------------------------------------')
            print(f'[INFO] Date: {work_date}')
        sequence_abs_range = options_c.get_sequence_range(work_date,True)
        if sequence_abs_range is not None:
            print(
                f'[INFO] Absolute sequence range: {dt.fromtimestamp(sequence_abs_range[0]).astimezone(timezone.utc).strftime("%Y-%m-%d %H:%M")}-{dt.fromtimestamp(sequence_abs_range[1]).astimezone(timezone.utc).strftime("%Y-%m-%d %H:%M")}')

        hdayfile = hday.get_hypernets_day_file(site, work_date)
        if hday.get_input_folder_date(site,work_date) is None:
            sequences_all = []
            sequences_no_data = []
        else:
            if hdayfile is None:
                sequences_all = hday.get_sequences_date(site, work_date)
                sequences_all = [x[:-2] if len(x) == 18 else x for x in sequences_all]
                sequences_no_data = sequences_all
            else:
                sequences_no_data, sequences_all = hday.get_sequences_info(site, work_date, hdayfile.get_sequences(),sequence_abs_range)

        print(f'[INFO] Total number of sequences with folders: {len(sequences_all)}')
        print(f'[INFO] Sequences without L2 data: {len(sequences_no_data)}')

        output_folder_date = hday.get_output_folder_date(site,work_date)
        if output_folder_date is None:
            print(f'[ERROR] Path image date could not be created in {output_path}. Please review permissions')
            work_date = work_date + timedelta(hours=interval)
            continue

        if hdayfile is None:
            print(f'[WARNING] HYPERNETS day file for date {work_date} is not available. Creating empty daily summary plot and skipping...')
            file_img = os.path.join(output_folder_date, f'{site}_{work_date.strftime("%Y%m%d")}_DailySummary.png')
            cf.create_empty_image(file_img, site, work_date)
            delete = False if args.nodelfiles else True
            for seq in sequences_all:
                files_img = hday.get_files_img_for_sequences_no_data(site, work_date, seq, True)
                hday.save_report_image_only_pictures(site, delete, args.overwrite, seq, files_img, output_folder_date)
            sequence_range = options_c.get_sequence_range(work_date, False)
            daily_sequences_summary = {
                'NTotal': 0,
                'NAvailable': 0,
                'start_time': dt.fromtimestamp(sequence_range[0]).astimezone(timezone.utc).strftime('%Y-%m-%d %H:%M'),
                'end_time': dt.fromtimestamp(sequence_range[1]).astimezone(timezone.utc).strftime('%Y-%m-%d %H:%M'),
                'expected_sequences': sequence_range[2]
            }
            work_date = work_date + timedelta(hours=interval)
            continue

        hdayfile.set_path_images_date(site, work_date)

        file_summary = None

        dir_img_summary = os.path.join(os.path.dirname(hdayfile.file_nc), 'SUMMARY')
        file_summary = os.path.join(os.path.dirname(hdayfile.file_nc),f'{site}_{work_date.strftime("%Y%m%d")}_DailySummary{hdayfile.format_img}')
        if os.path.exists(file_summary) and not args.overwrite:
            print(f'[WARNING] Summary file: {output_path} already exist. Skipping...')
            if start_date == end_date:  ##required daily sequences summary:
                print(f'[INFO] Retrieving daily sequences summary...')
                daily_sequences_summary = cf.plot_from_options(hdayfile, config_file_summary, dir_img_summary,sequences_no_data, True)
                if os.path.isdir(dir_img_summary):
                    for name in os.listdir(dir_img_summary):
                        os.remove(os.path.join(dir_img_summary, name))
                    os.rmdir(dir_img_summary)
        else:
            daily_sequences_summary = cf.plot_from_options(hdayfile, config_file_summary, dir_img_summary,sequences_no_data, False,verbose=args.verbose)
            hdayfile.save_report_summary_image(site, work_date, dir_img_summary, daily_sequences_summary)

        delete = False if args.nodelfiles else True
        for seq in sequences_all:
            isequence = sequences_all[seq]
            if isequence >= 0:
                hdayfile.isequence = isequence
                hdayfile.save_report_image(site, delete, args.overwrite)
            else:
                hday.set_rgb_refs(config_file_summary)
                files_img = hday.get_files_img_for_sequences_no_data(site, work_date, seq, use_seq_folders)
                hdayfile.save_report_image_only_pictures(site, delete, args.overwrite, seq, files_img)

        cf.create_daily_pdf_report(input_path, output_path, site, work_date, file_summary, sequences_all,overwrite=args.overwrite)

        work_date = work_date + timedelta(hours=interval)

    if start_date == end_date:
        folder_day = hday.get_output_folder_date(site, start_date)
        date_str = start_date.strftime("%Y%m%d")
        name_summary = f'{site}_{date_str}_DailySummary.png'
        name_pdf = f'Report_{site}_{date_str}.pdf'
        file_pdf = os.path.join(folder_day, name_pdf)
        file_qc_mail = os.path.join(output_path, site, 'QCMail.mail')

        print(f'[INFO] Creating e-mail file: {file_qc_mail}')
        extra_info = {
            'folder_day': folder_day,
            'name_summary': name_summary,
            'file_pdf': file_pdf,
            'file_log_disk_usage': hday.get_disk_usage_log_file(site),
            'file_log_last_sequence': None
        }


        daily_mail = cf.DailyMail(file_qc_mail)
        daily_mail.create_file(site, start_date, daily_sequences_summary, extra_info)

def make_create_day_files(input_path, output_path, site, start_date, end_date):
    if args.verbose:
        print(f'[INFO] Started creating files')

    config_file_summary = get_config_file(output_path,site)

    work_date = start_date.replace(hour=0, minute=0, second=0, microsecond=0)
    interval = 24
    hday = HYPERNETS_DAY_LAND(input_path, output_path)
    if config_file_summary is not None:
        hday.set_rgb_refs(config_file_summary)

    while work_date <= end_date:
        if args.verbose:
            print(f'--------------------------------------------------------------------------------------------------')
            print(f'[INFO] Date: {work_date}')

        hday.get_files_date(site, work_date)

        if len(hday.files_dates) == 0:
            print(f'[WARNING] No data files found for the sequences on date: {work_date}. Skipping...')
            work_date = work_date + timedelta(hours=interval)
            continue

        if args.verbose:
            print(f'[INFO] Number of sequences: {len(hday.files_dates)}')

        nseq = hday.start_file_date(site, work_date, args.overwrite)
        if nseq <= 0:
            work_date = work_date + timedelta(hours=interval)
            continue
        else:
            if args.verbose:
                print(f'[INFO] Sequences with L2 data: : {nseq}')

        hday.set_data(site, work_date)
        hday.close_file_data()
        work_date = work_date + timedelta(hours=interval)

def main():
    print(f'[INFO] Started Daily Checking for land sites!')
    start_date, end_date = cf.get_start_and_end_dates(args)
    if start_date is None:
        return
    #start_time, end_time = cf.get_start_and_end_times(args)

    site = 'JSIT'
    if args.site_name:
        site = args.site_name
    input_path = cf.check_input_path(args)
    if input_path is None:
        return
    output_path = args.output_path if args.output_path else args.input_path
    output_path = cf.check_output_path(output_path)
    if output_path is None:
        return

    if args.verbose:
            print(f'[INFO] Input path set to: {input_path}')
            print(f'[INFO] Output path set to: {output_path}')


    if args.mode == 'CREATEDAYFILES':
        make_create_day_files(input_path, output_path,site, start_date, end_date)

    if args.mode == 'REPORTDAYFILES':
        make_report_files(input_path, output_path, site, start_date, end_date)

# %%
if __name__ == '__main__':
    main()