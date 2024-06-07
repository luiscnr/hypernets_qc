#HYPERNETS-Daily Checking

Usage:

1. Make sun plots

python make_hypernets_qc.py -m SUNPLOTS -site site –i input_path –o output_path -sd myDate -ed myDate -ndw –v

2. Prepare sun plot e-mail

python make_hypernets_qc.py -m SUNMAIL  -site site -i input_path -o output_path -sd myDate

3. Create daily NetCDF files

python make_hypernets_qc.py -m CREATEDAYFILES -site site –i input_path –o output_path -sd myDate -ed myDate -ndw -v

4. Create daily plots and PDF

python make_hypernets_qc.py -m REPORTDAYFILES -site site –i input_path –o output_path -sd myDate -ed myDate -ndw -v


Step 4 requires the configuration file: /output_path/site/ConfigPlotSummary.ini

Example of bash script for launching daily check for a single site (path needs to be added):


myDate=$( date --date yesterday +%F )

##VEIT-START
python $code_path/make_hypernets_qc.py -m SUNPLOTS -site VEIT –i $input_path –o $output_path -sd $myDate -ed $myDate -ndw –v
python $code_path/make_hypernets_qc.py -m CREATEDAYFILES -site VEIT –i $input_path –o $output_path -sd $myDate -ed $myDate -ndw -v
python $code_path/make_hypernets_qc.py -m REPORTDAYFILES -site VEIT –i $input_path –o $output_path -sd $myDate -ed $myDate -ndw -v
sunfile=$(python $code_path/make_hypernets_qc.py -m SUNMAIL  -site VEIT -i $input_path -o $output_path -sd $myDate)

mailrcpt=luis.gonzalezvilas@artov.ismar.cnr.it,vittorioernesto.brando@cnr.it,ivan.farace@artov.ismar.cnr.it

mailfileqc=$output_path/VEIT/QCMail.mail
mailfilesun=$output_path/VEIT/SunMail.mail
fileSummary=$(grep -i Summary $mailfileqc | cut -c 15-)
cat $mailfilesun >> $mailfileqc
subject="[ESA-POP] QUALITY CONTROL - VEIT - $myDate"
cat $mailfileqc | mail -a $fileSummary -a $sunfile -s "$subject" "$mailrcpt"
##VEIT-END

