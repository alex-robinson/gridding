#!/bin/bash

search='YEARNOW'
templatename='submit_slurm'
filename='submit_slurm_now'

for year in {1979..2019}
do
    
#  # Insert the current year into the slurm submit script as an argument  
#  sed "s/$search/$year/g" $templatename > $filename
#
#  # Submit the job 
#  sbatch $filename
   #python get-era5-hourly-local.py $year
   python get-era5-monthly-single-levels.py $year
done
