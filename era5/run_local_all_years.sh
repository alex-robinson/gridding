#!/bin/bash

search='YEARNOW'
templatename='submit_slurm'
filename='submit_slurm_now'
var='2m_temperature'

# for year in {1979..2020}
# do
    
# #  # Insert the current year into the slurm submit script as an argument  
# #  sed "s/$search/$year/g" $templatename > $filename
# #
# #  # Submit the job 
# #  sbatch $filename
#    #python get-era5-hourly-local.py $year
#    ./get-era5-monthly-single-levels.py $var $year
# done

for v in {0:23}
do
for year in {1979..2020}
do
    
#  # Insert the current year into the slurm submit script as an argument  
#  sed "s/$search/$year/g" $templatename > $filename
#
#  # Submit the job 
#  sbatch $filename
   #python get-era5-hourly-local.py $year
   ./get-era5-monthly-single-levels.py $v $year
done
done
