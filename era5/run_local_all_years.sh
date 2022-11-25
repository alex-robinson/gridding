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

# Variables 0..23 are initial variables down through "total_column_water_vapour"
# Variables 24..34 are additional variables of interest 

#for v in {0..23}
for v in {24..34}
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

      #echo $v $year
   done
done
