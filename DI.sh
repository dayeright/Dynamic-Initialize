#!/bin/bash

echo "Please remember check and change di_info.txt for the information of typhoon that been investigated  "
echo "The first two lines are location of (t0-6)h (lat,lon)"
echo "The following three lines are location and strength of (t0)h (lat,lon,strength)"
source /etc/profile && yhq

(cat namelist.input | grep start_ | cut -f2 -d",")>di_time.txt

# The first step to relocate the center of the typhoon
ncl step1_judge.ncl

# First run 
yhbatch -n 120 -p TH_SR1 -J deyu_di ./sub.sh 
judge=$(yhq | grep "deyu_di")
echo $judge
while [ -n "$judge" ]
do 
sleep 2m
judge=$(yhq | grep "deyu_di")
done 

# Judge 
ncl step0_windloc.ncl 
differ1=$(cat di_judge.txt | tr -cd "[0-9],-,.")
differ=${differ1%.*}
echo $differ
n=1

while [[ $differ -le -2 || $n -lt 2 ]]
do 
ncl step1c_judge.ncl 

# Substitute 
ncl step2_substitute.ncl 

# Cycle run 
yhbatch -n 120 -p TH_SR1 -J deyu_di ./sub.sh 
judge=$(yhq | grep "deyu_di")
while [ -n "$judge" ]
do 
sleep 2m
judge=$(yhq | grep "deyu_di")
done 
ncl step0_windloc.ncl 
differ1=$(cat di_judge.txt | tr -cd "[0-9],-,.")
differ=${differ1%.*}
n=`expr $n + 1`
done 

ncl step1c_judge.ncl 
sed -i '3c  run_hours                           = 0,' namelist.input
sed -i '9c  start_hour                          = 12,   12,   12,' namelist.input 
sed -i '80,99d' namelist.input 
sed -i '1c 19.5' bogusvortex.txt 
sed -i '2c 128.6' bogusvortex.txt 

./real.exe
ncl step2_substitute.ncl
yhbatch -n 120 -p TH_SR1 -J deyu_dif ./sub.sh

echo "DI has been completed!"
