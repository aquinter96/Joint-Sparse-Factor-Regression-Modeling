#!/bin/bash

module unload r
module load r/4.3.2

array=( "Australia" "Austria" "Brazil" "Canada" "Croatia" "Germany" "Greece" "Indonesia" "Ireland" "Malaysia" "Netherlands" "Pakistan" "Panama" "Philippines" "Poland" "Portugal" "Serbia" "Slovakia" "South_Korea" "Spain" "Switzerland" "Turkey" "United_Kingdom" "United_States" )
for i in "${array[@]}"
do

sbatch -p general -N 1 --mem=32g --cpus-per-task=75 -n 1 -t 100:00:00 R CMD BATCH --vanilla --slave "--args $i" Real_Data_lavaan.R

done