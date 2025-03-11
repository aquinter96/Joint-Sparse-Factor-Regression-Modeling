#!/bin/bash

module unload r
module load r/4.1.0

for datset in {1..200}
do
for sampsize in 200 400 600 800
do

sbatch -p general -N 1 --mem=32g  -n 1 -t 168:00:00 --wrap="R CMD BATCH --vanilla --slave "--args $datset $sampsize" Fig4_lessSEM_Code.R

done
done