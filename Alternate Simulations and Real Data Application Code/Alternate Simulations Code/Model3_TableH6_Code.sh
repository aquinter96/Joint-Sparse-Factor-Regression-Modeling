#!/bin/bash

module unload r
module load r/4.1.0

for datset in {1..500}
do
for sampsize in 200 400 600 800
do

sbatch -p general -N 1 --mem=32g --cpus-per-task=25 -n 1 -t 50:00:00 --wrap="R CMD BATCH --vanilla --slave "--args $datset $sampsize" Model3_Simulation_Code.R

done
done