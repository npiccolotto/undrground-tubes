#!/bin/bash

# set the base directory
base='/home1/npiccolotto/ensemble-sets/results'

# set hard memory constraint
#hm=8G

# set required free memory
#rfm=8G

# set soft memory (not used)
#sm=7.8G

mem=32G

# set hard limit run time (seconds)
hlr=18000 # 4 calls to gurobi, 1h for each plus one hour for the rest = 5h


for i in $(seq 1 24); do
  j=$((i-1));
  dataset="ds_dataset${j}";
  for strategy in "heuristic" "opt"; do
    for weight in "0" "0.5" "1"; do
      jobname="essurvey-$dataset-$strategy-$weight";
      echo "submitting $jobname";
      qsub -N $jobname -l bc4 -l mem_free=$mem -l h_vmem=$mem -l h_rt=$hlr -e $base/logs/ -o $base/logs/ -r y run_surveygen.sh $base/$jobname $dataset $strategy $weight
    done
  done
done

