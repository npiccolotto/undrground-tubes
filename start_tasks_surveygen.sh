#!/bin/bash

# set the base directory
base='/home1/npiccolotto/ensemble-sets/results'

# set hard memory constraint
#hm=8G

# set required free memory
#rfm=8G

# set soft memory (not used)
#sm=7.8G

mem=64G

# set hard limit run time (seconds)
hlr=25200

for i in $(seq 1 4); do
  dataset="ds_dataset${i}";
  for strategy in "heuristic" "opt"; do
    for weight in "0" "0.5" "1"; do
      jobname="essurvey-$dataset-$strategy-$weight";
      echo "submitting $jobname";
      qsub -N $jobname -l bc4 -l mem_free=$mem -l h_vmem=$mem -l h_rt=$hlr -e $base/logs/ -o $base/logs/ -r y run_surveygen.sh $base/$jobname $dataset $strategy $weight
    done
  done
done
