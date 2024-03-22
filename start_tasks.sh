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


for i in $(seq 1 10); do
  j=$((i-1));
  s=$((i*20));
  dataset="ds_dataset${j}_${s}";
  for overlapper in "dgrid" "hagrid"; do
    for pipeline in "opt" "heuristic"; do
      for support_type in "steiner-tree" "path"; do
        for weight in "0" "0.5" "1"; do
          for run in $(seq 1 5); do
            jobname="esvis-$dataset-$weight-$support_type-$overlapper-$pipeline-$run";
            echo "submitting $jobname";
            qsub -N $jobname -l bc5 -l mem_free=$mem -l h_vmem=$mem -l h_rt=$hlr -e $base/logs/ -o $base/logs/ -r y run.sh $base/$jobname $support_type $overlapper $pipeline $dataset $weight
          done
        done
      done
    done
  done
done
