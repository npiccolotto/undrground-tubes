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
#hlr=144000
hlr=14400


for support_type in "steiner-tree" "path"; do
  for support_partition in "intersection-group" "set"; do
    for layouter in "mds" "qsap"; do
      for overlapper in "dgrid" "hagrid"; do
        for router in "opt" "heuristic"; do
          for i in $(seq 1 10); do
            j=$((i-1));
            s=$((i*20));
            dataset="dataset${j}_${s}";
            jobname="esvis-$dataset-$support_type-$support_partition-$layouter-$overlapper-$router";
            echo "submitting $jobname";
            qsub -N $jobname -l bc3 -l mem_free=$mem -l h_vmem=$mem -l h_rt=$hlr -e $base/logs/ -o $base/logs/ -r y run.sh $base/$jobname $support_type $support_partition $layouter $overlapper $router $dataset
          done
        done
      done
    done
  done
done

