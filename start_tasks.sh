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
hlr=46800 # 13 hours. 12h are all the ILPs if we're unlucky + 1h for the rest of the code.


for i in $(seq 1 10); do
  j=$((i-1));
  s=$((i*20));
  dataset="ds_dataset${j}_${s}";
  for layouter in "mds" "qsap"; do
    for overlapper in "dgrid" "hagrid"; do
      for connecter in "opt" "heuristic"; do
        for router in "opt" "heuristic"; do
          for support_type in "steiner-tree" "path"; do
            jobname="esvis-$dataset-$support_type-$layouter-$overlapper-$connecter-$router";
            echo "submitting $jobname";
            qsub -N $jobname -l bc3 -l longrun=1 -l mem_free=$mem -l h_vmem=$mem -l h_rt=$hlr -e $base/logs/ -o $base/logs/ -r y run.sh $base/$jobname $support_type $layouter $overlapper $connecter $router $dataset
          done
        done
      done
    done
  done
done

