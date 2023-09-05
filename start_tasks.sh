#!/bin/bash

# set the base directory
base='/home1/mwallinger/epb/path-bundling/epb_extension/code'

# set a job name
jobname=EnsembleSets

# set hard memory constraint
#hm=8G

# set required free memory
#rfm=8G

# set soft memory (not used)
#sm=7.8G

# set which queue is used to execute the job
cl=bc4

# set hard limit run time (seconds)
#hlr=144000
hlr=86400

  # for exp in 1 2 3 4 5 6 9 10 11 12 13
  # do
  #   qsub -N $jobname -l $cl -l mem_free=32G -l h_vmem=32G -l h_rt=$hlr -e $base/logs/ -o $base/logs/ -r y run.sh $exp
  #   sleep 5
  # done;

  # for exp in 8 15 16 17
  # do
  #   qsub -N $jobname -l $cl -l mem_free=64G -l h_vmem=64G -l h_rt=$hlr -e $base/logs/ -o $base/logs/ -r y run.sh $exp
  #   sleep 5
  # done;

  # for exp in 8 18
  # do
  #   qsub -N $jobname -l $cl -l mem_free=64G -l h_vmem=64G -l h_rt=$hlr -e $base/logs/ -o $base/logs/ -r y run.sh $exp
  #   sleep 5
  # done;

  # for exp in 22 23
  # do
  #   qsub -N $jobname -l $cl -l mem_free=32G -l h_vmem=32G -l h_rt=$hlr -e $base/logs/ -o $base/logs/ -r y run.sh $exp
  #   sleep 5
  # done;

  for exp in 24
  do
    qsub -N $jobname -l $cl -l mem_free=96G -l h_vmem=96G -l h_rt=$hlr -e $base/logs/ -o $base/logs/ -r y run.sh $exp
  done;
