#!/bin/bash

# set the base directory
base='/home1/npiccolotto/ensemble-sets/results'

# set a job name
jobname=EnsembleSets

# set hard memory constraint
#hm=8G

# set required free memory
#rfm=8G

# set soft memory (not used)
#sm=7.8G

mem=16G

# set hard limit run time (seconds)
#hlr=144000
hlr=86400

qsub -N $jobname -l bc3 -l mem_free=$mem -l h_vmem=$mem -l h_rt=$hlr -e $base/logs/ -o $base/logs/ -r y run.sh $base/test