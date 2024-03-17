#!/bin/bash

source /home1/npiccolotto/.bashrc

base='/home1/npiccolotto/ensemble-sets/code'
# adapt for local use
#base='/Users/npiccolotto/Projects/cvast/bssvis/ensemble-set-rendering'
#export TMPDIR='/tmp/garbage'

# test that all the binaries exist: pyenv, loom
if ! command -v pyenv 1>/dev/null 2>&1; then
  echo "pyenv not installed, can't continue." >&2
  exit 1
fi

if ! command -v loom 1>/dev/null 2>&1; then
  echo "loom not installed, can't continue." >&2
  exit 1
fi

# set python version
pyenv local 3.10

# set gurobi stuff
export GRB_LICENSE_FILE="/home1/share/gurobi/gurobi.lic"
export GUROBI_HOME="/home1/share/gurobi/gurobi/linux64"

dataset="$2"
strategy="$3"
weight="$4"

jobname="essurvey-$dataset-$strategy-$weight";

cd $base
python render2.py --strategy $strategy --write-dir $TMPDIR --overlap-remover dgrid --connecter auto --router auto --read-dir $base/eval/dataset --dataset $dataset --weight $weight --connect-objective joint --serialize-graph 0 --compute-metrics 0 -w 18 -h 18
cp -r $TMPDIR $1
