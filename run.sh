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

support_type="$2"
overlapper="$3"
connecter="$4"
router="$5"
dataset="$6"
weight="$7"

jobname="esvis-$dataset-$weight-$support_type-$layouter-$overlapper-$connecter-$router";

cd $base
python render2.py --write-dir $TMPDIR --support-type $support_type --overlap-remover $overlapper --connecter $connecter --router $router --read-dir $base/designspace-test --dataset $dataset --weight $weight
cp -r $TMPDIR $1
