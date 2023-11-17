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
layouter="$3"
overlapper="$4"
connecter="$5"
router="$6"
dataset="$7"

jobname="esvis-$dataset-$support_type-$layouter-$overlapper-$connecter-$router";

cd $base
python render2.py --write-dir $TMPDIR --support-type $support_type --layouter $layouter --overlap-remover $overlapper --connecter $connecter --router $router --read-dir $base/designspace-test --dataset $dataset
cp -r $TMPDIR $1
