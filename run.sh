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

support_type=$2
support_partition=$3
layouter=$4
overlapper=$5
router=$6
dataset=$7

jobname="esvis-$dataset-$support_type-$support_partition-$layouter-$overlapper-$router"

python $base/render2.py --read-dir $base/data --write-dir $TMPDIR --support-type $support_type --support-partition $support_partition --layouter $layouter --overlap-remover $overlapper --router $router --read-dir $base/designspace-test --dataset $dataset
cp -r $TMPDIR $1
