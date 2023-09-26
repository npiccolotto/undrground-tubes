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

python $base/render2.py --read-dir $base/data --write-dir $TMPDIR

cp -r $TMPDIR $base/$1
