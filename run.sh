#!/bin/bash

base='/home1/npiccolotto/ensemble-sets'
source $base/.venv/bin/activate

python -V
echo $1

python $base/code/render2.py --write-dir $TMPDIR --read-dir $base/code/data

cp -r $TMPDIR /home1/npiccolotto/ensemble-sets/results/$1
