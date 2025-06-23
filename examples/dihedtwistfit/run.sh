#!/bin/bash
set -e
set -u

ffpopt-DihedTwistWorkflow.py -p minimal_orig.parm7 -c minimal_orig.rst7 --bond=9,8 --bond=8,6 --bond=6,2 --model=qdpi2 > scan.sh

bash ./scan.sh


