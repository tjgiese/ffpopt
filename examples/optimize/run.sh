#!/bin/bash
set -e
set -u

ffpopt-Optimze.py -p minimal_orig.parm7 -c minimal_orig.rst7 --constrain "10,9,8,6" --oscan oscan.xyz
