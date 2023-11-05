#!/bin/bash
AA_DATA_REPO=/home/data_repo
export AA_DATA_REPO

# Only needed if using Mosek8
# MOSEKLM_LICENSE_FILE=/home/mosek/
# export MOSEKLM_LICENSE_FILE

python2 /home/programs/AmpliconArchitect-master/src/AmpliconArchitect.py $argstring
