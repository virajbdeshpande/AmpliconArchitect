#!/bin/bash
AA_DATA_REPO=/home/data_repo
export AA_DATA_REPO
MOSEKLM_LICENSE_FILE=/home/programs/mosek/8/licenses
export MOSEKLM_LICENSE_FILE

python2 /home/programs/AmpliconArchitect-master/src/AmpliconArchitect.py $argstring
