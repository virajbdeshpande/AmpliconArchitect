#!/bin/bash
AA_DATA_REPO=/home/data_repo
export AA_DATA_REPO
ls $AA_DATA_REPO
export MOSEKLM_LICENSE_FILE=$PWD/mosek/8/licenses

python programs/AmpliconArchitect-master/src/AmpliconArchitect.py --bed $BED_FILE --bam $BAM_FILE --out $OUT_PREFIX $OPTIONS