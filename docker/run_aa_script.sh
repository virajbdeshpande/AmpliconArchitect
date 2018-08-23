#!/bin/bash
export AA_DATA_REPO=/home/data_repo
python programs/AmpliconArchitect-master/src/AmpliconArchitect.py --bed $BED_FILE --bam $BAM_FILE --out $OUT_PREFIX $OPTIONS