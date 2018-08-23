#!/bin/bash
echo $AA_DATA_REPO
python programs/AmpliconArchitect-master/src/AmpliconArchitect.py --bed $BED_FILE --bam $BAM_FILE --out $OUT_PREFIX $OPTIONS