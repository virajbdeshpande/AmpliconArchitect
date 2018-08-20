#!/bin/bash
echo $AA_DATA_REPO
python programs/AmpliconArchitect/src/AmpliconArchitect.py --bed /home/input/$BED_FILE --bam /home/input/$BAM_FILE --out /home/output/$OUT_NAME $OPTIONS