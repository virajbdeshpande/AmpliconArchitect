#!/bin/bash
echo $AA_DATA_REPO
python programs/AmpliconArchitect-master/src/AmpliconArchitect.py --bed /home/input/$BED_FILENAME --bam /home/input/$BAM_FILENAME --out /home/output/$OUT_NAME $OPTIONS