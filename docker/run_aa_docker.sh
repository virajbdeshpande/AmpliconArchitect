#!/bin/bash
# BAM_FILE=`readlink -f "$1"`
# BED_FILE=`readlink -f "$2"`
# OUT_PREFIX=`readlink -f "$3"`
# OPTIONS=$4
# BAM_DIR=`dirname "$BAM_FILE"`
# BED_DIR=`dirname "$BED_FILE"`
# OUT_DIR=`dirname "$OUT_PREFIX"`
# BAM_FILENAME=`basename "$BAM_FILE"`
# BED_FILENAME=`basename "$BED_FILE"`
# OUT_NAME=`basename "$OUT_PREFIX"`

bam_set=False
bed_set=False
out_set=False
argstring=""
for i in "$@"; do
    if [ "$bam_set" == "True" ]
    then
        bam_set=False
        BAM_FILE=`readlink -f "$i"`
        BAM_DIR=`dirname "$BAM_FILE"`
        BAM_FILENAME=`basename "$BAM_FILE"`
        argstring="$argstring /home/bam_dir/$BAM_FILENAME"
    elif [ "$bed_set" == "True" ]
    then
        bed_set=False
        BED_FILE=`readlink -f "$i"`
        BED_DIR=`dirname "$BED_FILE"`
        BED_FILENAME=`basename "$BED_FILE"`
        argstring="$argstring /home/bed_dir/$BED_FILENAME"
    elif [ "$out_set" == "True" ]
    then
        out_set=False
        OUT_PREFIX=`readlink -f "$i"`
        OUT_DIR=`dirname "$OUT_PREFIX"`
        OUT_NAME=`basename "$OUT_PREFIX"`
        argstring="$argstring /home/output/$OUT_NAME"
    elif [ "$i" == "--bam" ]
    then
        bam_set=True
        argstring="$argstring $i"
    elif [ "$i" == "--bed" ]
    then
        bed_set=True
        argstring="$argstring $i"
    elif [ "$i" == "--out" ]
    then
        out_set=True
        argstring="$argstring $i"
    else
        argstring="$argstring $i"
    fi
done

echo $argstring


# docker run --rm -it -e AA_DATA_REPO=/home/data_repo -e OPTIONS="$OPTIONS" -e BED_FILE="/home/bed_dir/$BED_FILENAME" -e BAM_FILE="/home/bam_dir/$BAM_FILENAME" -e OUT_PREFIX="/home/output/$OUT_NAME" -v $AA_DATA_REPO:/home/data_repo -v $BAM_DIR:/home/bam_dir -v $BED_DIR:/home/bed_dir -v $OUT_DIR:/home/output aa sh /home/run_aa_script.sh -v $MOSEKLM_LICENSE_FILE:/home/programs/mosek/8/licenses


docker run --rm -it -e AA_DATA_REPO=/home/data_repo -e argstring="$argstring" -v $AA_DATA_REPO:/home/data_repo -v $BAM_DIR:/home/bam_dir -v $BED_DIR:/home/bed_dir -v $OUT_DIR:/home/output aa sh /home/run_aa_script.sh -v $MOSEKLM_LICENSE_FILE:/home/programs/mosek/8/licenses