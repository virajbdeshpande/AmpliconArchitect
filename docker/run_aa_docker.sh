#!/bin/bash

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
        BAM_ARG=$(printf %q "/home/bam_dir/$BAM_FILENAME")
        argstring="$argstring $BAM_ARG"
    elif [ "$bed_set" == "True" ]
    then
        bed_set=False
        BED_FILE=`readlink -f "$i"`
        BED_DIR=`dirname "$BED_FILE"`
        BED_FILENAME=`basename "$BED_FILE"`
        BED_ARG=$(printf %q "/home/bed_dir/$BED_FILENAME")
        argstring="$argstring $BED_ARG"
    elif [ "$out_set" == "True" ]
    then
        out_set=False
        OUT_PREFIX=`readlink -f "$i"`
        OUT_DIR=`dirname "$OUT_PREFIX"`
        OUT_NAME=`basename "$OUT_PREFIX"`
        OUT_ARG=$(printf %q "/home/output/$OUT_NAME")
        # OUT_ARG=\"/home/output/$OUT_NAME\"
        #echo "$OUT_ARG"
        #echo $argstring
        argstring="$argstring $OUT_ARG"
        #echo $argstring
        #echo "$argstring"
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

docker run --rm -e AA_DATA_REPO=/home/data_repo -e argstring="$argstring" -v $AA_DATA_REPO:/home/data_repo -v $BAM_DIR:/home/bam_dir -v $BED_DIR:/home/bed_dir -v $OUT_DIR:/home/output -v $MOSEKLM_LICENSE_FILE:/home/programs/mosek/8/licenses virajbdeshpande/ampliconarchitect bash /home/run_aa_script.sh
