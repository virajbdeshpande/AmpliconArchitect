#!/bin/bash
BAM_FILE=`readlink -f $1`
BED_FILE=`readlink -f $2`
OUT_PREFIX=`readlink -f $3`
OPTIONS=$4
BAM_DIR=`dirname $BAM_FILE`
BED_DIR=`dirname $BED_FILE`
OUT_DIR=`dirname $OUT_PREFIX`
BAM_FILENAME=`basename $BAM_FILE`
BED_FILENAME=`basename $BED_FILE`
OUT_NAME=`basename $OUT_PREFIX`

echo $BAM_DIR $BED_DIR $OUT_DIR


docker run --rm -it -e AA_DATA_REPO=/home/data_repo -e OPTIONS="$OPTIONS" -e BED_FILENAME=$BED_FILENAME -e BAM_FILENAME=$BAM_FILENAME -e OUT_NAME=$OUT_NAME -v $AA_DATA_REPO:/home/data_repo -v $BAM_DIR:/home/bam_dir -v $BED_DIR:/home/bed_dir -v $OUT_DIR:/home/output aa sh /home/run_aa_script.sh