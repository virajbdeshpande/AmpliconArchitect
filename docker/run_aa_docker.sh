#!/bin/bash
INPUT_DIR=$1
BAM_FILE=$2
BED_FILE=$3
OUT_DIR=$4
OUT_NAME=$5
OPTIONS=$6

#docker run --rm -it -e OPTIONS="$OPTIONS" -e BED_FILE=$BED_FILE -e BAM_FILE=$BAM_FILE -e OUT_NAME=$OUT_NAME -v $AA_DATA_REPO:/home/data_repo -v $INPUT_DIR:/home/input -v $OUT_DIR:/home/output aa sh /home/run_aa_script.sh

docker run --rm -it -e AA_DATA_REPO=/home/data_repo -e OPTIONS="$OPTIONS" -e BED_FILE=$BED_FILE -e BAM_FILE=$BAM_FILE -e OUT_NAME=$OUT_NAME -v $AA_DATA_REPO:/home/data_repo -v $INPUT_DIR:/home/input -v $OUT_DIR:/home/output aa sh /home/run_aa_script.sh