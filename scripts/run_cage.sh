#!/bin/bash

# ----------------------------------------------------------
# Edit these variables to specify input data and output path
export REF_FILE=<REFERENCE FILE>
export BAM_FILE=<BAM FILE>
export CONTIG=<CONTIG>
export START=<START POSITION, eg 60000>
export END=<END POSITION, eg 63000000>
OUTPUT_DIR=cage_output
CLASSIFY_CONFIG_FILE=scripts/config_example.txt
# ----------------------------------------------------------

# ----------------------------------------------------------
# Edit these variables to change parameters of changepoint detection
# Note: the defaults should work well for most scenarios
STEP=100
BETA=3.0
# ----------------------------------------------------------

# ----------------------------------------------------------
export SCRIPT_PATH="`dirname \"$0\"`"
mkdir -p $OUTPUT_DIR
( cat $REF_FILE; $SCRIPT_PATH/../bin/bamdump $BAM_FILE $CONTIG $START $END ) | $SCRIPT_PATH/../bin/cage $CONTIG $START $END $STEP $BETA $OUTPUT_DIR/$CONTIG.out -v -o $OUTPUT_DIR/"$CONTIG"_cage.vcf
$SCRIPT_PATH/classify.py $OUTPUT_DIR/$CONTIG.out $CONTIG $OUTPUT_DIR/$CONTIG $CLASSIFY_CONFIG_FILE
# ----------------------------------------------------------
