#!/bin/bash

# Edit these variables to specify input data and output path
export ref_file=<PATH TO REFERENCE FILE>
export bam_file=<PATH TO BAM FILE>
export contig=<CONTIG, e.g. chr20>
export start=<START LOCATION, e.g. 60000>
export end=<END LOCATION, e.g. 63000000>
output_dir=cage_output
classify_config_file=config_example.txt

# Edit these variables to change parameters of changepoint detection
# Note: the defaults should work well for most scenarios
step=100
beta=3.0

mkdir -p $output_dir
( cat $ref_file; ../bin/bamdump $bam_file $contig $start $end ) | ../bin/cage $contig $start $end $step $beta $output_dir/$contig.out -v -o $output_dir/"$contig"_cage.vcf
./classify.py $output_dir/$contig.out $contig $output_dir/$contig $classify_config_file
