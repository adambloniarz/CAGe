#!/bin/sh

# Dump reference fasta file in preparation for running CAGe.
# Strips all whitespace and creates file such that the nth
# byte of the file is the nth base of the reference genome

# Usage:
# ./dump_reference.sh <INPUT_FASTA> <CHROMOSOME> <OUTPUT_FILE>
# Example:
# ./dump_reference.sh na12878.fa chr20 na12878_chr20.txt

samtools faidx $1 $2 | tail -n +2 | perl -p -e 's/\n//g' > $3
