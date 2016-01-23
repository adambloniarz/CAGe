#!/bin/bash

# Convert VCF file to sqlite3 database in preparation for running cage
# Usage:
# ./vcf2sqlite3.sh <vcf> <out.db>

if [ "$#" -ne 2 ] 
then
    echo >&2 "Usage: $0 <vcf> <out.db>"
    exit 1
fi
rm -f $2
sqlite3 $2 "create table snps (chrom string, pos long, ref string, alt string); create index snps_index on snps (chrom, pos asc);"
egrep -v "^#" $1 | cut -f'1,2,4,5' -d$'\t' | sqlite3 -separator $'\t' $2 '.import /dev/stdin snps'
