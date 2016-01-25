#CAGe
####Changepoint Analysis for Efficient Variant Calling

##Installation
####Prerequisites
* C++ compiler with support for C++11.
* Bamtools library.  If not already installed, follow directions [here](https://github.com/pezmaster31/bamtools/wiki).

####Installation Steps
1. Edit the first two lines in the `Makefile`.
  1. Change the C++ compiler, if necessary (g++ by default).
  2. Provide path to the bamtools top-level directory (this should have `include` and `lib` as subdirectories). Include the final `/`.
2. Run `make`. This will place compiled binaries in the directory `bin`.

##Usage
The basic CAGe pipeline consists of the following steps:

1. Dump the reference to stdout.
2. Dump the reads to stdout (using `bin/bamdump`).
3. Run changepoint detection to segment the genome (using `bin/cage` with piped stdin from steps 1 and 2).
4. Classify segments as high-complexity or low-complexity based on the inferred parameters of the segment (using `scripts/classify.py`), resulting in two `.bed` files.

These steps can be run with `scripts/run_cage.sh`. Edit the variable definitions at the top of that file to provide your inputs. 

In step 1, the reference must be provided such that the nth byte corresponds exactly to the nth base. To create such a file from a fasta file, use `scripts/dump_reference.sh`.
In step 4, `scripts/classify.py` classifies regions using simple thresholds, which must be provided in a parameter file. An example, appropriate for na12878, is provided in `scripts/config_example.txt`. It is important to set the coverage and error rate thresholds to match the statistics of the sample.

##Notes
The changepoint detection in CAGe depends on a rough estimate of SNPs and indels in the sample. These can either be provided to CAGe through an sqlite3 database, or CAGe can call variants internally using a simple built-in variant caller (referred to as CAGe++ in the paper)

1. To use the internal variant caller, run CAGe with the flag `-o <output_vcf>`.  This will cause CAGe to call variants and output them to the corresponding vcf file.  The internal variant caller is very naive and its output should probably not be used for any downstream analyses. Note that

  * Variants are called at each location using hard-coded thresholds of coverage, mismatches, and strand bias.
  * No meta-information about the sample is written to the VCF.
  * No attempt is made to estimate quality.
  * No attempt is made to estimate zygosity.

2. To specify an input sqlite3 database, run CAGe with the flag `-s <filename>`.  To convert a vcf file to an sqlite3 database with the proper tables, use the provided `scripts/vcf2sqlite3.sh`.

The flag `-v` tells CAGe to print verbose output while running. This is recommended. Some notes on the output:

* Information is provided every 100k base pairs.
* If Rdstarts, Matches, and Depth are always zero, this means no data is being read.
* If Matches and Depth are very different, this means that the reference sequence does not match the sample.
* If Max Costs is growing, this means CAGe is not running in linear time. This is unusual. In this case, just run CAGe with `beta` set at 0. Then a changepoint will be called every `step` base pairs (`step` of 100 or 200 should work fine). The classification of regions may still be useful.

####Citation for CAGe:
Bloniarz, A., Talwalkar, A., Terhorst, J., Jordan, M. I., Patterson, D., Yu, B., & Song, Y. S. (2014, January). Changepoint Analysis for Efficient Variant Calling. In *Research in Computational Molecular Biology* (pp. 20-34). Springer International Publishing.

CAGe implements the change-point detection algorithm introduced in:

Killick, R., Fearnhead, P., & Eckley, I. A. (2012). Optimal detection of changepoints with a linear computational cost. *Journal of the American Statistical Association*, 107(500), 1590-1598.
