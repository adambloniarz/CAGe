#!/anaconda/bin/python

# Usage:
# ./classify.py <CAGE_OUTPUT> <CONTIG> <BEDFILE_PREFIX> <PARAMS_FILE>
#
# Example:
# ./classify.py chr20_output.txt chr20 chr20_classified chr20_params.txt

import sys
from configparser import RawConfigParser

# hardcoded thresholds determining 'low' vs 'high' complexity regions
# thresholds determined via manual inspection of clusters (see R code)

# thresholds for the human chr20 data
#COV_MIN = 0.45
#COV_MAX = 0.7
#SEQ_ERR_MAX = 0.02
#MUT_MAX = 0.003

# thresholds for the Venter chr20 subset data
config = RawConfigParser()
config_file = sys.argv[4]
config.read(config_file)

COV_MIN = config.getfloat("cutoffs", "cov_min")
COV_MAX = config.getfloat("cutoffs", "cov_max")
ERR_MAX = config.getfloat("cutoffs", "err_max")
MUT_MAX = config.getfloat("cutoffs", "mut_max")
INDEL_MAX = config.getfloat("cutoffs", "indel_max")
ZETA_MAX = config.getfloat("cutoffs", "zeta_max")

def add_region(all_regions, start_pos, end_pos):
  if len(all_regions) > 0:
    end_pos_last = all_regions[-1][1]
    if end_pos_last + 1 == start_pos: # regions on either side of changepoint have same designation 
      all_regions[-1][1] = end_pos
    else:
      all_regions.append([start_pos, end_pos])
  else:
    all_regions.append([start_pos, end_pos])

def write_list_to_file(regions, chr_string, out_file):
  f_out = open(out_file, 'w')
  for region in regions:
    f_out.write("%s\t%s\t%s\n" % (chr_string, region[0], region[1]+1)) 
  f_out.close() 

#-------------------
# MAIN PROGRAM
#-------------------

f_in = open(sys.argv[1])
regions_low = []
regions_high = []

def is_low_complexity(coverage, error_rate, mut_rate, indel_rate, zeta, length):
  # return (coverage > COV_MIN and coverage < COV_MAX and seq_err < SEQ_ERR_MAX and mut_rate < MUT_MAX) or \
  return (coverage > COV_MIN and coverage < COV_MAX and mut_rate < MUT_MAX 
          and error_rate < ERR_MAX and indel_rate < INDEL_MAX and zeta < ZETA_MAX) or \
          ((length > 20000) and (coverage < 5e-4))

for line in f_in:
  new_region = line.split()
  start_pos = int(new_region[0])
  end_pos = start_pos + int(new_region[1]) - 1
  coverage, error_rate, mut_rate, indel_rate, zeta = map(float, new_region[2:])
  length = end_pos - start_pos
  if is_low_complexity(coverage, error_rate, mut_rate, indel_rate, zeta, length):
    add_region(regions_low, start_pos, end_pos) 
  else:
    add_region(regions_high, start_pos, end_pos) 
f_in.close() 

chr_string = sys.argv[2]
write_list_to_file(regions_low, chr_string, sys.argv[3] + '_low.bed')
write_list_to_file(regions_high, chr_string, sys.argv[3] + '_high.bed')

