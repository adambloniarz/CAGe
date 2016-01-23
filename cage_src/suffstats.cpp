/**
 * Copyright 2015 Adam Bloniarz
 * 
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 *     http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
#include <stdlib.h>
#include "suffstats.h"

using namespace std;

namespace cage {
  SuffStats::SuffStats() : L(0), N(0), tot_bases(0), tot_errors(0), N_SNPs(0), indels(0), zero_mapq(0), inserts_sum(0), inserts_sumsq(0), n_inserts(0) {}

  SuffStats::SuffStats(const std::vector<SuffStats>& all_stats, const int& start, const int& end) {
    L = 0;
    N = 0;
    tot_bases = 0;
    tot_errors = 0;
    N_SNPs = 0;
    indels = 0;
    zero_mapq = 0;
    inserts_sum = 0;
    inserts_sumsq = 0;
    n_inserts = 0;

    for (int k = start; k < end; k++) {
      N += all_stats[k].N;  
      tot_bases += all_stats[k].tot_bases;  
      tot_errors += all_stats[k].tot_errors;  
      N_SNPs += all_stats[k].N_SNPs;  
      indels += all_stats[k].indels;  
      L += all_stats[k].L;  
      zero_mapq += all_stats[k].zero_mapq;
      n_inserts += all_stats[k].n_inserts;
      inserts_sum += all_stats[k].inserts_sum;
      inserts_sumsq += all_stats[k].inserts_sumsq;
    }
  }

  bool SuffStats::empty_region() {
    return (N == 0);
  }

  SuffStats::SuffStats(const Pileslist& piles, const int& start, const int& end, const int& zero_position) {

    L = end - start;
    N = 0;
    tot_bases = 0;
    tot_errors = 0;
    N_SNPs = 0;
    indels = 0;
    zero_mapq = 0;
    inserts_sum = 0;
    inserts_sumsq = 0;
    n_inserts = 0;

    int s;
    for ( int i = 0 ; i < (int)L; i++ ) {

      // Read starts
      s = piles.rdstarts[ i + zero_position ];	
      N += s;

      int bases, errors;
      if (piles.snp_indicators[i + zero_position] == 1) {
        N_SNPs++;
      }
      bases = piles.depth[i + zero_position];
      errors = bases - piles.matches[i + zero_position];
      indels += piles.indels[i + zero_position];
      tot_bases += bases;
      tot_errors += errors;
      zero_mapq += piles.zero_mapq[i + zero_position];

    }

    for (auto it = piles.inserts.lower_bound(start); it != piles.inserts.upper_bound(end); ++it) {
      for (int isize : it->second) {
        int asize = abs(isize);
        n_inserts++;
        inserts_sum += asize;
        inserts_sumsq += asize * asize; 
      }
    }
  }
}
