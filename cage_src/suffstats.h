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
#ifndef CAGE_SUFFSTATS
#define CAGE_SUFFSTATS

#include <vector>
#include "pileslist.h"

namespace cage {
  class SuffStats {

    public:
      int L;
      int N; // Read starts

      int tot_bases;
      int tot_errors;
      int N_SNPs;
      int indels;
      int zero_mapq;

      long inserts_sum;
      long long inserts_sumsq;
      int n_inserts;

      SuffStats(const Pileslist& piles, const int& start, const int& end, const int& zero_position);
      SuffStats(const std::vector<SuffStats>& all_stats, const int& start, const int& end);
      SuffStats();

      bool empty_region();
  };
}

#endif
