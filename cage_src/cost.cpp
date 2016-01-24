/**
 * Copyright 2015 Adam Bloniarz, Jonathan Terhorst, Ameet Talwalkar
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
#include <iostream>
#include <sstream>
#include <cmath>
#include <limits>
#include "cost.h"

#define PI 3.1415926535

using namespace std;

namespace cage {
  inline double binom_loglik(int N, int K,  double p) {
    if (p <= 0.0 || p >= 1.0) {
      return 0;
    } else {
      return(K * log(p) + (N - K) * log(1.0 - p));
    }
  }

  double cost(const SuffStats& stats) {

    if (stats.L == 0) {
      return 0;
    }

    long L = stats.L;

    /*
       -------------------------
       1) Read starts - lambda
       -------------------------
       */

    double llik_starts;

    // Deal with the case where there are no read starts at all
    if (stats.N == 0) {
      llik_starts = 0;
    } else {
      llik_starts = -(double)stats.N + stats.N * log((double)stats.N / L);
    }

    /*
       -------------------------
       2) Mutations - gamma
       3) Errors - epsilon
       -------------------------
       */

    // Mutations
    double gamma = (double)stats.N_SNPs / L;
    double llik_SNPs = binom_loglik(L, stats.N_SNPs, gamma);

    // Errors
    double epsilon = (stats.tot_bases == 0) ? 0.0 : (double)stats.tot_errors / stats.tot_bases;
    double llik_errors = binom_loglik(stats.tot_bases, stats.tot_errors, epsilon);

    // Indels
    double iota = (stats.tot_bases == 0) ? 0.0 : (double)stats.indels / stats.tot_bases;
    double llik_indels = binom_loglik(stats.tot_bases, stats.indels, iota);

    // Zero mapq
    double zeta = (stats.N == 0) ? 0.0 : (double)stats.zero_mapq / stats.N;
    double llik_zero_mapq = binom_loglik(stats.N, stats.zero_mapq, zeta);

    /*
       -------------------------
       4) Insert size calcs
       -------------------------
       */
    // loglik is (without constants which don't matter): 
    //
    //      -0.5 * (pow(x - mu, 2) / sigmasq + log(sigmasq));
    //      
    // loglik for all points when mu, sigmasq are set to MLE
    //      
    //      sum_x [ -0.5 * (pow(x - mu, 2) / sigmasq + log(sigmasq)) ]
    //    
    //    = -0.5 * ( sum_x [ pow(x - mu, 2) / sigmasq ] + n * log(sigmasq) )
    //    
    //    = -0.5 * ( n - 1 + n * log(sigmasq) )
    //
    //    = -0.5 * n * (1 + log(sigmasq)), again dropping constants.

    double llik_inserts;
    if (stats.n_inserts > 1) {
      double mu = (double)stats.inserts_sum / stats.n_inserts;
      // cout << stats.inserts_sumsq << " " << stats.inserts_sum << " " << stats.n_inserts << endl;
      double sigmasq = ((double)stats.inserts_sumsq - pow(mu, 2)) / (stats.n_inserts - 1);
      llik_inserts = -0.5 * stats.n_inserts * (1.0 + log(sigmasq));
    } else {
      llik_inserts = 0;
    }

    // Note: llik_inserts is not used
    
    return -1.0*(llik_starts + llik_SNPs + llik_indels + llik_errors + llik_zero_mapq);
  }
}
