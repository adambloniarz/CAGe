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
#include <algorithm>
#include <map>
#include <vector>
#include <set>
#include <sys/time.h>
#include <time.h>
#include <math.h>
#include <thread>
#include <cstring>
#include "cost.h"
#include "dataloader.h"

using namespace std;

namespace cage {

  template <class It> bool all_zero(It start, It end) {
    for (It i = start; i < end; i++) {
      if (*i != 0) {
        return false;
      }
    }
    return true;
  }

  double auto_pick_beta(Dataloader &dload, const int& start, const int& end, const int& step) {

    // Spawn the first data loading thread
    thread t(&Dataloader::loadfunc, &dload);

    // Load the data for this chunk
    cout << "Calculating costs for determining beta\n";

    double total_cost = cost(dload.get_region(start, end));
    printf("Total cost of %i to %i is %.4f\n", start, end, total_cost);

    for (int j = step; j < (int)(end - start - step); j += step) {
      double cost1 = cost(dload.get_region(start, start + j));
      double cost2 = cost(dload.get_region(start + j, end));
      double sumcost = cost1 + cost2;
      printf("Cost with changepoint at %i:\t%.4f, gain from whole chunk: %.4f\n", j, sumcost, total_cost - sumcost);

    }

    // End the data-loading thread
    if(t.joinable())
    {
      cout << "Data loading thread terminated\n" << endl;
      t.join();
    }

    return 150.0;

  }

  void pelt(Dataloader &dload, const CageParams& params) {

    // Copy parameters for convenience
    int start = params.dataloader_params.start;
    int end = params.dataloader_params.end;
    int step = params.dataloader_params.step;
    double beta = params.beta;
    double K = 0;
    bool verbose = params.dataloader_params.verbose;
    string output_file = params.output_file;

    // Allocate vector for storing subproblem results
    vector<double> f((end - start) / step + 1, 0);

    f[0] = -1.0 * beta;

    // Store last prior changepoint positions for subproblems
    vector<int> cps((end - start) / step + 1, 0);
    cps[0] = start;

    vector<int> R(1);
    R[0] = start;

    int i;
    if (verbose) 
      cout << "Progress:\n";

    /*
       Note on indexing:
       i is in genome coordinates
       R stores numbers as genome coordinates
       */

    int maxcosts = 0;
    int maxlength = 0;
    long long int total_length = 0;

    // Spawn the data loading thread
    thread t(&Dataloader::loadfunc, &dload);

    // Loop through the data
    for (i = start; i < end + 1; i+= step) {

      if (beta == 0) {
        R[0] = i;
        continue;
      }

      maxcosts = (int)R.size() > maxcosts ? R.size() : maxcosts;
      int length_now = (i - *min_element(R.begin(), R.end()));
      maxlength = length_now > maxlength ? length_now : maxlength;

      if (verbose) {
        if ((i - start) % 100000 == 0) {
          printf("Position: %i\tMax Costs: %i\tMax Length: %i\tTotal length: %lld\n", i, maxcosts, maxlength, total_length);
          maxcosts = 0;
          maxlength = 0;
          total_length = 0;
        }
      }

      // Deal with the centromere - area of no coverage
      if ((i > start) & (i < end - 1)) {	
        if (dload.get_region(i, min(i + step, end)).empty_region()) {
          f[(i - start) / step] = f[(i - step - start) / step];
          cps[(i - start) / step] = cps[(i - step - start) / step];
          continue;
        }
      }

      vector<float> costs(R.size());
      vector<float> Fs(R.size());

      for (int j = 0; j < (int)R.size(); j++) {
        costs[j] = cost(dload.get_region(R[j], i));

        total_length += i - R[j];
        Fs[j] = f[(R[j]-start) / step];
      }

      vector<float> vals(costs.size());

      for (int j = 0; j < (int)vals.size(); j++)
        vals[j] = costs[j] + Fs[j];

      auto min_ptr = min_element(vals.begin(), vals.end());
      int idx = distance(vals.begin(), min_ptr);
      f[(i - start) / step] = *min_ptr + beta;
      cps[(i - start) / step] = R[idx];

      vector<int> R_new(0);
      for (int j = 0; j < (int)vals.size(); j++) {
        if (vals[j] + K <= f[(i - start) / step]) {
          R_new.push_back(R[j]);
        }
      }

      R_new.push_back(i);
      R = move(R_new);
    }

    i -= step;	

    vector<int> cps_n;

    if (beta != 0) {

      int pos;
      cps_n.push_back(end);
      cps_n.push_back(i);

      pos = cps[(i - start) / step];

      while (pos != start) {
        cps_n.push_back(pos);
        pos = cps[(pos - start) / step];
      }

      cps_n.push_back(start); 
      reverse(cps_n.begin(), cps_n.end());

      // Take the unique elements
      auto it = unique(cps_n.begin(), cps_n.end());
      cps_n.resize(distance(cps_n.begin(), it));

    } else {
      printf("Calling changepoints every %i bases\n", step);
      for (int k = start; k < (int)end; k += step) {
        cps_n.push_back(k);
      }
      cps_n.push_back(end);
    }

    cout << "Num of changepoints: " << cps_n.size() << endl;
    // Get the parameter values
    vector<double> lambdas(cps_n.size()-1);
    vector<double> eps(cps_n.size()-1);
    vector<double> gammas(cps_n.size()-1);
    vector<double> iotas(cps_n.size()-1);
    vector<double> zetas(cps_n.size()-1);
    vector<double> mus(cps_n.size()-1);
    vector<double> sigmas(cps_n.size()-1);
    vector<int> lengths(cps_n.size()-1);

    int row = 0;

    for (auto i = cps_n.begin(); i != cps_n.end() - 1; ++i) {
      int left = *i;
      int right = *(i+1);

      SuffStats tot_stats = dload.get_region(left, right);
      lambdas[row] = (double)tot_stats.N / tot_stats.L;
      eps[row] = (tot_stats.tot_bases == 0) ? 0 : (double)tot_stats.tot_errors / tot_stats.tot_bases;
      gammas[row] = (double)tot_stats.N_SNPs / tot_stats.L;
      iotas[row] = (tot_stats.tot_bases == 0) ? 0 : (double)tot_stats.indels / tot_stats.tot_bases;
      zetas[row] = (tot_stats.N == 0) ? 0 : (double)tot_stats.zero_mapq / tot_stats.N;
      mus[row] = (tot_stats.n_inserts == 0) ? 0 : (double)tot_stats.inserts_sum / tot_stats.n_inserts;
      sigmas[row] = (tot_stats.n_inserts > 0) ? 0 : sqrt((tot_stats.inserts_sumsq - pow((long long)tot_stats.inserts_sum,2)) / (tot_stats.n_inserts - 1));

      lengths[row] = right - left;
      row++;
    }

    if(t.joinable())
    {
      t.join();
    }

    FILE *f_out;
    if (strcmp(output_file.c_str(), "") != 0) { 
      f_out = fopen(output_file.c_str(), "w");
      if (f_out == nullptr) {
        throw runtime_error("Could not open file " + output_file); 
      }
      for (row = 0; row < (int)cps_n.size() - 2; row++) {
        if ((lengths[row] != 0)) { 
          fprintf(f_out, "%i\t%i\t%0.8f\t%0.8f\t%0.8f\t%0.8f\t%0.8f\n", cps_n[row], 
              lengths[row], lambdas[row], eps[row], gammas[row], iotas[row], zetas[row]);
        }
      }
      fclose(f_out); 
    }
  }
}
