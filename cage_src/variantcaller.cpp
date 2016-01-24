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
#include <fstream>
#include <sstream>
#include <iterator>
#include "variantcaller.h"

using namespace std;

static int depth_min = 10;
static double var_depth_frac_min = .20;
static int var_depth_hard_min = 3;
static double strand_bias_min = 0.1;

namespace cage {

  string merge_string_vector(vector<string> strings) {
    const char* const delim = ",";

    if (strings.size() == 0)
      return(string());

    if (strings.size() == 1)
      return strings[0];

    ostringstream imploded;
    copy(strings.begin(), strings.end(), ostream_iterator<string>(imploded, delim));

    string ret = imploded.str();
    return(ret.substr(0, ret.size() - 1));

  }

  VariantCaller::VariantCaller(const int& length, const int& start_position, const string& full_ref_, const string& contig) :
    start_position(start_position), full_ref(full_ref_), contig(contig) {
      piles = vector<unique_ptr<alt_map>>(length);
    }

  void VariantCaller::increment_mismatch(const int& pos, const string& alt_bases, const bool& reverse_strand) {
    if (!piles[pos]) {
      piles[pos] = unique_ptr<alt_map>(new alt_map());
    }

    if (piles[pos]->find(alt_bases) == piles[pos]->end()) {
      if (reverse_strand)
        piles[pos]->insert(make_pair(alt_bases, alt_counts(0,1)));
      else
        piles[pos]->insert(make_pair(alt_bases, alt_counts(1,0)));
    } else {
      if (reverse_strand)
        piles[pos]->at(alt_bases).nr++;
      else
        piles[pos]->at(alt_bases).nf++;
    }
  }

vector<tuple<int, int, int>> VariantCaller::report_variants(const vector<int>& depths, const string& vcf_file) {
  vector<tuple<int, int, int>> positions;

  //TODO: deal with possible exception
  ofstream vcf_stream(vcf_file, ios::out | ios::app);

  for (int i = 0; i < (int)piles.size(); i++) {

    // No variants to consider
    if (!piles[i])
      continue;

    const map<string, alt_counts>& pile = *(piles[i].get());
    string ref = full_ref.substr(i + start_position, 1);

    vector<string> final_alts;
    vector<string> final_refs;

    if (depths[i] >= depth_min) {
      int var_depth_min = depths[i] * var_depth_frac_min;
      int max_altlen = 0;
      int max_reflen = 0;
      string ref_del;

      for (auto& alt : pile) {
        // TODO: check if this is necessary
        if (alt.first.size() == 0)
          continue;

        int count = alt.second.nf + alt.second.nr;

        if (min((double)alt.second.nf / count, (double)alt.second.nr / count) > strand_bias_min && 
            (count > max(var_depth_min, var_depth_hard_min)) &&
            (alt.first.compare(ref) != 0)) {

          // Look for a deletion
          if (alt.first.substr(0,1).compare("-") == 0) {
            if (i + alt.first.size() + 1 < piles.size()) {
              final_alts.push_back(ref);
              string ref_str = full_ref.substr(start_position + i, alt.first.size() + 1);
              final_refs.push_back(ref_str);
              max_reflen = max(max_reflen, (int)alt.first.size() + 1);
              max_altlen = max(max_altlen, 1);
              vcf_stream << contig << '\t' << i + start_position + 1 << "\t.\t" 
                << ref_str << '\t' << ref << "\t20\tPASS\t.\tGT\t1/1" << endl;
            }
          } else {
            final_alts.push_back(alt.first);
            final_refs.push_back(ref);
            max_reflen = max(max_reflen, 1);
            max_altlen = max(max_altlen, (int)alt.first.size());
            vcf_stream << contig << '\t' << i + start_position + 1 << "\t.\t" 
              << ref << '\t' << alt.first << "\t20\tPASS\t.\tGT\t1/1" << endl;
          }
        }
      }

      if (final_alts.size() == 0)
        continue;

      // Send call back to cagar where there is only one variant
      if (final_alts.size() == 1)
        positions.push_back(make_tuple(i + start_position, max_reflen, max_altlen));
    }
  }
  vcf_stream.close();
  return positions;
}
}
