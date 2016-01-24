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
#include <sstream>
#include <vector>
#include <iostream>
#include "alignment.h"

using namespace std;

namespace cage {
  vector<CigarOp> get_cigar_ops(string cigar) {

    size_t cig_idx = 0;
    vector<CigarOp> ops;
    ops.reserve(5);

    while (true) {
      size_t count_end = cigar.find_first_not_of("1234567890", cig_idx);
      if (count_end == string::npos) break;

      size_t count = stoi(cigar.substr(cig_idx, count_end - cig_idx));
      char op_type = cigar.at(count_end);
      ops.push_back(CigarOp(count, op_type));

      cig_idx = count_end + 1;
      if (cig_idx > cigar.size()) break;
    }
    return ops;

  }


  Alignment::Alignment(const string& al_string) {


    size_t al_pos, al_pos2, read_len;
    char delim = ' ';

    al_pos = al_string.find(delim, 0);
    Position = stoi(al_string.substr(0,al_pos));

    al_pos++;

    read_len = al_string.find(delim, al_pos) - al_pos;
    QueryBases = al_string.substr(al_pos, read_len);
    al_pos = al_pos + read_len + 1;

    al_pos2 = al_string.find(delim, al_pos);
    CigarData = get_cigar_ops(al_string.substr(al_pos, al_pos2 - al_pos));
    al_pos = al_pos2 + 1;

    Qualities = al_string.substr(al_pos, read_len);
    al_pos += read_len + 1;

    al_pos2 = al_string.find(delim, al_pos);
    MapQuality = stoi(al_string.substr(al_pos, al_pos2 - al_pos));
    al_pos = al_pos2 + 1;

    IsReverseStrand = stoi(al_string.substr(al_pos, 1)) == 1;


  }    
}
