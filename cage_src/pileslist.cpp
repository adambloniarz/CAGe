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
#include "pileslist.h"

using namespace std;

namespace cage {
  Pileslist::Pileslist( int length, string chr ) :
    rdstarts(length, 0), matches(length, 0), depth(length, 0), 
    indels(length, 0), chr(chr), snp_indicators(length, 0), zero_mapq(length, 0) {}

  vector<int> Pileslist::insertsInRange(int start, int end) {
    vector<int> ret;
    vector<pair<int, vector<int>>> b(inserts.lower_bound(start), inserts.upper_bound(end - 1));
    for (pair<int, vector<int>> p : b) {
      ret.insert(ret.end(), p.second.begin(), p.second.end());
    }
    return ret;
  }
}
