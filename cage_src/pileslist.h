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
#ifndef CAGE_PILESLIST
#define CAGE_PILESLIST

#include <vector>
#include <string>
#include <map>

namespace cage {
  // Struct of arrays for storing data
  class Pileslist {
    public:
      std::vector<int> rdstarts;
      std::vector<int> matches;
      std::vector<int> depth;
      std::vector<int> indels;
      std::string chr;
      std::map<int, std::vector<int>> inserts;

      std::vector<int> snp_indicators;
      std::vector<int> zero_mapq;

      Pileslist() {};
      Pileslist(int length, std::string chr);
      std::vector<int> insertsInRange(int start, int end);
  };
}

#endif
