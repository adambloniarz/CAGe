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
#ifndef CAGE_ALIGNMENT_H
#define CAGE_ALIGNMENT_H

#include <string>
#include <vector>

namespace cage {
  class CigarOp {
    public:
      int Length;
      char Type;
      CigarOp (int length, char type) : Length(length), Type(type) {}

  };

  class Alignment {
    public:
      int Position;
      std::string QueryBases;
      std::vector<CigarOp> CigarData;
      std::string Qualities;
      bool IsReverseStrand;
      int MapQuality;

      Alignment(const std::string& al_string);
  };
}
#endif
