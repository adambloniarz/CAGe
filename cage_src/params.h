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
#ifndef CAGE_PARAMS
#define CAGE_PARAMS

#include <string>

namespace cage {

  struct DataLoaderParams {
    std::string contig;
    int start;
    int end;
    int step;
    std::string output_VCF;
    std::string input_SNP_db;
    bool call_variants;
    bool verbose;
  };

  struct CageParams {
    DataLoaderParams dataloader_params;
    double beta;
    std::string output_file;
  };
}
#endif
