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
#ifndef CAGE_DATALOADER
#define CAGE_DATALOADER

#include <string>
#include <mutex>
#include <vector>
#include <condition_variable>
#include <memory>
#include "suffstats.h"
#include "params.h"
#include "snpdb.h"

namespace cage {

  class Dataloader {
    public:
      Dataloader(DataLoaderParams params, std::string&& ref_);

      // Forbid copying or assigning
      Dataloader(Dataloader const&)=delete;
      Dataloader& operator=(Dataloader const&)=delete;

      // Functions for data access
      SuffStats get_region(const int& seg_start, const int& seg_end);

      void loadfunc();

    private:
      cage::DataLoaderParams params;

      // Reference string
      std::string ref;

      // Determines the number of base pairs
      // to load at one time
      // Hard-coded for now
      int load_step = 1000000;

      // For keeping track of current data position
      int dataend;

      // Mutexes and condition variables for data loading / reading
      std::mutex loading_lock;
      std::condition_variable load_request;

      // Data storage
      std::vector<SuffStats> all_stats;

      // Variant calling
      std::unique_ptr<SNPdb> snp_db;
  };
}
#endif
