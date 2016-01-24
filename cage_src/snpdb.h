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
#ifndef CAGE_SNPDB
#define CAGE_SNPDB

#include <iostream>
#include <string>
#include <sqlite3.h>


namespace cage {
  class SNPdb {
    public:
      SNPdb(const std::string& SNP_db_file);
      std::vector<std::tuple<int, int, int>> report_variants(const std::string& contig, const int& start, const int& end);

      // Does nothing because we are using an existing dataset of SNPs
      void increment_mismatch(const int& pos, const std::string& alt_bases, const bool& reverse_strand) {}

      // Do not allow copying of SNPdb objects
      SNPdb(SNPdb const&) = delete; 
      void operator=(SNPdb const&) = delete;

    private:
      sqlite3* db;
      sqlite3_stmt *intervalSearch;
  };
}
#endif
