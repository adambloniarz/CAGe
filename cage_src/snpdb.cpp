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
#include <vector>
#include <mutex>
#include <iostream>
#include <sqlite3.h>
#include "snpdb.h"

using namespace std;

namespace cage {

  int max_token_length(string s, string delim) {
    int maxlen = 1;

    size_t pos = 0;
    string token;
    while (true) {
      pos = s.find(delim);
      token = s.substr(0, pos);

      maxlen = max(maxlen, (int)token.length());
      if (pos == string::npos) {
        break;
      } else {
        s.erase(0, pos + delim.length());
      }
    }
    return maxlen;
  }

  vector<tuple<int, int, int>> SNPdb::report_variants(const string& contig, const int& start, const int& end) {
    int result;
    sqlite3_reset(intervalSearch);
    vector<tuple<int, int, int>> ret;

    int rc;

    rc = sqlite3_bind_text(intervalSearch, 1, contig.c_str(), -1, SQLITE_STATIC);
    if(rc != SQLITE_OK) {
      sqlite3_close(db);
      throw runtime_error(string("Can't bind parameters: ") + string(sqlite3_errmsg(db)));
    }

    rc = sqlite3_bind_int(intervalSearch, 2, start + 1);
    if(rc != SQLITE_OK) {
      sqlite3_close(db);
      throw runtime_error(string("Can't bind parameters: ") + string(sqlite3_errmsg(db)));
    }

    rc = sqlite3_bind_int(intervalSearch, 3, end + 1);
    if(rc != SQLITE_OK) {
      sqlite3_close(db);
      throw runtime_error(string("Can't bind parameters: ") + string(sqlite3_errmsg(db)));
    }

    while(true) {
      result = sqlite3_step(intervalSearch);
      if (result == SQLITE_ROW) {
        int result_pos = sqlite3_column_int(intervalSearch, 0);
        string ref = string(reinterpret_cast<const char*>(sqlite3_column_text(intervalSearch, 1)));
        string alt = string(reinterpret_cast<const char*>(sqlite3_column_text(intervalSearch, 2)));

        // Determine the number of bases to ignore
        int ref_len = max_token_length(ref, ",");
        int alt_len = max_token_length(alt, ",");

        // Subtract one to adjust for 0-based indexing in CAGe
        ret.push_back(make_tuple(result_pos - 1, ref_len, alt_len));		
      } 
      else {
        break;
      }
    }
    //mtx.unlock();
    return ret;
  }

  SNPdb::SNPdb(const string& SNP_db_filename) {
    int rc;
    sqlite3 *db1;

    rc = sqlite3_open(":memory:", &db);
    if(rc != SQLITE_OK) {
      sqlite3_close(db);
      throw runtime_error(string("Can't open database: ") + string(sqlite3_errmsg(db)));
    }

    rc = sqlite3_open(SNP_db_filename.c_str(), &db1);
    if(rc != SQLITE_OK) {
      sqlite3_close(db);
      sqlite3_close(db1);
      throw runtime_error(string("Can't open database: ") + string(sqlite3_errmsg(db)));
    }

    sqlite3_backup* bak = sqlite3_backup_init(db, "temp", db1, "main");
    if (bak == NULL) {
      sqlite3_close(db);
      sqlite3_close(db1);
      throw runtime_error(string("Can't init backup: ") + string(sqlite3_errmsg(db)));
    }

    rc = sqlite3_backup_step(bak, -1);
    if (rc != SQLITE_DONE) {
      sqlite3_close(db);
      sqlite3_close(db1);
      throw runtime_error(string("Can't step backup: ") + string(sqlite3_errmsg(db)));
    }

    rc = sqlite3_backup_finish(bak);
    if (rc != SQLITE_OK) {
      sqlite3_close(db);
      sqlite3_close(db1);
      throw runtime_error(string("Can't finish backup: ") + string(sqlite3_errmsg(db)));
    }

    sqlite3_close(db1);

    rc = sqlite3_prepare_v2(db, "select pos, ref, alt from snps where chrom=? and pos between ? and ?", -1, &intervalSearch, NULL);
    if (rc != SQLITE_OK) {
      sqlite3_close(db);
      throw runtime_error(string("Can't prepare statement: ") + string(sqlite3_errmsg(db)));
    }
  }
}
