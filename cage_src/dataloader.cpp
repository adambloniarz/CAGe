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
#include <math.h>
#include <cctype>
#include <utility>
#include <vector>
#include <tuple>
#include <ctype.h>
#include <algorithm>
#include <queue>
#include <iostream>
#include "dataloader.h"
#include "variantcaller.h"
#include "alignment.h"
#include "pileslist.h"
#include "snpdb.h"

using namespace std;
using cage::DataLoaderParams;
using cage::Pileslist;

namespace cage {

  Dataloader::Dataloader (DataLoaderParams params, string&& ref_) :
    params(params), ref(ref_), dataend(params.start), 
    all_stats((params.end - params.start) / params.step) {

      // Make the ref string all uppercase, for consistency
      transform(ref.begin(), ref.end(), ref.begin(), ::toupper);

      // Create the variant reporter
      if (!params.call_variants) {
        snp_db = unique_ptr<SNPdb>(new SNPdb(params.input_SNP_db));
      }
    }

  // This will be run on a separate thread
  void Dataloader::loadfunc() {
    queue<string> next_overlap;
    queue<string> prev_overlap;
    string al_line;

    for (int block_start = params.start; block_start < params.end; block_start += load_step) {
      int block_end = min(block_start + load_step, params.end);

      // Allocate the piles for this section
      Pileslist piles_seg(block_end - block_start, params.contig);

      unique_ptr<VariantCaller> vc;

      if (params.call_variants) {
        // Allocate VariantCaller
        vc = unique_ptr<VariantCaller>(new VariantCaller(block_end - block_start, block_start, ref, params.contig));
      }

      while (true) {
        if (!prev_overlap.empty()) {
          al_line = prev_overlap.front();
          prev_overlap.pop();
        } else {
          getline(cin, al_line);
          if (al_line.size() == 0) {
            break;
          }
        }

        Alignment al = Alignment(al_line);

        int offset = 0;
        if (al.CigarData.size() > 0 && al.CigarData.front().Type == 'S') {
          offset = al.CigarData.front().Length;
        }

        if (al.Position >= block_end) {
          next_overlap.push(al_line);
          break;
        } 

        // This will take into account soft clipping
        int read_start = al.Position - offset;
        int pos_rel = read_start - block_start;

        if (pos_rel >= 0 && pos_rel < (int)piles_seg.rdstarts.size()) {
          piles_seg.rdstarts[pos_rel]++;
          if (al.MapQuality == 0) {
            piles_seg.zero_mapq[pos_rel]++;
          }
        }

        if (al.CigarData.size() == 0) {
          continue;
        }

        int j = 0; // Will store position in genome coordinates
        int k = 0; // Will store position in QueryBases coordinates

        for (CigarOp op : al.CigarData) { 
          switch (op.Type) {
            case ('M'):
            case ('='):
            case ('X'):
              for (int l = 0; l < op.Length; l++) {
                if (al.MapQuality != 0 && al.Qualities.at(k) > '%' && j + pos_rel >= 0 && j + pos_rel < (int)(block_end - block_start)) {
                  // Only count up mismatches; do not consider variants in soft-clipped region
                  piles_seg.depth[j + pos_rel]++;
                  if (al.QueryBases.at(k) == ref.at(j + pos_rel + block_start)) {
                    piles_seg.matches[j + pos_rel]++;
                  } else {
                    if (params.call_variants) {
                      vc->increment_mismatch(j + pos_rel, al.QueryBases.substr(k, 1), al.IsReverseStrand);
                    }
                  }
                }
                j++; k++;
              }
              break;

            case ('S'): // Soft-clipping
              for (int l = 0; l < op.Length; l++) {
                if (al.MapQuality != 0 && al.Qualities.at(k) > '%' && j + pos_rel >= 0 && j + pos_rel < (int)(block_end - block_start)) {
                  // Only increment the depth and number of matches; do not consider variants in soft-clipped region
                  piles_seg.depth[j + pos_rel]++;

                  if (al.QueryBases.at(k) == ref.at(j + pos_rel + block_start)) {
                    piles_seg.matches[j + pos_rel]++;
                  } else {
                    //  cout << "Mismatch at " << j + pos_rel + block_start << seq[j] << memblock[j + pos_rel] << endl;
                  }
                }
                j++; k++;
              }
              break;

            case ('I'):
              if (al.MapQuality != 0 && j + pos_rel - 1 >= 0 && j + pos_rel - 1 < (int)(block_end - block_start)) {
                piles_seg.indels[j + pos_rel - 1]++;
                string alt_bases;
                if (k == 0) {
                  alt_bases = ref.at(j + pos_rel + block_start - 1) + al.QueryBases.substr(k, op.Length);
                } else {
                  alt_bases = al.QueryBases.substr(k - 1, op.Length + 1);
                }
                if (params.call_variants) {
                  vc->increment_mismatch(j + pos_rel - 1, alt_bases, al.IsReverseStrand);
                }
              }
              k += op.Length;
              break;

            case ('D'):
              if (al.MapQuality != 0 && j + pos_rel - 1 >= 0 && j + pos_rel < (int)(block_end - block_start)) {
                piles_seg.indels[j + pos_rel]++;
                if (params.call_variants) {
                  vc->increment_mismatch(j + pos_rel - 1, string(op.Length, '-'), al.IsReverseStrand);
                }
              }
              for (int l = 0; l < op.Length; l++) {
                piles_seg.depth[j + pos_rel]++;
                j++;
              }
              break;

              // invalid CIGAR op-code
            default:
              const string message = string("invalid CIGAR operation type: '") + op.Type + "'";
              cout << message << endl;
              throw message;
          }
        }
        if (al.Position - offset + k >= block_end) {
          next_overlap.push(al_line);
        }
      }

      prev_overlap = move(next_overlap);
      next_overlap = queue<string>(); // Reset next_overlap

      // Get the variants
      vector<tuple<int, int, int>> snps;
      if (params.call_variants) {
        snps = vc->report_variants(piles_seg.depth, params.output_VCF);
      } else {
        snps = snp_db->report_variants(params.contig, block_start, block_end);
      }


      for (const tuple<int, int, int>& p : snps) {
        int i = get<0>(p);
        int ref_len = get<1>(p);
        int alt_len = get<2>(p);

        piles_seg.snp_indicators[ i - block_start ] = 1;

        // In case of a deletion
        if (ref_len > 1) {

          for (int k = i + 1; k < i + ref_len; k++) {
            if (k < (int)block_end) {
              piles_seg.indels[ k - block_start ] = 0;
              piles_seg.matches[ k - block_start ]++;
            }
          }

          // Make sure number of matches does not exceed depth
          if (piles_seg.matches[i - block_start] > piles_seg.depth[i - block_start])
            piles_seg.matches[i - block_start] = piles_seg.depth[i - block_start];
        }

        // In case of an insertion
        if (alt_len > 1) {
          piles_seg.indels[i - block_start] = 0;
        }

        // In case of a SNP
        if (ref_len == 1 && alt_len == 1) {
          piles_seg.matches[i - block_start] = piles_seg.depth[i - block_start];
        }
      }

      // Output progress
      if (params.verbose) {
        for (int i = block_start; i < (int)block_end; i+= 100000) {
          int idx = i - block_start;
          printf("Position: %i\tRdstarts: %i\tMatches: %i\tDepth: %i\n", i, piles_seg.rdstarts[idx], piles_seg.matches[idx], piles_seg.depth[idx]);
        }
      }

      // Insert sufficient statistics into the data store 
      for (int pos = block_start; pos < (int)block_end; pos += params.step) {
        all_stats[(pos - params.start) / params.step] = SuffStats(piles_seg, pos, pos + params.step, pos - block_start);
      }

      // Signal that loading of current block is completed
      dataend = block_end;
      load_request.notify_one();
    }
  }

  SuffStats Dataloader::get_region(const int& seg_start, const int& seg_end) {
    // The case where we haven't loaded the requested data yet
    if (seg_end > dataend) {
      unique_lock<mutex> l_lk(loading_lock);

      // Wait for the loading thread to actually load the data
      load_request.wait(l_lk, [&] { return (dataend >= seg_end); });
    }
    return SuffStats(all_stats, (seg_start - params.start) / params.step, (seg_end - params.start) / params.step);
  }
}
