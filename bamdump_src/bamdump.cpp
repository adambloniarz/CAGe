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
#include <tclap/CmdLine.h>
#include "api/BamReader.h"

using namespace std;
using namespace BamTools;

string get_cigar(vector<CigarOp> cig_ops) {

  string cigar;

  for (auto cig_op : cig_ops) {
    cigar += to_string(cig_op.Length);
    cigar += cig_op.Type;
  }
  return cigar;
}

int main(int argc, char** argv) {

  /* Command line arguments:
   *  1) Bamfile
   *  2) Contig
   *  3) Start
   *  4) End
   */

  TCLAP::CmdLine cmd("Bamdump - dump the contents of a bam file to stdout, to be piped to CAGe");

  TCLAP::UnlabeledValueArg<string> bamfile_arg("bamfile", "bam file", true, "", "bamfile", cmd);
  TCLAP::UnlabeledValueArg<string> contig_arg("contig", "contig name", true, "", "contig", cmd);
  TCLAP::UnlabeledValueArg<int> start_arg("start", "start position", true, 0, "start", cmd);
  TCLAP::UnlabeledValueArg<int> end_arg("end", "end position", true, 0, "end", cmd);
  cmd.parse(argc, argv);

  BamReader reader;
  string bamfile = bamfile_arg.getValue();
  reader.Open(bamfile);
  reader.OpenIndex(bamfile + ".bai");

  string contig = contig_arg.getValue();
  int refID = reader.GetReferenceID(contig);

  int start = start_arg.getValue();
  int end = end_arg.getValue();

  reader.SetRegion(refID, start, refID, end);

  BamAlignment al;

  while (reader.GetNextAlignment(al)) {
    printf("%i %s %s %s %i %i\n", al.Position, al.QueryBases.c_str(), get_cigar(al.CigarData).c_str(), al.Qualities.c_str(), al.MapQuality, al.IsReverseStrand() ? 1 : 0);
  }

  reader.Close();
  return 0;
}
