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

int main( int argc, char *argv[] ) {


    /* Command line arguments:
     *  1) Bamfile
     *  2) Contig
     *  3) Start
     *  4) End
     */
    
    BamReader reader;
    string bamfile = argv[1];
    reader.Open(bamfile);
    reader.OpenIndex(bamfile + ".bai");
    int refID = reader.GetReferenceID( argv[2] );
    reader.SetRegion( refID, atoi(argv[3]), refID, atoi(argv[4]) );

    BamAlignment al;

    while ( reader.GetNextAlignment(al) ) {
        printf("%i %s %s %s %i %i\n", al.Position, al.QueryBases.c_str(), get_cigar(al.CigarData).c_str(), al.Qualities.c_str(), al.MapQuality, al.IsReverseStrand() ? 1 : 0);
    }

    reader.Close();

    return 0;
}
