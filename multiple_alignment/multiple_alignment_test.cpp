//-------------------------------------------------------------------------------
// 
// MultipleAlignment - Class for constructing and managing a multiple alignment
// constructed from a set of pairwise overlaps
//
// Copyright (C) 2011 Jared Simpson (jared.simpson@gmail.com)
// 
// Permission is hereby granted, free of charge, to any person obtaining a copy of
// this software and associated documentation files (the "Software"), to deal in
// the Software without restriction, including without limitation the rights to
// use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
// of the Software, and to permit persons to whom the Software is furnished to do
// so, subject to the following conditions:
// 
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
// ------------------------------------------------------------------------------
#include "multiple_alignment.h"
#include "../overlap/overlapper.h"
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <cstdlib>

int main(int argc, char** argv)
{
    std::cout << "Testing set1\n";
    std::string name0 = "read0";
    std::string name1 = "read1";
    std::string name2 = "read2";
    std::string name3 = "read3";
    std::string name4 = "read4";
    std::string name5 = "read5";
    std::string name6 = "read6";
    std::string name7 = "read7";
    std::string name8 = "read8";

    std::string sequence0 = "ATAGAGATACAACGTAAGACTACGAT";
    std::string sequence1 =  "ATACAACTAATGACTACGAT";
    std::string sequence2 =   "TACAACTAATGACTACGAT";
    std::string sequence3 =    "ACAACTAATGACCTACGAT";
    std::string sequence4 =     "CAACTAATGACCCTACGAT";
    std::string sequence5 =     "CAACTAATGACCCCTACGAT";
    std::string sequence6 =  "ATAGAGAATACAACTAATGACTACGAT";

    SequenceOverlap overlap01 = Overlapper::computeOverlap(sequence0, sequence1);
    SequenceOverlap overlap02 = Overlapper::computeOverlap(sequence0, sequence2);
    SequenceOverlap overlap03 = Overlapper::computeOverlap(sequence0, sequence3);
    SequenceOverlap overlap04 = Overlapper::computeOverlap(sequence0, sequence4);
    SequenceOverlap overlap05 = Overlapper::computeOverlap(sequence0, sequence5);
    SequenceOverlap overlap06 = Overlapper::computeOverlap(sequence0, sequence6);

    MultipleAlignment ma;
    ma.addBaseSequence(name0, sequence0);
    ma.addSequenceClipped(name1, sequence1, name0, overlap01);
    ma.addSequenceClipped(name2, sequence2, name0, overlap02);
    ma.addSequenceClipped(name3, sequence3, name0, overlap03);
    ma.addSequenceClipped(name4, sequence4, name0, overlap04);
    ma.addSequenceClipped(name5, sequence5, name0, overlap05);
    ma.addSequenceClipped(name6, sequence6, name0, overlap06);

    std::cout << "Final MA\n";
    ma.print();

    std::cout << "Testing multiple alignment from file\n";
    if(argc < 2) {
        std::cerr << "No filename given\n";
        exit(EXIT_FAILURE);
    }

    std::vector<std::string> sequences;
    
    std::string filename = argv[1];
    std::cout << "Reading sequences from " << filename << "\n";
    std::ifstream input_stream(filename.c_str());
    std::string tmp;
    while(getline(input_stream, tmp))
        sequences.push_back(tmp);

    std::cout << "Read " << sequences.size() << " sequences from file\n";

    std::cout << "Aligning all sequences against sequence[0]\n";

    std::vector<std::string> sequence_names;
    std::vector<SequenceOverlap> sequence_overlaps;
    sequence_names.push_back("root");

    MultipleAlignment real_ma;
    real_ma.addBaseSequence(sequence_names[0], sequences[0]);

    for(size_t i = 1; i < sequences.size(); ++i) {
        std::stringstream namer;
        namer << "seq-" << i;
        sequence_names.push_back(namer.str());
        SequenceOverlap overlap = Overlapper::computeOverlap(sequences[0], sequences[i]);

        real_ma.addSequenceClipped(sequence_names[i], sequences[i], sequence_names[0], overlap);
    }

    real_ma.print();
}
