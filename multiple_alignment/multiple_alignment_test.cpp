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
#include <cstdlib>

int main(int argc, char** argv)
{
    std::cout << "Testing multiple alignment\n";
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
    for(size_t i = 1; i < sequences.size(); ++i) {
        printf("Alignment %d/%zu\n", 0, i);
        SequenceOverlap overlap = Overlapper::computeOverlap(sequences[0], sequences[i]);
        overlap.printAlignment(sequences[0], sequences[i]);
    }
}
