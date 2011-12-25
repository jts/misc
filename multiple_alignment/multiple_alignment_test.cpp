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
#include "multiple_alignment_builder.h"
#include "overlapper.h"
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <cstdlib>

int main(int argc, char** argv)
{
    std::cout << "Testing set1\n";
    std::vector<std::string> test_sequences_0;

    test_sequences_0.push_back("ATAGAGATACAACGTAAGACTACGAT");
    test_sequences_0.push_back("ATACAACTAATGACTACGAT");
    test_sequences_0.push_back("TACAACTAATGACTACGAT");
    test_sequences_0.push_back("ACAACTAATGACCTACGAT");
    test_sequences_0.push_back("CAACTAATGACCCTACGAT");
    test_sequences_0.push_back("CAACTAATGACCCCTACGAT");
    test_sequences_0.push_back("ATAGAGAATACAACTAATGACTACGAT");

    MultipleAlignment ma_0 = MultipleAlignmentBuilder::alignToFirst(test_sequences_0);
    ma_0.print();
    ma_0.printPileup();

    std::cout << "\nTesting multiple alignment from file\n";
    if(argc < 2) {
        std::cerr << "No filename given\n";
        exit(EXIT_FAILURE);
    }

    std::vector<std::string> test_sequences_1;
    
    std::string filename = argv[1];
    std::cout << "Reading sequences from " << filename << "\n";
    std::ifstream input_stream(filename.c_str());
    std::string tmp;
    while(getline(input_stream, tmp))
        test_sequences_1.push_back(tmp);

    std::cout << "Read " << test_sequences_1.size() << " sequences from file\n";

    MultipleAlignment ma_1 = MultipleAlignmentBuilder::alignToFirst(test_sequences_1);
    ma_1.print();
    ma_1.printPileup();
}
