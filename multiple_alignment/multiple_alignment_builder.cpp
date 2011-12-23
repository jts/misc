//-------------------------------------------------------------------------------
// 
// MultipleAlignmentBuilder - Utility functions to generate multiple alignments
// from a set of sequences
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
#include "multiple_alignment_builder.h"
#include "overlapper.h"
#include <sstream>
#include <assert.h>
#include <iostream>

namespace MultipleAlignmentBuilder
{

//
MultipleAlignment alignToFirst(const std::vector<std::string>& sequences)
{
    assert(sequences.size() > 0);

    std::string base_name = "sequence-0";

    MultipleAlignment out_alignment;
    out_alignment.addBaseSequence(base_name, sequences[0]);

    for(size_t i = 1; i < sequences.size(); ++i) {
        SequenceOverlap overlap = Overlapper::computeOverlap(sequences[0], sequences[i]);
        std::stringstream namer;
        namer << "sequence-" << i;
        out_alignment.addSequenceClipped(namer.str(), sequences[i], base_name, overlap);
    }

    return out_alignment;
}

//
MultipleAlignment greedyExtended(const std::vector<std::string>& sequences)
{
    (void)sequences;
    MultipleAlignment out_alignment;
    return out_alignment;
}

}
