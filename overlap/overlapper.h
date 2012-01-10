//-------------------------------------------------------------------------------
// 
// overlapper - Functions to calculate overlaps between pairs of strings
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
#ifndef OVERLAPPER_H
#define OVERLAPPER_H

#include <string>
#include <ostream>
#include <assert.h>

// A start/end coordinate pair representing
// a subsequence. The end coordinate is
// the index of the last base aligned.
struct SequenceInterval
{
    // functions
    
    // Check that the interval is valid
    bool isValid() const { return start <= end; }
    
    // Change the interval to represent the same
    // set of bases but on the opposite strand.
    void flipStrand(int sequence_length)
    {
        assert(isValid());
        int tmp = sequence_length - start - 1;
        start = sequence_length - end - 1;
        end = tmp;
        assert(isValid());
    }

    // data
    int start;
    int end; // inclusive
};

// Data structure to hold the result of
// an overlap calculation
struct SequenceOverlap
{
    // Functions
    
    // Check that the record is properly formed
    bool isValid() const;

    // Print the alignment with padding characters
    void printAlignment(const std::string& s1, const std::string& s2) const;
    
    // Return the percent identity which we define to be
    // the number of matching columns divided by the total number of columns
    double getPercentIdentity() const;

    // Returns the length of the overlap, defined to be the 
    // number of columns in the alignment
    int getOverlapLength() const { return total_columns; }

    //
    friend std::ostream& operator<<(std::ostream& out, const SequenceOverlap& overlap);

    // Data

    // The coordinates of the matching portion of each string
    // The end coordinate are the index of the last base matched
    SequenceInterval match[2];
    
    // The length of the input sequences
    int length[2];

    //
    int score;
    int edit_distance;
    int total_columns;

    // The cigar string follows the sam convention with s1 being the "reference":
    // I is an insertion into s1
    // D is a deletion from s1

    // A-C s1
    // AAC s2
    // C: 1M1I1M
    //
    // ATC s1
    // A-C s2
    // C: 1M1D1M
    std::string cigar;

};

//
namespace Overlapper
{

// Compute the highest-scoring overlap between s1 and s2.
// This is a naive O(M*N) algorithm
SequenceOverlap computeOverlap(const std::string& s1, const std::string& s2);

// Extend a match between s1 and s2 into a full overlap using banded dynamic programming
SequenceOverlap extendMatch(const std::string& s1, const std::string& s2, SequenceOverlap seed_match, int bandwidth);

}

#endif
