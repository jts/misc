//-------------------------------------------------------------------------------
// 
// overlapper - String-string overlap algorithm 
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
    friend std::ostream& operator<<(std::ostream& out, const SequenceOverlap& overlap)
    {
        out << "[" << overlap.start_1 << " " << overlap.end_1 << "] ";
        out << "[" << overlap.start_2 << " " << overlap.end_2 << "] ";
        out << "C:" << overlap.cigar;
        return out;
    }

    // Data

    // The coordinates of the matching portion of each string
    // The end coordinate are the index of the last base matched
    int start_1;
    int end_1;
    int start_2;
    int end_2;

    //
    int score;
    int edit_distance;
    int total_columns;
    std::string cigar;

};

//
namespace Overlapper
{

// Compute the highest-scoring overlap between s1 and s2.
SequenceOverlap computeOverlap(const std::string& s1, const std::string s2);

}

#endif
