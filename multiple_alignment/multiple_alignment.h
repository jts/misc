//-------------------------------------------------------------------------------
// 
// MultipleAlignment - Class for progressively constructing and managing
// a multiple alignment from a set of pairwise overlaps
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
#ifndef MULTIPLE_ALIGNMENT_H_
#define MULTIPLE_ALIGNMENT_H_

#include "overlapper.h"
#include <vector>

struct MultipleAlignmentElement
{
    // Functions
    MultipleAlignmentElement(const std::string& _name, const std::string& _sequence, size_t offset);

    // Returns the sequence with all padding characters removed
    std::string getUnpaddedSequence() const;

    // Returns the position in the padded string of the base at index idx of
    // the unpadded sequence.
    // Precondition: idx is less than the number of sequence bases
    int getPaddedPositionOfBase(size_t idx) const;

    // Insert a new gap before the specified column
    void insertGapBeforeColumn(size_t column_index);
    
    // Data
    std::string name;
    std::string padded_sequence;
    std::string padded_quality;

    // The number of columns in the multiple alignment before/after
    // the sequence data for this sequence starts
    size_t leading_columns;
    size_t trailing_columns;

};

class MultipleAlignment
{
    public:
        MultipleAlignment() {}

        // Add the first element to the multiple alignment
        void addBaseSequence(const std::string& name, const std::string& sequence);

        // Add a new sequence to the multiple alignment, which overlaps
        // the first sequence added to the multiple alignment (the base sequence).
        // Preconditions are that the base sequence exists and the overlap
        // that is passed in refers to the overlap between the base sequence (as
        // the first set of coordinates) and the incoming sequence. 
        void addOverlap(const std::string& incoming_name,
                        const std::string& incoming_sequence,
                        const SequenceOverlap& reference_incoming_overlap);

        // Add a new sequence to the multiple alignment by extending the last sequence added.
        // This function allows the progressive construction of a multiple alignment for a series
        // of sequences that are laid out into a contig. The SequenceOverlap object must be
        // defined such that the coordinates for the sequence already in the multiple alignment
        // are defined first, and the coordinates for the incoming sequence are defined second.
        void addExtension(const std::string& incoming_name,
                          const std::string& incoming_sequence,
                          const SequenceOverlap& previous_incoming_overlap);

        // Print the alignment to stdout
        void print() const;
    
    private:
     
        // Internal function for performing the addition of a new sequence. Called by
        // addSequenceClipped/addSequenceExtend
        void _addSequence(const std::string& name, 
                          const std::string& sequence, 
                          MultipleAlignmentElement* template_element, 
                          const SequenceOverlap& overlap);
        
        // Insert a new gap into all sequences in the multiple alignment
        // after the given column
        void insertGapBeforeColumn(size_t column_index);

        // Expand a cigar string by having one symbol per event instead
        // of run length encoding
        std::string expandCigar(const std::string& cigar);

        // Data
        std::vector<MultipleAlignmentElement> m_sequences;
};

#endif  // MULTIPLE_ALIGNMENT_H_
