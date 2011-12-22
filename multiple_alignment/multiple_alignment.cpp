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
#include <iostream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <assert.h>
#include <stdio.h>

//#define MA_DEBUG 1

//
// MultipleAlignmentElement
//

MultipleAlignmentElement::MultipleAlignmentElement(const std::string& _name, 
                                                   const std::string& _sequence,
                                                   size_t offset) : name(_name), 
                                                                    sequence(_sequence),
                                                                    column_offset(offset)
{
}

//
int MultipleAlignmentElement::getPaddedPositionOfBase(size_t idx) const
{
    size_t unpadded_count = 0;
    for(size_t i = 0; i < sequence.size(); ++i) {
        if(sequence[i] != '-') {
            if(unpadded_count == idx)
                return i;
            else
                unpadded_count += 1;
        }
    }
    std::cerr << "Base index out of bounds: " << idx << "\n";
    assert(false);
    return -1;
}

//
void MultipleAlignmentElement::insertGapBeforeColumn(size_t column_index)
{
    // Check if the column to insert the gap falls within the offset region
    // If so, just increase the offset to account for the shift in columns
    // If the column index is one greater than the offset, then we want to 
    // insert a gap before the first base. This is equivalent to just
    // extending the offset.
    if(column_index <= column_offset + 1) {
        column_offset += 1;
    }
    else {
        assert(column_index > column_offset);
        sequence.insert(column_index - column_offset, 1, '-');
    }
}

//
std::string MultipleAlignmentElement::getUnpaddedSequence() const
{
    std::string out;
    for(size_t i = 0; i < sequence.size(); ++i) {
        if(sequence[i] != '-')
            out.push_back(sequence[i]);
    }
    return out;
}

//
// MultipleAlignment
//

//
void MultipleAlignment::addBaseElement(const std::string& name, const std::string& sequence)
{
    m_sequences.push_back(MultipleAlignmentElement(name, sequence, 0));
}

// See header
void MultipleAlignment::addSequenceClipped(const std::string& incoming_name,
                                           const std::string& incoming_sequence,
                                           const std::string& reference_name,
                                           const SequenceOverlap& reference_incoming_overlap)
{
    //
    MultipleAlignmentElement* template_element = NULL;
    for(size_t i = 0; i < m_sequences.size(); ++i) {
        if(m_sequences[i].name == reference_name)
            template_element = &m_sequences[i];
    }

    _addSequence(incoming_name, incoming_sequence, template_element, reference_incoming_overlap);
}

//
void MultipleAlignment::addSequenceExtendLast(const std::string& incoming_name,
                                              const std::string& incoming_sequence,
                                              const std::string& reference_name,
                                              const SequenceOverlap& previous_incoming_overlap)
{
    // This function cannot be called before a base element has been added
    assert(!m_sequences.empty());

    // Get an iterator to the sequence already existing in the multiple alignment
    // that we use to detemine how to pad the new sequence.
    MultipleAlignmentElement* template_element = &m_sequences.back();
    _addSequence(incoming_name, incoming_sequence, template_element, previous_incoming_overlap);
}

// Adds a new string into the multiple alignment using the 
// padded sequence of an existing sequence to set up the new padded
// string
void MultipleAlignment::_addSequence(const std::string& name, 
                                     const std::string& sequence, 
                                     MultipleAlignmentElement* template_element, 
                                     const SequenceOverlap& overlap)
{
    // Ensure the name of the incoming element is unique
    for(size_t i = 0; i < m_sequences.size(); ++i) {
        if(m_sequences[i].name == name) {
            std::cerr << "Error in constructing multiple alignment, non-unique name " << name << "\n";
            exit(EXIT_FAILURE);
        }
    }

    // Get the padded sequence for the template element
    const std::string& template_padded = template_element->sequence;

    // The output padded sequence for the incoming
    std::string padded_output;

    // Expand the cigar for easier parsing
    std::string expanded_cigar = expandCigar(overlap.cigar);
    assert(!expanded_cigar.empty());

#ifdef MA_DEBUG
    std::cout << "Cigar: " << expanded_cigar << "\n";
    std::cout << "template: " << template_padded << "\n";
    std::cout << "incoming: " << sequence << "\n";
#endif

    // Iterate over the cigar string and the padded sequence of the template element
    // to determine where to insert gap symbols
    size_t cigar_index = 0;
    size_t template_index = template_element->getPaddedPositionOfBase(overlap.start_1);
    size_t incoming_index = 0;
    size_t template_offset = template_element->column_offset;
    size_t incoming_offset = template_index + template_offset;

#ifdef MA_DEBUG
    std::cout << "Pairwise:\n";
    overlap.printAlignment(template_element->getUnpaddedSequence(), sequence);
#endif

    assert(template_index < template_padded.size());
    assert(template_padded[template_index] != '-');

    while(cigar_index < expanded_cigar.size()) {

        // Check if we are in an existing template gap. This must be handled
        // seperately
        bool in_template_gap = template_padded[template_index] == '-';
        if(in_template_gap) {
            // If we are in a incoming sequence insertion
            // (cigar D) then we are adding a base into a known
            // gap. Add the current incoming base to the output
            if(expanded_cigar[cigar_index] == 'D') {
                padded_output.push_back(sequence[incoming_index++]);
                cigar_index += 1;
                template_index += 1;
            } else { 
                // This is an insertion that is in some other sequence
                // in the multiple alignment. Add a gap to the padded output
                padded_output.push_back('-');
                
                // Increment the template index
                template_index += 1;    
            }
        } else {
            // Not a template gap
            switch(expanded_cigar[cigar_index]) {
                case 'M':
                    padded_output.push_back(sequence[incoming_index]);
                    incoming_index += 1;
                    template_index += 1;
                    cigar_index += 1;
                    break;
                case 'D':
                    insertGapBeforeColumn(template_index + template_offset);
                    padded_output.push_back(sequence[incoming_index]);
                    incoming_index += 1;
                    cigar_index += 1;
                    template_index += 1; // skip the newly introduced gap
                    break;
                case 'I':
                    padded_output.push_back('-');
                    cigar_index += 1;
                    template_index += 1;
                    break;
            }
        }
    }

    m_sequences.push_back(MultipleAlignmentElement(name, padded_output, incoming_offset));
}

void MultipleAlignment::print() const
{
    for(size_t i = 0; i < m_sequences.size(); ++i) {
        std::string padding = std::string(m_sequences[i].column_offset, ' ');
        printf("\t%s%s\t%s\n", padding.c_str(),
                               m_sequences[i].sequence.c_str(), 
                               m_sequences[i].name.c_str());
    }
}

//
void MultipleAlignment::insertGapBeforeColumn(size_t column_index)
{
    for(size_t i = 0; i < m_sequences.size(); ++i) {
        m_sequences[i].insertGapBeforeColumn(column_index);
    }
}

//
std::string MultipleAlignment::expandCigar(const std::string& cigar)
{
    std::string out;
    std::stringstream parser(cigar);
    int length;
    char symbol;
    while(parser >> length >> symbol)
        out.append(length, symbol);
    return out;
}