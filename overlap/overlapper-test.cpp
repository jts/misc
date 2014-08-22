//-------------------------------------------------------------------------------
// 
// overlapper - Basic string-string overlap algorithms 
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
#include "overlapper.h"
#include <iostream>
#include <string>

int main(int argc, char** argv)
{
    //
    std::cout << "Full length match test\n";
    std::string s1 = "AGATGATACGTAGACACGACTG";
    std::string s2 = "AGATGATACGTACGACACGACTG";
    SequenceOverlap full_12 = Overlapper::computeOverlap(s1, s2);
    std::cout << full_12 << "\n\n";
    full_12.printAlignment(s1, s2);

    SequenceOverlap full_21 = Overlapper::computeOverlap(s2, s1);
    std::cout << full_21 << "\n";
    full_21.printAlignment(s2, s1);
    
    //
    std::cout << "\nOverlap test\n";
    std::string s3 = "AGATGATACGTAGACACGACTG";
    std::string s4 = "GATGATACGTAGACACGACTGA";
    SequenceOverlap overlap_34 = Overlapper::computeOverlap(s3, s4);
    std::cout << overlap_34 << "\n\n";
    overlap_34.printAlignment(s3, s4);

    SequenceOverlap overlap_43 = Overlapper::computeOverlap(s4, s3);
    std::cout << overlap_43 << "\n";
    overlap_43.printAlignment(s4, s3);

    //
    std::cout << "\nOverlap test with indel\n";
    std::string s5 = "AGATGATACGTAGAACGACTG";
    std::string s6 = "GATGATACGGACACGACTGA";
    SequenceOverlap overlap_56 = Overlapper::computeOverlap(s5, s6);
    std::cout << overlap_56 << "\n\n";
    overlap_56.printAlignment(s5, s6);

    SequenceOverlap overlap_65 = Overlapper::computeOverlap(s6, s5);
    std::cout << overlap_65 << "\n";
    overlap_65.printAlignment(s6, s5);

    std::cout << "\nContainment test\n";
    std::string s7 = "AGATGATACGGTAGAACGACTG";
    std::string s8 = "ATGATACGTAGAACGACT";
    SequenceOverlap contain_78 = Overlapper::computeOverlap(s7, s8);
    std::cout << contain_78 << "\n\n";
    contain_78.printAlignment(s7, s8);

    SequenceOverlap contain_87 = Overlapper::computeOverlap(s8, s7);
    std::cout << contain_87 << "\n";
    contain_87.printAlignment(s8, s7);
    
    std::cout << "\nReal read test\n";
    std::string a, b;
    a = "GGGTCGTTACGCTGCAGTCGCCGCCTTCATCTTCAACACCACCGAATTTGTCCCTGTTGGCCTGCTCTC"
        "TGACATTGCGCAAAGTTTTCACATGCAAACCGCTCAGGTCGGCATCATGTTGACCATTTACGCATGGGT"
        "AGTAGCGCTAATGTCATTGCCTTTTATGTTAATGACCAGTCAGGTTGAACGGCGCAAATTACTGATCTG"
        "CCTGTTTGTGGTGTTTATTGCCA";

    b = "GCCTTCATCTTCAACACCACCGAATTTGTCCCTGTTGGCCTGCTCTCTGACATTGCGCAAAGTTTTCAC"
        "ATGCAAACCGCTCAGGTCGGCATCATGTTGACCATTTACGCATGGGTAGTAGCGCTAATGTCATTGCCT"
        "TTTATGTTAATGACCAGTCAGGTTGAACGGCGCAAATTACTGATCTGCCTGTTTGTGGTGTTTATTGCC"
        "AGCCACGTACTGTCGTTTTTGTCGTGGAGCTTTACCGTTCTGGTGATA";

    SequenceOverlap ovr_ab = Overlapper::computeOverlap(a, b);
    std::cout << ovr_ab << "\n\n";
    ovr_ab.printAlignment(a, b);

    a = "TTCTACTAATTGTATGTGCCACAGTATATTTATCGTTAGTTTTTAAAACCACAATACTAAATACCTTGTTAAGTCCTGATCAATTTGGTCTGGACATCT";
    b = "TCTACTAATTGTATGTGCCACAGTATATTTATCGTTAGTTTTTAAAAACCACAATACTAAATACCTTGTTAAGTCCTGATCAATTTGGTCTGACATCCA";

    ovr_ab = Overlapper::computeOverlap(a, b);
    std::cout << ovr_ab << "\n\n";
    ovr_ab.printAlignment(a, b);

    a = "CAATGCCCCCCCCCCCGCAGGGGGGGAGGGGGGGGGGCTGCTTTTTTTTTTAATAAAAAAAAAATTTTT"
        "TTTTTACGTCTGTTTTTTTCACGGGGGGGGGTTTTTATTTTTCTCCCCCCCCGTGAGTCAGCGTTATCT"
        "TTTTTCTTTTTTGCTGACCGACCCCCCAGTTCACGCACCCCCGCCCATCACGCTGACACCCCGGCAACC"
        "CCAACCCCAACCCCAACACCCCGCCCCGAGCCAAACCCAATTAGCGCCCCTGCAGAGACCGTTGCCCCG"
        "CCCCTCCTTACAGTTAACATCCCGGTTTTCTCC";
    b = "GACTGACAGACGCAGCGCGCCCGTCAGCAATGCCCCGCAGGAGGCTGCTTTAATAAAATTACGTCTGTT"
        "CACGGTTATTCTCCCCGTGAGTCAGCGTTATCTTTCTTTTGCTGACGACCCAGTTCACGCACCGCCATC"
        "ACGCTGACACCGCAACCAACCCAACCACAACCGCCGAGCCAACCAATCAGGCAGCCCTAGACTCGCCAT"
        "CGAGAGCCGC";

    ovr_ab = Overlapper::extendMatch(a, b, 117, 81, 20);
    ovr_ab.printAlignment(a, b);

    std::cout << "Testing mismatch at end\n";
    a = "TAGTAGAGATGGGGTTTCACCATGTTGGCTAGGATGGTCTCCATCTCCTGACCTTGTGATCCTCCCGCCTCGGCCTCCCA";
    b = "TGCTTCAACTTCCCGAGTAGCTGGGATTACAGGTGTCTGCCACCACACCTGGTTAATTTTTGTATTTTCAGTAGAGATGGGGTTTCACCATGTTGGCTAGG";
    
    ovr_ab = Overlapper::computeOverlap(a, b);
    std::cout << ovr_ab << "\n\n";
    ovr_ab.printAlignment(a, b);

    std::cout << "Testing affine alignment\n";
    a = "AAGAGGTATAGTCTGTCATCCTAGCTCCTACTACATAACCACATTTAGAG"
        "CACATGGTTAATGCGTACCTTAAATGTATAAATAGATTGTTAGACATATA"
        "TTTTGATAAATACTTCAGAGACAAATAGTATTATGACAAATCAATGCTAC"
        "ATTGGAGGGGAGGGGAATTAAATGGGTTCCATGTGTATTTTTCTTCCACT"
        "GGCCCATAGATAATTATTATTTCTGGAATAAAAGATTAGGGATAGATCAC"
        "TGAATGCTTCAGAACCTTATCTTAATTCTATTTGATTTTTGGTGATATAA"
        "AGCCAGCCCAAGAGTAATCTTCTGAGGGTAGCTATCTTTATAACTTGTTT"
        "GCTTCTTGGCAAACAGAAATATAGAGCTTGCATCATTTGAGAACAAAAGA"
        "GGGAAACCCAATGATTAATAATTCAGAAATCTCTCCTGAGGTCTTCTCAG"
        "CTATATTTCTGGACTTTGCCTGAATACTATTAACTGAACATTTCTCTTAC"
        "TGAGCTAATTGTCTCTTGCTCCTTCTCCTCCTCCTATTTTGGCCAACCAT"
        "GTACTTTAGCAGAAACCTAGGAATTCATCCTGCATGGTGACTGCCACTGC"
        "CCTAGTTCAGGCCATCACTGTTTTTGGTCTAGAGGGTGTCAACAGCTGCC"
        "TGACTCACACCGCCTTC";
    
    b = "TTTTGATAAATACTTCAGAGACAAATAGTATTATGACAAATCAATGCTAC"
        "ATTGGAGGGGAGGGGAATTAAATGGGTTCCATGTGTATTTTTCTTCCACT"
        "T"
        "TCTTCTCAGCTATATTTCTGGACTTTGCCTG"
        "AATACTATTAACTGAACATTTCTCTTACTGAGCTAATTGTCTCTTGCTCC"
        "TTCTCCTCCTCCTATTTTGGCCAACCATGTACTTTAGCAGAAACC";

    OverlapperParams params = affine_default_params;
    
    params.gap_penalty = -100;
    params.gap_ext_penalty = 0;
    params.type = ALT_CONTAINMENT;
    
    ovr_ab = Overlapper::computeAlignmentAffine(a, b, params);
    std::cout << ovr_ab << "\n\n";
    ovr_ab.printAlignment(a, b);

    a = "AGCACTACTTTATAGGGCCAACCCTATTGTTTTTAATGACCTGACACTATGATTTGGGGT"
        "TCTGGGGTGTTAAAGAAGCACTCTGGCTACTGGTACAGAGGATGTCATTAAGAAGGAGGG"
        "CAGCCATTGGCCCACCTTCCACCAATGGGAAGCTTTAGCTAGAGGGGACATTACACCTCA"
        "CAATAGAGACATAATGGCTCCTCCCTTCCTTGGGTCCCCTCAGCCTCTGAGAGAAACTTC"
        "TGCTGAAATGTAAATAATAATGTATTAGTAATACTACTAGTAATAGAAGTGATGACTTCT"
        "AGAGCCCTCTCAAGTTCTTCTCCTAAGGACATTTCAGCACAAGTTTCTCTCAGAGGCTGA"
        "GGGGACCCAAGGGAGGGAAAAACCATAATGTCTCCATTGTGAGCCATAATATCCTCTCTA"
        "GCTAAAGCTTCCCATTGGTGGAAGGTGGGCCAATGGCTGCCGTCCTTCTTAATGACATCC"
        "TCTGTAACAGCACAAGGTGCTTAACACCTCAGATTCCGAAATAAAGGGTTGAACCATTAA"
        "AAGAAACATGGGTGGCCCTAGAAATTAGTGTTCTCTCCCTCCCA";

    b = "CTACTTTATAGGGCCAACCCTATTGTTTTTAATGACCTGACACTATGATTTGGGGT"
        "TCTGGGGTGTTAAAGAAGCACTCTGGCTACTGGTACAGAGGATGTCATTAAGAAGGTGGG"
        "CCAATGGCTGCCGTCCTTCTTAATGACATCCTCTGTAACAGCACAAGGTGCTTAACACCT"
        "CAGATTCCGAAATAAAGGGTTGAACCATTAAAAGAAACATGGGTGGCCCTAGAAATTAGT"
        "GTTCTCTCCCT";
    
    params.gap_ext_penalty = 0;
    params.type = ALT_CONTAINMENT;
    ovr_ab = Overlapper::computeAlignmentAffine(a, b, params);
    std::cout << ovr_ab << "\n\n";
    ovr_ab.printAlignment(a, b);
}
