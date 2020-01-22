#include "Alignment.h"

#include <string>
#include <iostream>
#include <array>
#include <algorithm>
#include <cassert>

void printHelp() {
    std::cout << "Usage: Alignment <seq1> <seq2>\n" 
              << "seq is either a DNA or RNA seq\n"
              << "Parameters:\n"
              << "-l Use local alignment (Smith-Waterman)"
              << std::endl;
}

void test() {
    // RNA-RNA
    std::string s1 = "UACGAUGAGAUU";
    std::string s2 = "UAAAAACGAUGAGAAU";

    Align::Alignment ali(s1,s2);
    ali.align(false);
    std::cout << ali << "\n" << std::endl;

    // DNA-DNA
    s1 = "TACGAGGATA";
    s2 = "TACGATGATA";
    ali = Align::Alignment(s1,s2);
    ali.align(false);
    std::cout << ali << std::endl;
}

int main(int argc, char** argv) {
    // handle arguments
    if (argc < 3) {
        printHelp();
        // test();
        return -1;
    }
    bool useLocalAlignment = false;
    std::array<std::string, 2> seqs;
    size_t seqIdx = 0;
    for (size_t i = 1; i < argc; ++i) {
        if (std::string(argv[i]) == "-l") {
            useLocalAlignment = true;
        } else {
            assert(seqIdx < seqs.size());
            std::string s(argv[i]);
            // only upper-case letters are allowed
            std::transform(s.begin(), s.end(), s.begin(), ::toupper);
            seqs[seqIdx++] = s;
        }
    }
    Align::Alignment ali(seqs[0], seqs[1]);
    bool localAlignment = false;
    if (useLocalAlignment) {
        localAlignment = true;
    }
    ali.align(localAlignment);
    std::cout << ali << "\n" << std::endl;
}
