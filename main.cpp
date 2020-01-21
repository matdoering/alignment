#include "Alignment.h"

#include <string>
#include <iostream>

int main(int argc, char** argv) {

    // RNA-RNA
    std::string s1 = "UACGAUGAGAUU";
    std::string s2 = "UAAAAACGAUGAGAAU";

    Align::Alignment ali(s1,s2);
    ali.align();
    std::cout << ali << "\n" << std::endl;

    // DNA-DNA
    s1 = "TACGAGGATA";
    s2 = "TACGATGATA";
    ali = Align::Alignment(s1,s2);
    ali.align();
    std::cout << ali << std::endl;
}
