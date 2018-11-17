#include <cstdlib>
#include <iostream>
#include <vector>
#include "ancestry.hh"

int main(int argc, char **argv) {
    if (argc != 3) {
        std::cout << "Call " << argv[0] << " <path to config file> ";
        std::cout << "<number of replicates>" << std::endl;
        return 1;
    }
    Ancestry ancestry(argv[1]);
    for (int i = 0; i < atoi(argv[2]); i++) {
        ancestry.simulate();
        ancestry.print_normalised_sfs();
        //ancestry.print_fixed_s(50);
        //print_normalised_branch_lengths();
        ancestry.reset(argv[1]);
    }
    return 1;
}
