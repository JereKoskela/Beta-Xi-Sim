#include "ancestry.hh"
#include <cstdlib>
#include <iostream>

int main(int argc, char **argv) {
  if (argc != 3) {
    std::cout << "Call " << argv[0] << " <path to config file> ";
    std::cout << "<number of replicates>" << std::endl;
    return 1;
  }
  Ancestry ancestry(argv[1]);
  for (int i = 0; i < atoi(argv[2]); i++) {
    ancestry.simulate();
    ancestry.print_normalised_sfs_by_island();
    // ancestry.print_normalised_sfs();
    // ancestry.print_normalised_sfs_with_number_of_sites();
    // ancestry.print_fixed_s(50);
    // ancestry.print_normalised_branch_lengths();
    ancestry.reset();
  }
  return 1;
}
