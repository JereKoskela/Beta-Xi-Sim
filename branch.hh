#ifndef BRANCH
#define BRANCH

#include <cstdlib>
#include <iostream>
#include <vector>

struct Branch {

  Branch()
      : locus(), island(), selective_type(), virtual_flag(), leaf_epoch(),
        leaf_time(), parents(), children(), incoming(), checking(),
        ancestral_blocks(), order() {}

  Branch(const int locus_, const int island_, const int epoch,
         const double left, const double right, double leaf_time_)
      : locus(locus_), island(island_), selective_type(-1), virtual_flag(0),
        leaf_epoch(epoch), leaf_time(leaf_time_), parents(), children(),
        incoming(), checking(), ancestral_blocks(2, left), order() {
    ancestral_blocks[1] = right;
  }

  int contains(const double x) const {
    // returns 1 if point x lies inside ancestral material (denoted by an
    // interval (ancestral_blocks[k], ancestral_blocks[k + 1]) for some
    // even k. returns 0 otherwise.
    int ret = 0;
    if (x > ancestral_blocks[0] && x < ancestral_blocks.back()) {
      int k = 0;
      while (x > ancestral_blocks[k]) {
        k++;
      }
      ret = k % 2;
    }
    return ret;
  }

  int locus, island, selective_type, virtual_flag, leaf_epoch;
  double leaf_time;
  std::vector<int> parents, children, incoming, checking;
  std::vector<double> ancestral_blocks;
  std::vector<std::vector<int>> order;
};

#endif
