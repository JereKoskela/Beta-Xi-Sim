# Beta-Xi-Sim
An efficient simulator for ancestries of multiple unlinked loci driven by Beta-Xi-coalescents with exponential population growth. Valid for any number of loci and sample sizes &lt; 10 000.

# User interface
Parameters can be specified via a text file interface (see provided dev.cfg for an example).
More powerful customisation, such as arbitrary fitness functions and recombination maps, can be done by specifying the user-defined functions at the top of the Ancestry struct in ancestry.hh.

# Output
The software outputs a space-delimited, normalised site frequency spectrum (SFS) for each locus and each iteration. Only one SFS is output per locus regardless of the number of islands - population subdivision is not reflected in the implemented default output. If loc_i_sfs_j denotes the jth entry of the SFS for locus i, then the output format is

  loc_1_sfs_1 loc_1_sfs_2 ... loc_1_sfs_n-1
  
  loc_2_sfs_1 loc_2_sfs_2 ... loc_2_sfs_n-1
  
  ...
  
  loc_m_sfs_1 loc_m_sfs_2 ... loc_m_sfs_n-1

when the simulation consists of m loci and a total sample size (across all islands) of n. Methods for returning similar lists of normalised branch lengths (instead of normalised SFSs) as well as SFSs simulated with a fixed number of mutations (as opposed to a fixed mutation rate) have also been implemented as ancestry.print_normalised_branch_lengths and ancestry.fixed_s, respectively.

# Genetic forces
Currently implemented are crossover recombination within loci, spatial structure with migration on an arbitrary finite graph, and weak natural selection with fitness depending on the full genotype of both parents, as well as their geographical location. Simulating more than one island without migration will result in an infinite loop as the most recent common ancestor of the whole sample cannot be reached, unless only one island has a positive sample size.

# Epochs
The change_times parameter lets the user specify times (in generations) for jump discontinuities in simulation parameters. Relative population sizes across islands, population growth rates, weak selection parameters, and migration rates can vary in this way. All of these parameters must be specified for all epochs (even when they don't change between two epochs) as outlined in dev.cfg. Phylogenetic species trees can be simulated by specifying multiple islands, initially with migration rates set to zero, provided that all islands communicate with positive rate eventually. 
