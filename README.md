# Beta-Xi-Sim
An efficient simulator for ancestries of multiple unlinked loci driven by Beta-Xi-coalescents with exponential population growth. Valid for any number of loci and sample sizes &lt; 10 000.

# Compilation and usage
In Unix-like environments, compile Beta-Xi-Sim by navigating into the directory to which you've downloaded the source files source files, and calling `make`. The programme can then be called as `./simulate path_to_config_file`. An example config file is provided in `dev.cfg`. Any custom config files must conform to the same format!

# User interface
Parameters can be specified via a text file interface (see provided `dev.cfg` for an example).
More powerful customisation, such as arbitrary fitness functions and recombination maps, can be done by specifying the user-defined functions at the top of the Ancestry struct in `ancestry.hh`.

# Output
The software outputs a space-delimited, normalised site frequency spectrum (SFS) for each island and each locus. If isl_i_loc_j_sfs_k denotes the kth entry of the SFS for locus j, on island i, then the output format is

  `isl_1_loc_1_sfs_1 isl_1_loc_1_sfs_2 ... isl_1_loc_1_sfs_n-1`    
  `isl_1_loc_2_sfs_1 isl_1_loc_2_sfs_2 ... isl_1_loc_2_sfs_n-1`  
  `...`  
  `isl_1_loc_m_sfs_1 isl_1_loc_m_sfs_2 ... isl_1_loc_m_sfs_n-1`  
  `isl_2_loc_1_sfs_1 isl_2_loc_1_sfs_2 ... isl_2_loc_1_sfs_n-1`  
  `...`  
  `isl_l_loc_m_sfs_1 isl_l_loc_m_sfs_2 ... isl_l_loc_m_sfs_n-1`

when the simulation consists of `m` loci, `l` islands, and a total sample size (across all islands) of `n`. Methods for returning similar lists of pooled (across islands) normalised SFSs, branch lengths, as well as SFSs simulated with a fixed number of mutations (as opposed to a fixed mutation rate) have also been implemented as `ancestry.print_normalised_sfs`, `ancestry.print_normalised_branch_lengths`, and `ancestry.fixed_s`, respectively.

# Genetic forces
Currently implemented are crossover recombination within loci, spatial structure with migration on an arbitrary finite graph, and weak natural selection with fitness depending on the full genotype of both parents, as well as their geographical location. Simulating more than one island without migration will result in an infinite loop as the most recent common ancestor of the whole sample cannot be reached, unless only one island has a positive sample size.

# Epochs
The change_times parameter lets the user specify times (in generations) for jump discontinuities in simulation parameters. Relative population sizes across islands, population growth rates, weak selection parameters, and migration rates can vary in this way. All of these parameters must be specified for all epochs (even when they don't change between two epochs) as outlined in dev.cfg. Phylogenetic species trees can be simulated by specifying multiple islands, initially with migration rates set to zero, provided that all islands communicate with positive rate eventually. 

# Citation
If you use Beta-Xi-Sim in your work or wish to refer to it, please cite the following accompanying article:

`@article{Koskela:2018:SAGMB:20170011,`  
  `author = "Jere Koskela",`  
  `title = "Multi-locus data distinguishes between population growth and multiple merger coalescents,`  
  `journal = "Statistical Applications in Genetics and Molecular Biology",`  
  `volume = "17",`  
  `number = "3",`  
  `pages = "20170011",`  
  `year = "2018",`  
  `URL = "https://www.degruyter.com/view/j/sagmb.2018.17.issue-3/sagmb-2017-0011/sagmb-2017-0011.xml"`  
`}`
