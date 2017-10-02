# Beta-Xi-Sim
An efficient simulator for ancestries of multiple unlinked loci driven by Beta-Xi-coalescents, as well as Kingman coalescents with algebraic or exponential population growth. Valid for any number of loci and sample sizes &lt; 10 000.

# User interface
Parameters can be specified via a text file interface (see provided dev.cfg for an example).
More powerful customisation, such as arbitrary fitness functions and recombination maps, can be done by specifying the user-defined functions at the top of the Ancestry struct in ancestry.hh.

# Genetic forces
Currently implemented are crossover recombination within loci, spatial structure with migration on an arbitrary finite graph, and weak natural selection with fitness depending on the full genotype of both parents, as well as their geographical location.
