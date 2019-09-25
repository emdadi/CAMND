# CAMND
Comparative Analysis of Metabolic Network Decomposition

 These codes provide a package for decomposing of metabolic networks with 10 diﬀerent methods, 9 datasets and the ability of computing
 12 criteria, to evaluate and compare the results of diﬀerent methods using ten previously deﬁned and two new criteria which are based
on Chebi ontology and Co-expression of Enzymes information.

The list of methods for decomposition of metabolic networks is as follows:
Ding,
Muller2,
Schuster,
Holme,
Guimera,
Newman,
Poolman,
Sridharan,
Muller,
Verwoerd. The files corresponding to each method includes the implemented codes of them.
The list of criteria for evaluate methods is as follow:
cohesion_coupling
efficacy,
go_distance_bp_F,
go_distance_bp_G,
go_distance_cc_F,
go_distance_cc_G,
go_distance_mf_F,
go_distance_mf_G,
modularity,
module_count,
coexpression_of_enzymes,
chebi_distance. These criteria implemented in "critera" file.
