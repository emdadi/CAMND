SEPARATOR program downloaded from http://pinguin.biologie.uni-jena.de/bioinformatik/networks/index.html

initial external metabolites should be chosen to start!

when started, tries to convert some internal metabolites to external and find CCs

also proper threshold can be chosen which adds all metabolites with degree higher
than the given threshold are considered pseudo-externals (and added to external
 metabolites) and run again


process: convert high degree internal metabolites to external metabolites and find
    connected components THEN WHY MODULES WITH ONLY ONE REACTION AND WITH NO INTERNAL
    METABOLITE? HOW IS A MODULE DEFINED?

implement schilling/schuster (what threshold/initial external metabolites should be used?)

STRANGE:
Sink needed to allow 4-hydroxy-benzoate to leave system
is a reaction for metabolite M_4hba_c which is like exchange but without _e and is not
 named exchange and is irreversible! five of these in ecoli_iar1260
 M_4hba_c_BND M_5drib_c_BND M_aacald_c_BND M_hmfurn_c_BND M_oxam_c_BND



if removing currency metabolites at the start of decomposition is added to the pipeline, then these list can be
used as additional External Metabolites in Schuster decompostion


NOTE ABOUT SCHUSTER PREPROCESSINGS ITSELF:

Applying the decomposition procedure outlined in the previous section and cancelling all subnetworks that only involve external metabolites, 19 subnetworks are obtained (Pfeiffer, 1999). It has turned out that the biochemical interpretability of results is considerably improved by slightly editing the automated classification, based on biochemical knowledge. For example, in some subsys-tems, it is useful to combine the reactions using NADP with those using NAD, in order to avoid combinatorial explosion. Moreover, for enzymes with broad specificity, the functions with low activities have been deleted. For example, uridine kinase (EC 2.7.1.48) can use a wide spectrum of nucleotide triphosphates. However, CTP, dCTP, UTP and dUTP are very poor substrates. 