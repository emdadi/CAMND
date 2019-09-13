[m, mh] = Shred_Network_plos2011(model, false)
m: each row one module. of course many modules are 1-react, 1-metab. do not be afraid
mh: hierarchical relations between modules defined in m:

TRUE instead of FALSE as the second argument, can take regulatory information in calculations.
geryk needs them, I think there are others that can use. seems that there is a way to extract (which
is noted in sridharan paper) so why not use?


Dear Abolfazl,

I apologize for the delay. I was on vacation and I simply forgot to get back to you. I have not finished completing the documentation in an 'official manner', but as you said I'll just tell you how to use it. 

You first want to create a struct in matlab for your model; call it 'Model'. Let's say your Model has N reactions and M metabolites. The elements of the struct are 

Model.S = stoichiometric matrix (metabolite rows and reaction columns) which is a 'double' variable of MXN dimensions. 
Model.rxns = NX1 cell array with the names of the reactions
Model.rev = 1XN double vector whose entries are either 0 or 1 for whether or not the reaction is reversible
Model.reg = regulatory matrix MXN whose entries are 0 or 1. 1 if the metabolite allosterically regulates the reaction. 
Model.mets = MX1 cell array of the metabolite names. 

Unzip the files I have sent and you want to run [Module, module_hierarchy] = Shred_Network_plos11(Model, true) - the true is for whether or not you want to include regulatory information. 

The Module output is a struct with information about each module that is generated in the hierarchy. These include the reactions present, the mets present, the module S matrix, the b-matrix, the modularity Q score, and the approximate eigenvector. If the module has more than one connected component, then it just breaks apart into those components and that module will obviously not have a Q score or b-matrix. 

The module_hierarchy is a connectivity matrix for the modules, so you can find the parent and children of each module in the hierarchy. 

Let me know if this works and if you have any further questions! Good luck.


CUT_HEIGHT: when a tree is cut at height H you should assume that all nodes with height less than H are washed out and then
 remaining leaves are the answer. 1 will be height of leaves so that cutting at height 0.5 will give all leaf modules.
 cutting at TREE_HEIGHT - 0.5 will give one module containing all nodes. cutting at TREE_HEIGHT+ will give nothing!





 Dear Gautham,

 I'm now executing network partitioning algorithms I've obtained so far. I have question about your method: should I remove biomass reaction before execution or it is better not to be removed? Does your method rely on biomass reaction?

 Regards.



 Hi Abolfazl,

 It really depends on the application. Biomass reactions involve many substrates and if you want to include the interactions between biomass and central carbon metabolism, you should include it. As long as you clearly state which reactions you include, you should be ok.

 Regards,
 Gautham