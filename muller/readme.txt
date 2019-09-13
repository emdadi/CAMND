Dear Abolfazl,

currently there exist two implementations for the flux modules computations.

1. http://sourceforge.net/projects/fluxmodules/
This is implemented in matlab. It uses a more efficient algorithm than the one described in the paper that you have read. I attached the publication describing the faster algorithm.

2. http://sourceforge.net/projects/cbmpy/
This is actually a metabolic network analysis toolbox written in python. You can find the fluxmodules code in the package fluxmodules. This implementation is done in rational arithmetic, since I encountered a number of numerical inprecisions in the other version. This implementation contains even more features and algorithms for computing k-modules, a generalization of the flux-modules concept. For more details I refer to my arxiv preprint: 
http://arxiv.org/abs/1404.5584

Best,
Arne