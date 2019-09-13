import pickle
import math

import scipy as sp
import scipy.linalg as linalg
import numpy as np


# matplotlib.use("PDF")
from scipy import matrix
import scipy.cluster.hierarchy as cluster
from util import my_constants
from util import importer
from util import my_util

CONSIDER_ZERO_THRESHOLD = 1e-10


def read_hierarchical_decomposition_poolman(species, file_path):
    tree = pickle.load(open(file_path, 'r'))

    method_dir = r'%s/method/poolman' % my_constants.basePath
    cts = eval(open('%s/cut_heights.txt' % method_dir, 'r').read())
    thresholds = cts[species]

    return tree, tree[-1][2], thresholds


def methods_list_getter(species):
    method_dir = r'%s/method/poolman' % my_constants.basePath
    cts = eval(open('%s/cut_heights.txt' % method_dir, 'r').read())
    return ['poolman%s' % l for l in cts[species]]


def append_zero(m):
    return np.append(m, matrix(np.zeros(m.shape[1]), dtype=m.dtype), axis=0)


def print_s(s):
    i = 0
    while True:
        if i < len(s):
            x = str(s[i]) + '                      '
            print x[:24],
            i += 1
            if i % 5 == 4:
                print
        else:
            break


def null(m, eps):
    print 'WARNING: proper value for eps should be set to find orthogonal null-space. look at values of S to choose :('

    m = append_zero(m)
    u, s, vh = linalg.svd(m)
    # print_s(s)
    null_mask = (s <= eps)
    null_space = sp.compress(null_mask, vh, axis=0)
    return sp.transpose(null_space)


# def gram_schmidt(x, inplace=False):
# # Returns the Gramm-Schmidt orthogonalization of matrix X
# if not inplace:
# v = [row[:] for row in x]  # make a copy.
# else:
# v = x
# k = len(x[0])  # number of columns.
# n = len(x)  # number of rows.
#
# for j in range(k):
# for i in range(j):
# # D = < Vi, Vj>
# d = sum([v[p][i] * v[p][j] for p in range(n)])
#
# for p in range(n):
# # Note that the Vi's already have length one!
# # Vj = Vj - <Vi,Vj> Vi/< Vi,Vi >
#                 v[p][j] -= (d * v[p][i])
#
#         # Normalize column V[j]
#         invnorm = 1.0 / math.sqrt(sum([(v[p][j]) ** 2 for p in range(n)]))
#         for p in range(n):
#             v[p][j] *= invnorm
#     return v


def gram_schmidt(x, inplace=False):
    # Returns the Gramm-Schmidt orthogonalization of matrix X
    if not inplace:
        v = x.copy()  # make a copy.
    else:
        v = x
    k = x.shape[1]  # number of columns.
    n = x.shape[0]  # number of rows.

    for j in range(k):
        for i in range(j):
            # d = < Vi, Vj>
            d = sum([v[p, i] * v[p, j] for p in range(n)])

            for p in range(n):
                # Note that the Vi's already have length one!
                # Vj = Vj - <Vi,Vj> Vi/< Vi,Vi >
                v[p, j] -= (d * v[p, i])

        # Normalize column V[j]
        invnorm = 1.0 / math.sqrt(sum([(v[p, j]) ** 2 for p in range(n)]))
        for p in range(n):
            v[p, j] *= invnorm
    return v


def orthogonal_null(m, eps):
    nm = null(m, eps)
    # TODO: surely some transpose needed :(
    onm = gram_schmidt(nm)
    return onm


def build_distance(stoich, eps):
    rc = len(stoich[0])
    stoich = matrix(stoich, dtype=float)
    phi = sp.zeros((rc, rc), dtype=float)
    on = orthogonal_null(stoich, eps)
    for i in xrange(rc):
        for j in xrange(i + 1, rc):
            phi[i, j] = math.fabs(on[i].dot(on[j].T) / (math.sqrt(on[i].dot(on[i]) * on[j].dot(on[j]))))

    phi[phi < CONSIDER_ZERO_THRESHOLD] = CONSIDER_ZERO_THRESHOLD
    dist = 1 / phi      # TODO:LATER: dist(i,j) = arccos( phi(i,j) ) then why we use 1 / phi(i,j) instead?

    # NOTE: this is a dummy code for compensating all-zero rows in "on" which should not exist!
    for i in xrange(rc):
        if on[i].dot(on[i]) < CONSIDER_ZERO_THRESHOLD:  # all-zero row in "on"
            for j in xrange(i + 1, rc):     # for all elements that are for that row (reaction) do (of course, only over diagonal)
                dist[i, j] = 1 / CONSIDER_ZERO_THRESHOLD    # put infinite distance to other rows (reactions)

    return dist


# def wpgma(dist):
#
#     tree_nodes = {}  # map of label: {other_node_label: dist}
#
#     internal = 0
#     while len(tree_nodes) > 1:
#         i, j = nearest_neighbors(tree_nodes)
#
#         delta = tree_nodes[i][j]
#
#         new_node_dist = {}
#         for l, ds in tree_nodes.iteritems():
#             if l in [i, j]:
#                 continue
#             new_node_dist[l] = (tree_nodes[i][l] + tree_nodes[j][l]) / 2
#
#         tree_nodes['p%d' % internal] = new_node_dist
#         internal += 1
#
#         del tree_nodes[i]
#         del tree_nodes[j]
#

def wpgma(dist):
    dist_condense = []
    for i in xrange(len(dist)):
        dist_condense += dist[i, i + 1:].tolist()
    z = cluster.linkage(dist_condense, 'weighted')
    return z


def go(species, only_cut_dendogram=False):
    method_dir = r'%s/method/poolman' % my_constants.basePath
    out_dir = r'%s/%s/poolman' % (my_constants.resultPath, species)
    my_util.mkdir_p(out_dir)

    S, mets, rxns, revs, met_names, rxn_names, biomass, met_comparts = importer.sbmlStoichiometricMatrix(r'%s/dataset/networks/%s' % (my_constants.basePath, my_constants.species_sbml[species]), True, read_species_compart=True, normalize_stoich=True, remove_biomass=True)

    if not only_cut_dendogram:
        epses = eval(open('%s/EPSes.txt' % method_dir, 'r').read())
        dist = build_distance(S, epses[species])
        tree = wpgma(dist)

        # show dendogram
        # plt.title("dendogram")
        # cluster.dendrogram(tree, show_leaf_counts=True)
        ## plt.show()
        # plt.savefig('%s/dendogram.png' % out_dir)

        # outputs: wpgma tree, cut at proper height?
        out = open('%s/dendogram.py' % out_dir, 'w')
        out.write("#print 'A 4 by (n-1) matrix Z is returned. At the i-th iteration, clusters with indices Z[i, 0] and Z[i, 1]'\n")
        out.write("#print 'are combined to form cluster n+i. A cluster with an index less than n corresponds to one of the'\n")
        out.write("#print 'n original observations. The distance between clusters Z[i, 0] and Z[i, 1] is given by Z[i, 2].'\n")
        out.write("#print 'The fourth value Z[i, 3] represents the number of original observations in the newly formed cluster'\n")
        out.write("\n")
        out.write("tree = " + str(tree.tolist()))
        out.close()

        pf = open('%s/pickled_tree.obj' % out_dir, 'w')
        pickle.dump(tree, pf)
        pf.close()

    else:
        tree = pickle.load(open('%s/pickled_tree.obj' % out_dir, 'r'))

    tree_height = tree[-1][2]   # height of root

    cts = eval(open('%s/cut_heights.txt' % method_dir, 'r').read())
    for l in cts[species]:

        # TODO: very bad job
        # iter = int(l * len(levels))
        thr = l * tree_height

        res = {}
        rxns_idx_to_mod = cluster.fcluster(tree, thr, criterion='distance')
        for i in range(len(rxns_idx_to_mod)):
            if rxns_idx_to_mod[i] not in res:
                res[rxns_idx_to_mod[i]] = []
            res[rxns_idx_to_mod[i]].append(rxns[i])

        out = open('%s/final_modules_%s.txt' % (out_dir, l), 'w')
        for mod in res.itervalues():
            out.write(' '.join(mod))
            out.write('\n')
        out.close()


if __name__ == '__main__':
    # go('toy_model')
    # go('ecoli_core')
    # go('helico_iit341')
    # go('ecoli_iaf1260')
    # go('saccaro_ind750')
    # go('ecoli_ijo1366')
    # go('homo_recon1')
    go('arabidopsis_irs1597')

    # dist = np.array([
    #     [0,1,100,100],
    #     [1,0,100,100],
    #     [100,100,0,1],
    #     [100,100,1,0]
    # ], dtype=float)
    # print wpgma(dist)