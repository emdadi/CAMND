import gc
import copy
import os
import shutil
from operator import itemgetter

from util import my_constants
from util import importer
from util import my_util
from util import dendogram_util


def read_hierarchical_decomposition_holme(species, file_path):
    # file_path = r'e:\wintemp\yyy.txt'

    s = open(file_path, 'r').read()
    tree_idx = s.index('tree')
    st_idx = s.find('[', tree_idx)
    ed_idx = s.rfind(']', st_idx)
    tree = my_util.read_standard_written_2d_list(s, st_idx, ed_idx, int)   # eval(s[s.index('=', tree_idx) + 1:])

    add_first_all_one_row = False
    for f in tree[0]:
        if f != 1:
            add_first_all_one_row = True
    if add_first_all_one_row:
        tree.insert(0, [1 for i in range(len(tree[0]))])

    # build super_sub_map and sub_super_map for the final tree to show. label of each node is (depth, number)
    # temp_X_Y_map means initial map
    # semi_normalized_X_Y_map means normalized but simple paths are present
    # X_Y_map means final version!

    temp_sub_super_map = {}
    for depth in range(1, len(tree)):
        row = tree[depth]
        for ri, r in enumerate(row):
            temp_sub_super_map[(depth, r)] = (depth - 1, tree[depth - 1][ri])

    temp_super_sub_map = {}
    for sub, super in temp_sub_super_map.iteritems():
        if super not in temp_super_sub_map:
            temp_super_sub_map[super] = []
        temp_super_sub_map[super].append(sub)

    print '\tgeting compressed paths ...',
    gc.collect()

    semi_normalized_sub_super_map = {}

    indirect_real_super = copy.deepcopy(temp_sub_super_map)
    sorted_keys = sorted(temp_super_sub_map.keys(), key=itemgetter(0))  # sort by level
    for super in sorted_keys:
        subs = temp_super_sub_map[super]
        if len(subs) == 1 and super in temp_sub_super_map:
            indirect_real_super[subs[0]] = indirect_real_super[super]

    for sub, super in temp_sub_super_map.iteritems():
        # semi_normalized_sub_super_map[sub] = get_super_in_compressed_tree(sub, temp_sub_super_map, temp_super_sub_map)
        semi_normalized_sub_super_map[sub] = indirect_real_super[sub] if sub in indirect_real_super else super

    del temp_sub_super_map
    gc.collect()

    print 'done!'

    semi_normalized_super_sub_map = {}
    for sub, super in semi_normalized_sub_super_map.iteritems():
        if super not in semi_normalized_super_sub_map:
            semi_normalized_super_sub_map[super] = []
        semi_normalized_super_sub_map[super].append(sub)

    sub_super_map = {}
    for super, subs in semi_normalized_super_sub_map.iteritems():
        for sub in subs:
            if sub not in temp_super_sub_map or len(temp_super_sub_map[sub]) > 1:
                sub_super_map[sub] = super

    del semi_normalized_sub_super_map
    del semi_normalized_super_sub_map
    del temp_super_sub_map
    gc.collect()

    super_sub_map = {}
    for sub, super in sub_super_map.iteritems():
        if super not in super_sub_map:
            super_sub_map[super] = []
        super_sub_map[super].append(sub)

    # deals with special case of a simple path to root!
    root_subs = super_sub_map[(0, 1)]
    if len(root_subs) == 1:
        sub_super_map.pop(root_subs[0])
        super_sub_map.pop((0, 1))

    for s in super_sub_map.keys():
        super_sub_map[s] = sorted(super_sub_map[s], key=itemgetter(1))  # sort children of each module by number (depth of all are the same)

    print 'build sub<->super map done!'

    # now collapse more than two way branches to two ways by adding more clusters!
    super_sub_map = dendogram_util.collapse_more_than_two_way_branches(super_sub_map)

    # build merge_occurrences according to matplotlib.dendogram specification. each point in super_sub_map
    # exactly specifies one branch point.
    len_map = {}
    cluster_id_map = {}
    leaves_depth = len(tree) - 1
    observation_count = len(tree[-1])

    for s in tree[leaves_depth]:
        cluster_id_map[(leaves_depth, s)] = s
        len_map[(leaves_depth, s)] = 1

    merge_occurrences = dendogram_util.build_merge_occurrences(super_sub_map, cluster_id_map, len_map, observation_count)

    final_results = []
    for i in range(len(merge_occurrences)):
        f = merge_occurrences[(observation_count + 1) + i]
        final_results.append([float(f[0] - 1), float(f[1] - 1), float(f[2]), float(f[3])])

    method_dir = r'%s/method/holme' % my_constants.basePath
    cts = eval(open('%s/cut_iterations.txt' % method_dir, 'r').read())
    thresholds = cts[species]

    return final_results, final_results[-1][2], thresholds


def methods_list_getter(species):
    method_dir = r'%s/method/holme' % my_constants.basePath
    cts = eval(open('%s/cut_iterations.txt' % method_dir, 'r').read())
    return ['holme%s' % l for l in cts[species]]


def read_int(f):
    rb = f.read(4)
    if rb == '':
        return None
    return ord(rb[0]) + ord(rb[1]) * 256 + ord(rb[2]) * 256 * 256 + ord(rb[3]) * 256 * 256 * 256


def go(species, only_cut_dendogram=False):
    method_dir = r'%s/method/holme' % my_constants.basePath
    out_dir = r'%s/%s/holme' % (my_constants.resultPath, species)
    my_util.mkdir_p(out_dir)

    S, mets, rxns, revs, met_names, rxn_names, biomass, met_comparts = importer.sbmlStoichiometricMatrix(r'%s/dataset/networks/%s' % (my_constants.basePath, my_constants.species_sbml[species]), True, read_species_compart=True, remove_biomass=True)

    if not only_cut_dendogram:
        # graph is a bipartite representation of the network. should be written to file with format:
        # (1.) each line represents a directed link: 'from' 'to'
        # (2.) substances are enumerated 1, 2, . . .
        # (3.) reaction nodes are enumerated 1000000, 1000001, . . .
        METAB_OFFS = 1
        REACT_OFFS = 1000000

        lines = []
        for i, row in enumerate(S):
            for j, col in enumerate(S[i]):
                if col < 0:  # metab i is consumed by react j
                    lines.append((i + METAB_OFFS, j + REACT_OFFS))
                elif col > 0:  # metab i is produced by react j
                    lines.append((j + REACT_OFFS, i + METAB_OFFS))

        my_util.mkdir_p('inp')
        shutil.copy('%s/inp/cmds' % method_dir, 'inp/')
        shutil.copy('%s/orgnames' % method_dir, './')

        my_util.mkdir_p('cell')

        lines.sort()
        inf = open('cell/ho.dat', 'w')
        for l in lines:
            inf.write('%d %d\n' % (l[0], l[1]))
        inf.close()

        inf = open('cell/ho.nam', 'w')
        for l in mets:
            inf.write('%s\n' % l)
        inf.close()

        # hi cmd_file_name
        # cmd_file_name will be a filename located in inp/
        my_util.mkdir_p('data')
        res = os.system('hi cmds')
        shutil.copy('data/cellho', out_dir)

    resf = open('%s/cellho' % out_dir, 'rb')
    levels = []
    finish = False
    while not finish:
        levels.append([])
        for i in range(len(mets)):
            x = read_int(resf)
            if x is None:
                finish = True
                break
            levels[-1].append(x)
    levels.pop()
    resf.close()

    if not only_cut_dendogram:
        # outputs: wpgma tree, cut at proper height?
        out = open('%s/dendogram.py' % out_dir, 'w')
        out.write("#print tree is 2d list. each entry is the result of algorithm in one iteration. for each iteration there is a list of cluster-index for each metabolite\n")
        out.write("#print e.g. for the first level all values are 1 (meaning that no split is still done) and in the last level metabolites are numbered from 1 to len(mets)\n")
        out.write("\n")
        out.write("tree = " + str(levels))
        out.close()

    # TODO: very bad job!!!!
    dummy_tree, tree_height, dummy_thresholds = read_hierarchical_decomposition_holme(species, '%s/dendogram.py' % out_dir)

    cts = eval(open('%s/cut_iterations.txt' % method_dir, 'r').read())
    for l in cts[species]:

        # TODO: very bad job
        # iter = int(l * len(levels))
        iter = int(tree_height - int(l * tree_height) + 1)

        level_iter = levels[iter]
        modules = {}
        for i, mmod in enumerate(level_iter):
            if mmod not in modules:
                modules[mmod] = []
            modules[mmod].append(i)

        out = open('%s/final_modules_%s.txt' % (out_dir, l), 'w')
        for midx in range(1, len(modules) + 1):
            out.write(' '.join([mets[s] for s in modules[midx]]))
            out.write('\n')
        out.close()

if __name__ == '__main__':
    go('ecoli_core')
    # go('toy_model')
