import shutil

from util import dendogram_util
from util import my_constants
from util import my_util
from util import importer
from util.my_util import write_line


def read_hierarchical_decomposition_sridharan(species, file_path, mets, rxns, is_mets=True):
    # each mh in module_hierarchy is [super_idx, [sub_idxs], depth]
    module_hierarchy, mods_mets, mods_rxns, tree_height = read_module_hierarchy(file_path, mets, rxns)

    raw_observations = mets if is_mets else rxns
    modules = mods_mets if is_mets else mods_rxns

    # first correct module ids to start from 'observation_count + 1' instead of '0'
    raw_observations_count = len(raw_observations)
    corrected_module_hierarchy = {}
    corrected_modules = {}
    for mi, mh in enumerate(module_hierarchy):
        mh_super, mh_subs, mh_depth = mh[0], mh[1], mh[2]

        corrected_id = mi + raw_observations_count + 1
        corrected_super = mh_super + raw_observations_count + 1 if mh_super is not None else None
        corrected_subs = [s + raw_observations_count + 1 for s in mh_subs]

        corrected_module_hierarchy[corrected_id] = (corrected_super, corrected_subs, mh_depth)
        corrected_modules[corrected_id] = modules[mi]
    module_hierarchy = corrected_module_hierarchy
    modules = corrected_modules

    super_sub_map = {}
    observations = []
    for mi, mh in module_hierarchy.iteritems():
        if len(mh[1]) > 0:  # has any child module
            super_sub_map[(mh[2], mi)] = [(module_hierarchy[mhs][2], mhs) for mhs in mh[1]]
        else:
            obser_clust = []
            for o in modules[mi]:
                obser = (mh[2] + 1, raw_observations.index(o) + 1)  # observations should be numbered from 1
                obser_clust.append(obser)
                observations.append(obser)
            super_key = (mh[2], mi)
            if super_key not in super_sub_map:
                super_sub_map[super_key] = []
            super_sub_map[super_key].extend(obser_clust)

    # compute reverse relation: sub_super_map
    sub_super_map = {}
    for super, subs in super_sub_map.iteritems():
        for sub in subs:
            sub_super_map[sub] = super

    # compress simple paths: they only occur for leaf modules (which contain observations)
    for obser in observations:
        super = sub_super_map[obser]
        if len(super_sub_map[super]) == 1:
            real_super_key = sub_super_map[super]
            super_sub_map[real_super_key].remove(super)
            super_sub_map[real_super_key].append(obser)
            del super_sub_map[super]

    if len(observations) != len(raw_observations):
        print 'different number of observations and raw_observations when treating %s (is_mets: %s)' % (file_path, str(is_mets))
        exit(-1)

    # now collapse more than two way branches to two ways by adding more clusters!
    super_sub_map = dendogram_util.collapse_more_than_two_way_branches(super_sub_map)

    # build merge_occurrences according to matplotlib.dendogram specification. each point in super_sub_map
    # exactly specifies one branch point.
    len_map = {}
    cluster_id_map = {}
    observation_count = len(observations)

    for si, s in enumerate(sorted(observations)):
        cluster_id_map[s] = si + 1
        len_map[s] = 1

    merge_occurrences = dendogram_util.build_merge_occurrences(super_sub_map, cluster_id_map, len_map, observation_count)

    final_results = []
    for i in range(len(merge_occurrences)):
        f = merge_occurrences[(observation_count + 1) + i]
        final_results.append([float(f[0] - 1), float(f[1] - 1), float(f[2]), float(f[3])])

    method_dir = r'%s/method/sridharan' % my_constants.basePath
    cts = eval(open('%s/cut_heights.txt' % method_dir, 'r').read())
    thresholds = cts[species]

    return final_results, final_results[-1][2], thresholds


def methods_list_getter(species):
    method_dir = r'%s/method/sridharan' % my_constants.basePath
    cts = eval(open('%s/cut_heights.txt' % method_dir, 'r').read())
    return ['sridharan%s' % l for l in cts[species]]


def convert_to_int_recursive(nested_list):
    if type(nested_list) == list:
        ret = []
        for l in nested_list:
            ret.append(convert_to_int_recursive(l))
        return ret
    else:
        return int(nested_list)


def get_path_to_root_recursive(modules_hierarchy, module_index):
    if modules_hierarchy[module_index][0] is None:
        return [module_index]
    else:
        return [module_index] + (get_path_to_root_recursive(modules_hierarchy, modules_hierarchy[module_index][0]))


def read_module_hierarchy(file_path, mets, rxns):
    res_vars = my_util.try_load_matlab_result(file_path)
    raw_modules = res_vars['mods']
    raw_modules_hierarchy = res_vars['mods_hier']

    mods_rxns = []
    mods_mets = []
    for x, raw_mod in enumerate(raw_modules):
        raw_mod = raw_mod[0]
        mod_rxns_idxs = raw_mod[0][:, 0]
        mod_mets_idxs = raw_mod[1][:, 0]
        # mod_metname[0] = str(raw_mod[2][:, 0][1][0])
        mod_s = raw_mod[3].toarray().tolist()
        # mod_rxnname[0] = str(raw_mod[4][:, 0][1][0])
        mod_bmatrix = raw_mod[5].tolist() if len(raw_mod[5][0]) > 0 else None
        mod_q = raw_mod[6][0][0] if len(raw_mod[6][0]) > 0 else None
        mod_s_eig = raw_mod[7][:, 0] if len(raw_mod[7][0]) > 0 else None

        mod_rxns = [rxns[i - 1] for i in mod_rxns_idxs]
        mod_mets = [mets[i - 1] for i in mod_mets_idxs]

        mods_rxns.append(mod_rxns)
        mods_mets.append(mod_mets)

    raw_modules_hierarchy = convert_to_int_recursive(raw_modules_hierarchy.toarray().tolist())

    modules_hierarchy = [[None, [], -1] for i in range(len(raw_modules))]  # map of module index to (super_index, [sub_indexes], depth). note that depth of root is 1
    for ri, rm in enumerate(raw_modules_hierarchy):
        for di, d in enumerate(rm):
            if d != 0:
                modules_hierarchy[ri][1].append(di)
                modules_hierarchy[di][0] = ri

    # note that tree_height - depth + 1 will give the height of each module. height of deepest modules will be 1! root will be the highest.
    tree_height = -1
    for mi, mh in enumerate(modules_hierarchy):
        depth = len(get_path_to_root_recursive(modules_hierarchy, mi))
        mh[2] = depth
        if depth > tree_height:
            tree_height = depth

    # NOTE: this important hack causes leaf modules (which are real outputs of sridharan be placed on the bottom of the
    # dendrogram, nearest place to their members (i.e. observations)
    leaf_modules_idxs = []
    for mi, mh in enumerate(modules_hierarchy):
        if len(mh[1]) == 0:
            leaf_modules_idxs.append(mi)

    for mi in leaf_modules_idxs:
        modules_hierarchy[mi][2] = tree_height

    return modules_hierarchy, mods_mets, mods_rxns, tree_height


def read_module_hierarchy_incorrect(file_path, mets, rxns):
    res_vars = my_util.try_load_matlab_result(file_path)
    raw_modules = res_vars['mods']
    raw_modules_hierarchy = res_vars['mods_hier']

    mods_rxns = []
    mods_mets = []
    for x, raw_mod in enumerate(raw_modules):
        raw_mod = raw_mod[0]
        mod_rxns_idxs = raw_mod[0][:, 0]
        mod_mets_idxs = raw_mod[1][:, 0]
        # mod_metname[0] = str(raw_mod[2][:, 0][1][0])
        mod_s = raw_mod[3].toarray().tolist()
        # mod_rxnname[0] = str(raw_mod[4][:, 0][1][0])
        mod_bmatrix = raw_mod[5].tolist() if len(raw_mod[5][0]) > 0 else None
        mod_q = raw_mod[6][0][0] if len(raw_mod[6][0]) > 0 else None
        mod_s_eig = raw_mod[7][:, 0] if len(raw_mod[7][0]) > 0 else None

        mod_rxns = [rxns[i - 1] for i in mod_rxns_idxs]
        mod_mets = [mets[i - 1] for i in mod_mets_idxs]

        mods_rxns.append(mod_rxns)
        mods_mets.append(mod_mets)

    raw_modules_hierarchy = convert_to_int_recursive(raw_modules_hierarchy.toarray().tolist())

    modules_hierarchy = [[None, [], -1] for i in range(len(raw_modules))]  # map of module index to (super_index, [sub_indexes], depth). note that depth of root is 1
    for ri, rm in enumerate(raw_modules_hierarchy):
        for di, d in enumerate(rm):
            if d != 0:
                modules_hierarchy[ri][1].append(di)
                modules_hierarchy[di][0] = ri

    # note that tree_height - depth + 1 will give the height of each module. height of deepest modules will be 1! root will be the highest.
    tree_height = -1
    for mi, mh in enumerate(modules_hierarchy):
        depth = len(get_path_to_root_recursive(modules_hierarchy, mi))
        mh[2] = depth
        if depth > tree_height:
            tree_height = depth

    return modules_hierarchy, mods_mets, mods_rxns, tree_height


def go(species, only_cut_dendogram=False):
    need_biomass_removal = False and my_constants.species_artificial_biomass[species]

    method_dir = r'%s/method/sridharan' % my_constants.basePath
    out_dir = r'%s/%s/sridharan' % (my_constants.resultPath, species)
    my_util.mkdir_p(out_dir)

    source_file = '%s/dataset/networks/%s' % (my_constants.basePath, my_constants.species_sbml[species])

    S, mets, rxns, revs, met_names, rxn_names, biomass, met_comparts = importer.sbmlStoichiometricMatrix(source_file, True, read_species_compart=True, remove_biomass=need_biomass_removal, normalize_stoich=False)

    if not only_cut_dendogram:
        f = open('sridharan.m', 'w')
        write_line(f, 'addpath %s/code' % method_dir)
        write_line(f, "model = readCbModel('%s')" % source_file)
        if need_biomass_removal:
            for l in my_util.get_remove_reaction_matlab_script('model', biomass):
                write_line(f, l)
        write_line(f, '[mods, mods_hier] = Shred_Network_plos2011(model, false);')
        write_line(f, 'save sridharan.mat mods mods_hier')

        f.close()

        my_util.prepare_matlab_file_and_exec_and_wait_finish('sridharan', 'sridharan.mat', False)

        shutil.copy('sridharan.mat', out_dir)

    modules_hierarchy, mods_mets, mods_rxns, raw_tree_height = read_module_hierarchy('%s/sridharan.mat' % out_dir, mets, rxns)

    # TODO: very bad job!!!!
    # NOTE: mets are removed because their modules share metabolites!!!
    dummy_tree, tree_height, dummy_thresholds = read_hierarchical_decomposition_sridharan(species, '%s/sridharan.mat' % out_dir, mets, rxns, is_mets=False)

    def dendogram_height(hierarchy_height):
        return tree_height - hierarchy_height + 1

    cut_heights = eval(open('%s/cut_heights.txt' % method_dir, 'r').read())
    for thrr in cut_heights[species]:

        # TODO: very bad job
        thr_height = thrr * tree_height

        all_modularized_observations = set()
        cut_mods = []
        # to cut, find all modules above below cut height whose super above cut height
        # also all other observations are added as remaining modules (observations are not in module_hierarchy as one-sized-module and so should be specially treated)
        for m_idx, m_sup_sub in enumerate(modules_hierarchy):
            module_height = m_sup_sub[2]
            super_module_height = modules_hierarchy[m_sup_sub[0]][2] if m_sup_sub[0] is not None else tree_height + 1
            if dendogram_height(module_height) <= thr_height < dendogram_height(super_module_height):
                cut_mods.append(mods_rxns[m_idx])
                all_modularized_observations.update(mods_rxns[m_idx])

        for r in rxns:
            if r not in all_modularized_observations:
                cut_mods.append([r])

        outr = open('%s/react_modules_%s.txt' % (out_dir, thrr), 'w')
        # outm = open('%s/metab_modules_%s.txt' % (out_dir, thrr), 'w')
        # for rmod, mmod in zip(cut_mods_rxns, cut_mods_mets):
        for rmod in cut_mods:
            write_line(outr, ' '.join(rmod))
            # write_line(outm, ' '.join(mmod))
        outr.close()
        # outm.close()

        shutil.copy('%s/react_modules_%s.txt' % (out_dir, thrr), '%s/final_modules_%s.txt' % (out_dir, thrr))  # TODO: which file r/m should be selected as final modules?


if __name__ == '__main__':
    go('ecoli_core')
    # go('helico_iit341')
    # go('ecoli_ijr904')
    # go('ecoli_iaf1260')
    # go('ecoli_ijo1366')
    # go('saccaro_ind750')
    # go('homo_recon1')

    # dist = np.array([
    # [0,1,100,100],
    # [1,0,100,100],
    #     [100,100,0,1],
    #     [100,100,1,0]
    # ], dtype=float)
    # print wpgma(dist)