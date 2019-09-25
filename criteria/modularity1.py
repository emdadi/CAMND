#__author__ = 'fatemeh'
import eval_util
import my_constants
import importer
import my_util
import eval_constants

def preprocess(inf,is_mmod,species):
    sbml_file = r'%s' % (my_constants.species_sbml[species])
    S, mets, rxns, revs, met_names, rxn_names, biomass, met_comparts \
        = importer.sbmlStoichiometricMatrix(sbml_file, True, read_species_compart=True)
    mo=compute_modularity(inf,not is_mmod,False,True,S,mets,rxns)
    return mo

# formula: Sigma(s=1..NumOfModules) [ ls / L - ( ds / 2L )^2 ] where 'ls' is number of links inside module 's' and 'ds' is sum of degree of
# nodes in module 's' and 'L' is the total number of links in the network
# NOTE: ls and ds can be computed normally for rmod_converted networks but L should have a new definition because of overlap between
# metabolites of different modules
def compute_modularity(inf, is_rmod, do_bipartite_explode, is_partial_module, S, mets, reacts):
    mods = eval_util.read_mmods(inf)

    metab_edges = {}

    # the following IF merely fills metabmod_orgranized_reacts
    # TODO: do_bipartite_explode + is_rmod not supported yet! do_bipartite_explode is simply ignored in the case
    if is_rmod:
        mods, metabmod_orgranized_reacts = eval_util.convert_reactmod_to_overlapping_metabmod(mods, is_partial_module, S, mets, reacts)
    elif do_bipartite_explode:
        metabmod_orgranized_reacts = {}

        graph_edges = my_util.graph_by_explode_reactions_to_complete_bipartite(S, mets)

        # the following computes metabmod_orgranized_reacts
        metab_to_module = {v: k for k in mods.keys() for v in mods[k]}
        for e in graph_edges:
            r_mods = {metab_to_module[e[0]], metab_to_module[e[1]]}
            if len(r_mods) == 1:
                mod = r_mods.pop()
                if mod not in metabmod_orgranized_reacts:
                    metabmod_orgranized_reacts[mod] = ([], [], [])
                metabmod_orgranized_reacts[mod][0].append(e)
            else:
                for mod in r_mods:
                    if mod not in metabmod_orgranized_reacts:
                        metabmod_orgranized_reacts[mod] = ([], [], [])
                    metabmod_orgranized_reacts[mod][2].append(e)

    else:
        metabmod_orgranized_reacts = {}

        metab_to_module = {v: k for k in mods.keys() for v in mods[k]}

        if is_partial_module:
            done_metabs = set(metab_to_module.keys())
            for m in mets:
                if m not in done_metabs:
                    metab_to_module[m] = eval_constants.EXTERNAL_MODULE_ID

        for r in reacts:
            r_metabs = eval_util.get_metabolites_of_reaction_idx(S, mets, reacts.index(r))
            r_mods = set([metab_to_module[m] for m in r_metabs])

            # TODO: with the following if, reactions where all nodes are in either one module or external space, are counted as outside-reaction
            if len(r_mods) == 1:  # inner-edge detected
                mod = r_mods.pop()

                if is_partial_module and mod == eval_constants.EXTERNAL_MODULE_ID:  # TODO:PARTIAL: skip reactions completely outside all modules
                    continue

                if mod not in metabmod_orgranized_reacts:
                    metabmod_orgranized_reacts[mod] = ([], [], [])
                metabmod_orgranized_reacts[mod][0].append(r)
            else:  # outer-edge detected
                for mod in r_mods:
                    if mod not in metabmod_orgranized_reacts:
                        metabmod_orgranized_reacts[mod] = ([], [], [])
                    metabmod_orgranized_reacts[mod][2].append(r)

    # compute metab_edges
    if do_bipartite_explode:
        for e in graph_edges:
            for ee in e:
                if ee not in metab_edges:
                    metab_edges[ee] = []
                metab_edges[ee].append(e)
    else:
        for ri, r in enumerate(reacts):
            for re_idx in range(len(S)):
                if S[re_idx][ri] == 0:
                    continue

                m = mets[re_idx]
                if m not in metab_edges:
                    metab_edges[m] = []
                metab_edges[m].append(r)

    # compute total_links
    if do_bipartite_explode:
        total_links = len(graph_edges)
    else:
        total_links = 0
        for m_row in S:  # TODO: what about metabolites that are counted more than one time because of overlapping modules?
            m_degree = sum([1 for ri in m_row if ri != 0])
            total_links += m_degree
        total_links /= 2.0

    modularity = 0.0
    for mod_name, mod_metabs in mods.iteritems():
        ls = len(metabmod_orgranized_reacts[mod_name][0])  # + len(metabmod_orgranized_reacts[mod_name][1])   # TODO: wrongly-inner edges should be counted?

        ds = 0
        for m in mod_metabs:
            # metab_row = S[mets.index(m)]
            # ds += sum([1 for ri in metab_row if ri != 0])
            ds += len(metab_edges[m])  # TODO: for hyperarcs this causes a reaction to be counted more than once for a module

        sum_term = ls * 1.0 / total_links - pow(ds / (2.0 * total_links), 2)
        if mod_name == eval_constants.EXTERNAL_MODULE_ID:
            modularity -= sum_term
        else:
            modularity += sum_term

    return modularity


if __name__ == '__main__':
    # species = 'toy_model'
    import cPickle as pickle
    for species in my_constants.species_sbml.keys():
        # species = 'ecoli_iaf1260'

        out_dir = r'%s/%s/newman' % (my_constants.resultPath, species)

        print species
        S, mets, rxns, revs, met_names, rxn_names, biomass, met_comparts = importer.sbmlStoichiometricMatrix(r'%s/dataset/networks/%s' % (my_constants.basePath, my_constants.species_sbml[species]), True, read_species_compart=True)

        # pickle.dump([S, mets, rxns, revs, met_names, rxn_names, biomass, met_comparts], open(r'D:\University\DiseaseSim\zsh\MSB\%s.pkl' % species, 'wb'))
        # aS, amets, arxns, arevs, amet_names, arxn_names, abiomass, amet_comparts = pickle.load(open(r'D:\University\DiseaseSim\zsh\MSB\ecoli_iaf1260.pkl', 'rb'))

        # print S == aS, mets == amets, revs == arevs, amet_names == met_names, abiomass == biomass, amet_comparts == met_comparts, arxns == rxns
        # print compute_modularity('%s/final_modules.txt' % out_dir, False, True, False, S, mets, rxns)