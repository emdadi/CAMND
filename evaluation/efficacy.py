import math

import eval_util
from util import my_constants
from util import importer


alpha = 1.0


def efficacy_empirical_p(N):
    return 0.25 * math.sqrt(N)


def efficacy_f(n, alpha, p):
    return alpha * math.pow(n, p)


def efficacy_F(n_vec, vec_sz, alpha, p):
    return efficacy_f(vec_sz, alpha, p) + sum([efficacy_f(n, alpha, p) for n in n_vec]) * 1.0 / vec_sz


def compute_efficacy(inf, is_rmod, is_partial_module, S, mets, reacts):
    if is_rmod:
        mods = eval_util.read_mmods(inf)
    else:
        mods = eval_util.read_rmods(inf)

    sizes = []
    for mod_name, mod_metabs in mods.iteritems():
        sizes.append(len(mod_metabs))


    if is_partial_module:
        pass
        # TODO:PARTIAL: should create one big module of unused ones!? NO! DO NOTHING JUST IGNORE THEM!

    N = sum(sizes)

    p = efficacy_empirical_p(N)

    if len(sizes) != 0:
        efficacy = 1.0 * (math.log(efficacy_f(N, alpha, p)) - math.log(efficacy_F(sizes, len(sizes), alpha, p))) * 1.0 / (math.log(efficacy_f(N, alpha, p)) - math.log(2 * efficacy_f(math.sqrt(N), alpha, p)))
    else:
        efficacy = 0.0

    return efficacy


    # metab_edges = {}
    #
    #     metabmod_orgranized_reacts = {}
    #
    #     metab_to_module = {v: k for k in mods.keys() for v in mods[k]}
    #
    #     if is_partial_module:
    #         done_metabs = set(metab_to_module.keys())
    #         for m in mets:
    #             if m not in done_metabs:
    #                 metab_to_module[m] = eval_constants.EXTERNAL_MODULE_ID
    #
    #
    #
    #
    #
    #     for ri, r in enumerate(reacts):
    #         for re_idx in range(len(S)):
    #             if S[re_idx][ri] == 0:
    #                 continue
    #
    #             m = mets[re_idx]
    #             if m not in metab_edges:
    #                 metab_edges[m] = []
    #             metab_edges[m].append(r)
    #
    #
    #
    #
    #
    # # compute total_links
    #     total_links = 0
    #     for m_row in S:  # TODO: what about metabolites that are counted more than one time because of overlapping modules?
    #         m_degree = sum([1 for ri in m_row if ri != 0])
    #         total_links += m_degree
    #     total_links /= 2.0
    #
    #
    #
    #
    # modularity = 0.0
    # for mod_name, mod_metabs in mods.iteritems():
    #     ls = len(metabmod_orgranized_reacts[mod_name][0])  # + len(metabmod_orgranized_reacts[mod_name][1])   # TODO: wrongly-inner edges should be counted?
    #
    #     ds = 0
    #     for m in mod_metabs:
    #         # metab_row = S[mets.index(m)]
    #         # ds += sum([1 for ri in metab_row if ri != 0])
    #         ds += len(metab_edges[m])  # TODO: for hyperarcs this causes a reaction to be counted more than once for a module
    #
    #     modularity += ls * 1.0 / total_links - pow(ds / (2.0 * total_links), 2)
    #
    # return modularity


if __name__ == '__main__':
    # species = 'toy_model'
    species = 'ecoli_core'

    out_dir = r'%s/%s/newman' % (my_constants.resultPath, species)

    S, mets, rxns, revs, met_names, rxn_names, biomass, met_comparts = importer.sbmlStoichiometricMatrix(r'%s/dataset/networks/%s' % (my_constants.basePath, my_constants.species_sbml[species]), True, read_species_compart=True)

    print compute_efficacy('%s/final_modules.txt' % out_dir, False, False, S, mets, rxns)
