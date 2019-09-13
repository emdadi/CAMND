import copy

import eval_constants
from util import my_constants
from util import importer
from util import my_util
from util.my_util import write_line
import eval_util


def mat_rows(l, lj=0):
    m = l[0]
    if lj > 0:
        l[0] = "'%s'" % l[0][:lj].ljust(lj)
        x = reduce(lambda x, y: str(x) + "; '" + (str(y)[:lj].ljust(lj)) + "'", l)
    else:
        l[0] = str(l[0])
        x = reduce(lambda x, y: str(x) + '; ' + str(y), l)
    l[0] = m
    return x


couplings_cache = []
# couplings: 0 - uncoupled , 1 - fully coupled , 2 - partially coupled , 3 - reaction i is directionally coupled to j , 4 - reaction j is directionally coupled to i
# blocks: 1: blocked, 0: ok
def compute_couplings(S, mets, rxns, revs, REACT_MAX_LEN, METAB_MAX_LEN):
    global couplings_cache

    for cc in couplings_cache:
        if S == cc[0] and mets == cc[1] and rxns == cc[2] and revs == cc[3] and REACT_MAX_LEN == cc[4] and METAB_MAX_LEN == cc[5]:
            return cc[6], cc[7]
    couplings_cache.append([S, mets, rxns, revs, REACT_MAX_LEN, METAB_MAX_LEN])

    f = open('cohesion_coupling.m', 'w')

    s_str = 'thenet.stoichiometricMatrix = zeros(%d,%d);' % (len(mets), len(rxns))
    write_line(f, s_str)

    for i in xrange(len(S)):
        s_str = 'thenet.stoichiometricMatrix(%d,:) = %s;' % (i + 1, str(S[i]))
        write_line(f, s_str)

    rev_str = 'thenet.reversibilityVector = [' + mat_rows(revs) + '];'
    write_line(f, rev_str)

    react_str = 'thenet.Reactions = [' + mat_rows(rxns, REACT_MAX_LEN) + '];'
    write_line(f, react_str)

    met_str = 'thenet.Metabolites = [' + mat_rows(mets, METAB_MAX_LEN) + '];'
    write_line(f, met_str)

    # mod_path = r'addpath %s/evaluation/fcf/ffca'
    # func = 'FFCA'
    mod_path = r'addpath %s/evaluation/fcf/f2c2'
    func = 'F2C2'
    fcf_tool_path = mod_path % my_constants.basePath
    write_line(f, fcf_tool_path)

    compute_str = "[cpls, blk] = %s('clp', thenet);" % func
    write_line(f, compute_str)

    save_str = 'save cohesion_coupling.mat cpls blk'
    write_line(f, save_str)

    f.close()

    my_util.prepare_matlab_file_and_exec_and_wait_finish('cohesion_coupling', 'cohesion_coupling.mat', False)

    res_vars = my_util.try_load_matlab_result('cohesion_coupling.mat')
    # if matlab has failed, this will throw exception!
    couplings = res_vars['cpls'].tolist()
    blocks = [b[0] for b in res_vars['blk'].tolist()]

    couplings_cache[-1].append(couplings)
    couplings_cache[-1].append(blocks)
    return couplings, blocks


def compute_module_coupling_table(couplings, blocks, mods, reacts):
    # fully, partially, i-to-j, j-to-i, uncoupled
    module_coupling_table = [[[0, 0, 0, 0, 0] for i in range(len(mods))] for j in range(len(mods))]
    module_order_in_header = sorted(list(mods.keys()))

    # TODO:BUG: does not work for overlapping modules, a reaction may be in several modules
    react_to_module = {v: module_order_in_header.index(k) for k in mods.keys() for v in mods[k]}

    for i in range(len(reacts)):
        for j in range(i + 1, len(reacts)):
            mi = react_to_module[reacts[i]]
            mj = react_to_module[reacts[j]]

            if couplings[i][j] == -1:  # blocked
                pass
            elif couplings[i][j] == 1:  # fully
                if mi == mj:
                    module_coupling_table[mi][mi][0] += 1
                else:
                    module_coupling_table[mi][mj][0] += 1
                    module_coupling_table[mj][mi][0] += 1
            elif couplings[i][j] == 2:  # partially
                if mi == mj:
                    module_coupling_table[mi][mi][1] += 1
                else:
                    module_coupling_table[mi][mj][1] += 1
                    module_coupling_table[mj][mi][1] += 1
            elif couplings[i][j] == 3:  # i to j
                module_coupling_table[mi][mj][2] += 1
                module_coupling_table[mj][mi][3] += 1
            elif couplings[i][j] == 4:  # j to i
                module_coupling_table[mi][mj][3] += 1
                module_coupling_table[mj][mi][2] += 1
            elif couplings[i][j] == 5:  # uncoupled
                if mi == mj:
                    module_coupling_table[mi][mj][4] += 1
                else:
                    module_coupling_table[mi][mj][4] += 1
                    module_coupling_table[mj][mi][4] += 1

    return module_coupling_table, module_order_in_header


def OLD_compute_simple_aggregate_inter_intra_cc_ration(module_coupling_table, module_order_in_header):
    # res will be [internal_coupling_count_list, external_coupling_count_list] where a list is '# fully, partially, i-to-j, j-to-i'
    res = {}
    for mi, mcr in enumerate(module_coupling_table):
        mcv = [module_coupling_table[mi][mi], [0, 0, 0, 0]]
        for mj, mce in enumerate(mcr):
            if mi == mj:
                continue
            mcv[1][0] += mce[0]
            mcv[1][1] += mce[1]
            mcv[1][2] += mce[2]
            mcv[1][3] += mce[3]
        res[module_order_in_header[mi]] = mcv

    sum_intra = 0
    sum_inter = 0
    for mod_cpls in res.itervalues():
        sum_intra += sum(mod_cpls[0])
        sum_inter += sum(mod_cpls[1])

    # print '\n'.join(['%s       %s' % (str(z[0]), str(z[1])) for z in res.itervalues()])
    if sum_inter == 0:
        sum_inter = 0.5
    return sum_intra * 1.0 / sum_inter


def compute_coupling_uncoupling_score(module_coupling_table, module_order_in_header):
    res = {}
    for mi, mcr in enumerate(module_coupling_table):
        intra_cpl_piece = module_coupling_table[mi][mi][0:4]
        intra_uncpl_piece = module_coupling_table[mi][mi][4]

        inter_cpl_piece = [0, 0, 0, 0]
        inter_uncpl_piece = 0
        for mj, mce in enumerate(mcr):
            if mi == mj:
                continue
            inter_cpl_piece[0] += mce[0]
            inter_cpl_piece[1] += mce[1]
            inter_cpl_piece[2] += mce[2]
            inter_cpl_piece[3] += mce[3]

            inter_uncpl_piece += mce[4]

        # XYZ_cpl_piece is 'fully, partially, i-to-j, j-to-i'
        res[module_order_in_header[mi]] = [intra_cpl_piece, intra_uncpl_piece, inter_cpl_piece, inter_uncpl_piece]

    sum_nominator = 0
    sum_denominator = 0
    for mod_name, mod_cpls in res.iteritems():
        if mod_name == eval_constants.EXTERNAL_MODULE_ID:
            pass
        else:
            sum_nominator += sum(mod_cpls[0]) + mod_cpls[3]
            sum_denominator += mod_cpls[1]
            # sum_nominator += sum(mod_cpls[0]) + mod_cpls[3] - mod_cpls[1]
            # sum_denominator += sum(mod_cpls[2])

    # print '\n'.join(['%s       %s' % (str(z[0]), str(z[1])) for z in res.itervalues()])
    # if sum_denominator == 0:
    #     sum_denominator = 0.5
    # return sum_nominator * 1.0 / sum_denominator
    return sum_nominator - sum_denominator

def compute_cohesion_coupling(inf, is_mmod, is_partial_module, S, mets, reacts, revs):
    couplings, blocks = compute_couplings(S, mets, reacts, revs, max(map(lambda s: len(s), reacts)) + 1, max(map(lambda s: len(s), mets)) + 1)

    if is_mmod:
        mmods = eval_util.read_mmods(inf)
        mods = eval_util.convert_metabmod_to_reactmod(mmods, is_partial_module, S, mets, reacts)
    else:
        mods = eval_util.read_rmods(inf)

        # TODO:PARTIAL: all reactions that are not counted in modules (NO: are assigned to external module
        # and) -1 (like blocked reactions) is placed for couplings (in associated rows and columns)
        if is_partial_module:
            done_reacts = set()
            for mod_n, mod_rs in mods.iteritems():
                done_reacts.update(mod_rs)

            not_done_reacts = set(reacts) - set(done_reacts)
            if len(not_done_reacts) > 0:
                mods[eval_constants.EXTERNAL_MODULE_ID] = not_done_reacts

            couplings = copy.deepcopy(couplings)
            for r in not_done_reacts:
                r_idx = reacts.index(r)
                for i in range(len(reacts)):
                    couplings[i][r_idx] = -1
                    couplings[r_idx][i] = -1

    module_coupling_table, module_order_in_header = compute_module_coupling_table(couplings, blocks, mods, reacts)

    return compute_coupling_uncoupling_score(module_coupling_table, module_order_in_header)


if __name__ == '__main__':
    # species = 'toy_model'
    species = 'ecoli_core'

    out_dir = r'%s/%s/newman' % (my_constants.resultPath, species)

    # S, mets, rxns, revs, met_names, rxn_names, biomass, met_comparts = importer.sbmlStoichiometricMatrix(r'D:\Programs\Python\sfw\fca\test\files\networks\pp108.sbml', True, readSpeciesCompart=True)
    S, mets, rxns, revs, met_names, rxn_names, biomass, met_comparts = importer.sbmlStoichiometricMatrix(r'%s/dataset/networks/%s' % (my_constants.basePath, my_constants.species_sbml[species]), True, read_species_compart=True)

    print compute_cohesion_coupling('%s/final_modules.txt' % out_dir, True, S, mets, rxns, revs)