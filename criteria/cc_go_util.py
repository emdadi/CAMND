from copy import deepcopy
import math

import eval_util
import eval_constants


# enz_grp: a list of enzyme_group where each entry is an iso-enzyme. structure of an
# enzyme_group is [group_tag, sub_enzyme_group1 [, ...]] where group_tag is either
# OR/AND and each sub_enzyme_group is of type enzyme_group
def accumulate_all_enzymes_altogether(enz_grp):
    if type(enz_grp) != list:
        return [enz_grp]

    res = []
    for i in range(1, len(enz_grp)):
        if type(enz_grp[i]) == list:
            res.extend(accumulate_all_enzymes_altogether(enz_grp[i]))
        else:
            res.append(enz_grp[i])

    return res


def recursive_and_expression_compare(a, b, simtbl):
    vec = sorted(list(set(a + b)))

    aval, bval = [], []
    for v in vec:
        if v in a:
            va = 1.0
        else:
            va = max([recursive_orand_expression_compare(v, v2, simtbl) for v2 in a])

        if v in b:
            vb = 1.0
        else:
            vb = max([recursive_orand_expression_compare(v, v2, simtbl) for v2 in b])

        aval.append(va)
        bval.append(vb)

    anorm = pow(sum([pow(av, 2) for av in aval]), 0.5)
    bnorm = pow(sum([pow(bv, 2) for bv in bval]), 0.5)

    # print aval, bval
    x = sum([av * bv for av, bv in zip(aval, bval)]) * 1.0 / (anorm * bnorm)
    if abs(x) > 1:
        if abs(x) - 1 > 0.00001:
            print 'recursive_and_expression_compare resulted in math.acos(%.10f)! exiting...' % x
            exit(1)
        else:
            x = 1.0 if x > 0 else -1.0
    z = math.acos(x)
    return 1 - z * 2.0 / math.pi


def recursive_or_expression_compare(a, b, simtbl):
    sop = 0
    for av in a:
        for bv in b:
            sop += recursive_orand_expression_compare(av, bv, simtbl)

    x = sop * 1.0 / (len(a) * len(b))
    return x


def recursive_orand_expression_compare(grp1, grp2, simtbl):
    oper1 = grp1[0] if isinstance(grp1, tuple) else ''
    oper2 = grp2[0] if isinstance(grp2, tuple) else ''

    if oper1 == '' and oper2 == '':
        if grp1 in simtbl and grp2 in simtbl:
            return simtbl[grp1][grp2]
        else:
            print 'MISSED WHEN COMPUTING SIMILARITY OF AND/OR, SHOULD NOT HAPPEN, EXITING ... %s ... %s' % (str(grp1), str(grp2))
            exit(1)

    if oper1 != oper2:
        if oper1 != '':
            grp2 = (oper1, grp2)
        else:   # if oper2 != '':
            grp1 = (oper2, grp1)

    if grp1[0] == 'and':
        return recursive_and_expression_compare(grp1[1:], grp2[1:], simtbl)
    else:   # if grp1[0] == 'or':
        return recursive_or_expression_compare(grp1[1:], grp2[1:], simtbl)


def recursive_remove_unmapped_genes_and_tupalize(enz_grp, simtbl):
    if not isinstance(enz_grp, list):
        if enz_grp not in simtbl:
            return None
        else:
            return enz_grp

    if enz_grp[0] != 'and' and enz_grp[0] != 'or':
        print 'ENZYME GROUP WITHOUT STARTING and/or, SHOULD NEVER HAPPEN!'
        exit(1)

    i = 1
    while i < len(enz_grp):
        enz_grp[i] = recursive_remove_unmapped_genes_and_tupalize(enz_grp[i], simtbl)
        if enz_grp[i] is None:
            enz_grp.pop(i)
        else:
            i += 1

    if len(enz_grp) == 1:   # all but initial 'and'/'or' are removed!
        return None
    elif len(enz_grp) == 2:     # [and/or elem] is converted to 'elem'
        return enz_grp[1]
    else:
        return tuple(enz_grp)


# z = 0
def get_even_ecnumber_annotation_term_ids_for_reaction(r, species, criteria, go_loaded_model, moredata_loaded_model, stats_mode=False):
    if stats_mode:
        mapping_state, enz_grp1 = eval_util.get_annotation_term_ids_for_reaction(r, species, criteria, go_loaded_model, stats_mode=True)
        if enz_grp1 is not None:
            return mapping_state, enz_grp1
    else:
        enz_grp1 = eval_util.get_annotation_term_ids_for_reaction(r, species, criteria, go_loaded_model, stats_mode=False)
        if enz_grp1 is not None:
            return enz_grp1

    r_moredata = eval_util.get_moredata_for_reaction(r, moredata_loaded_model)

    if stats_mode:
        if r_moredata is None:
            return mapping_state, None
    else:
        if r_moredata is None:
            return None

    ecs = r_moredata['EC Numbers']
    if len(ecs) >= 2:
        mapping_state2 = -5
        res = ['or'] + ecs
    elif len(ecs) == 1:
        mapping_state2 = -5
        res = ecs[0]
    else:
        if stats_mode:
            mapping_state2 = mapping_state
        res = None

    if stats_mode:
        return mapping_state2, res      # originally not enzymatic/undefined gpd/missed but instead we have used ec number!
    else:
        return res


def filter_out_reactions_by_type(mods, filter_rules):
    mods_res = {}
    for mod_name, mod in mods.iteritems():
        mod_res = []
        for m in mod:
            for fr in filter_rules:
                if not ((fr[0] == 0 and m.startswith(fr[1])) or (fr[0] == 1 and m.endswith(fr[1]))):
                    mod_res.append(m)
        if len(mod_res) > 0:
            mods_res[mod_name] = mod_res
    return mods_res


def compute_similarity_of_gene_pair(cached_reaction_similarities, r1, r2, sim_table, species, criteria, factor1, factor2, go_loaded_model, moredata_loaded_model):
    if (r1, r2, species, factor1, factor2) in cached_reaction_similarities:  # and not stats_mode:
        return cached_reaction_similarities[(r1, r2, species, factor1, factor2)]
    else:
        # global z
        # if z % 100000 == 0:
        #     print 'compute_similarity_of_gene_pair %s %s %d %d' % (species, factor1, z, len(cached_reaction_similarities))
        # z+= 1
        enz_grp1 = get_even_ecnumber_annotation_term_ids_for_reaction(r1, species, criteria, go_loaded_model, moredata_loaded_model)
        enz_grp2 = get_even_ecnumber_annotation_term_ids_for_reaction(r2, species, criteria, go_loaded_model, moredata_loaded_model)

        # TODO:MISSING: if no gene association with a reaction is found, similarity is zero!
        # reason1: reaction really is not enzymatic
        # reason2: all geneProducts associated with reaction are improperly defined
        # reason3: genes cannot be mapped to UniprotIDs
        # reason4: mapped UniprotID cannot be found in GOA file (so there is no similarity information)
        if enz_grp1 is None or enz_grp2 is None:
            cached_reaction_similarities[(r1, r2, species, factor1, factor2)] = 0
            return 0

        # remove all missed genes in mapping to GO/CC here, (NOT: and compute count of missed) (NULL lists are removed)
        enz_grp1, enz_grp2 = deepcopy(enz_grp1), deepcopy(enz_grp2)
        enz_grp1 = recursive_remove_unmapped_genes_and_tupalize(enz_grp1, sim_table)
        enz_grp2 = recursive_remove_unmapped_genes_and_tupalize(enz_grp2, sim_table)

        if enz_grp1 is None or enz_grp2 is None:
            cached_reaction_similarities[(r1, r2, species, factor1, factor2)] = 0
            return 0

        sim = recursive_orand_expression_compare(enz_grp1, enz_grp2, sim_table)

        cached_reaction_similarities[(r1, r2, species, factor1, factor2)] = sim

        return sim


def general_compute_distance(cached_reaction_similarities, similarity_table, table_gene_order, inf, species, criteria, factor1, factor2, is_mmod, is_partial_module, S, mets, reacts, revs, go_loaded_model, moredata_loaded_model):
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

    if species == 'arabidopsis_irs1597':
        mods = filter_out_reactions_by_type(mods, [(0, 'R_EX')])
    else:
        mods = filter_out_reactions_by_type(mods, [(0, 'R_EX_')])

    mods_similarity_table = {}
    for mod1_name, mod1 in mods.iteritems():
        mods_similarity_table[mod1_name] = {}
        for mod2_name, mod2 in mods.iteritems():
            sim = 0
            for r1 in mod1:
                for r2 in mod2:
                    if mod1_name != mod2_name or r1 != r2:  # designed to work for overlapping modules. only filters out a gene-similarity-with-itself for intra-module
                        sim += compute_similarity_of_gene_pair(cached_reaction_similarities, r1, r2, similarity_table, species, criteria, factor1, factor2, go_loaded_model, moredata_loaded_model)
            mods_similarity_table[mod1_name][mod2_name] = sim

    # print 'mod simetable done'

    # TODO:BUG: for overlapping modules this fails because (len(rxns) - len(mods[mod_name])) is wrong
    intra_score = 0
    inter_score = 0
    for mod_name, mods_sim in mods_similarity_table.iteritems():
        intra_piece = mods_sim[mod_name] * 1.0 / pow(len(mods[mod_name]), 2)
        inter_piece = sum([mods_sim[other_mod_name] * 1.0 / (len(mods[mod_name]) * len(mods[other_mod_name])) for other_mod_name in mods_similarity_table.iterkeys() if other_mod_name != mod_name])
        if mod_name == eval_constants.EXTERNAL_MODULE_ID:
            intra_score -= intra_piece
            inter_score -= inter_piece
        else:
            intra_score += intra_piece
            inter_score += inter_piece

    # print 'done'

    return intra_score - inter_score

    # if inter_score == 0:
    #     return 0  # this happens in cases where decomposition is not really done and there is only one module!!!
    # else:
    #     return intra_score * 1.0 / inter_score


