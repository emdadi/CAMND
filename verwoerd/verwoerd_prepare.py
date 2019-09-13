import os

from util import my_constants
from util import my_util
from util import importer
from util.my_util import write_line


def methods_list_getter(species):
    method_dir = r'%s/method/verwoerd' % my_constants.basePath
    cts = eval(open('%s/thresholds.txt' % method_dir, 'r').read())
    return ['verwoerd%s' % l for l in cts[species]]


def read_result_file(resf_path):
    resf = open(resf_path, 'r')

    parts = [[], [], [], []]
    part_idx = -1
    for l in resf:
        l = l.strip()

        if l.startswith('-C'):
            break
        elif l.startswith('-E') or l.startswith('-M'):
            part_idx += 1
        elif part_idx >= 0:
            parts[part_idx].extend([z for z in l.split(' ') if len(z) > 0])

    rmod = parts[0] + parts[1]
    mmod = parts[2]

    resf.close()

    return rmod, mmod


def run_for_all_template_files(template, func):
    i = 1
    while os.path.isfile(template % i):
        try:
            func(template % i, i)
        except TypeError:
            func(template % i)
        i += 1


def go(species, pseudo_species=None):
    method_dir = r'%s/method/verwoerd' % my_constants.basePath
    out_dir = r'%s/%s/verwoerd' % (my_constants.resultPath, species)
    my_util.mkdir_p(out_dir)

    source_file = '%s/dataset/networks/%s' % (my_constants.basePath, my_constants.species_sbml[species])

    S, mets, rxns, revs, met_names, rxn_names, biomass, met_comparts = importer.sbmlStoichiometricMatrix(source_file, True, read_species_compart=True, remove_biomass=True, normalize_stoich=True)

    external_strategy = eval(open('%s/external_strategy.txt' % method_dir, 'r').read())

    # both may either be with/without boundary reactions which is specified by the strategy
    with_boundary_reacts = {}
    with_boundary_mets = list(mets)
    with_boundary_comparts = list(met_comparts)

    next_boundary_met_idx = len(mets)
    if external_strategy[species] == 'KEGG_PSEUDO_NETWORK':     # this implies networks are built from kegg by merging pathways
        subs_to_reacts = {}
        prods_to_reacts = {}
        for j, r in enumerate(rxns):
            subs = []
            prods = []
            for i in range(len(S)):
                if S[i][j] > 0:
                    prods.append([S[i][j], i])    # [stoichiometry, met_name]
                elif S[i][j] < 0:
                    subs.append([-S[i][j], i])    # [stoichiometry, met_name]
            with_boundary_reacts[r] = [subs, prods]

            has_subs, has_prods = False, False
            for i in range(len(S)):
                if S[i][j] > 0 or (S[i][j] < 0 and revs[j] != 0):
                    if mets[i] not in prods_to_reacts:
                        prods_to_reacts[mets[i]] = set()
                    prods_to_reacts[mets[i]].add(r)
                    has_prods = True
                if S[i][j] < 0 or (S[i][j] > 0 and revs[j] != 0):
                    if mets[i] not in subs_to_reacts:
                        subs_to_reacts[mets[i]] = set()
                    subs_to_reacts[mets[i]].add(r)
                    has_subs = True

            if not has_subs or not has_prods:
                raise Exception('reaction without substrate/products in KEGG pseudo networks! should never happen!')

        # checks whether the reaction is internal in presence of reversiblity
        def is_internal(met):
            if met not in prods_to_reacts or met not in subs_to_reacts:
                return False
            if prods_to_reacts[met] == subs_to_reacts[met] and len(prods_to_reacts[met]) == 1:
                the_react = list(prods_to_reacts[met])[0]
                the_react_subs = [mets[ss[1]] for ss in with_boundary_reacts[the_react][0]]
                the_react_prods = [mets[ss[1]] for ss in with_boundary_reacts[the_react][1]]
                if not (met in the_react_subs and met in the_react_prods):  # not polymeric
                    return False
            return True

        initial_external_metabolites_idxs = [mets.index(m) for m in with_boundary_mets if not is_internal(m)]
        # internal_metabolites = [m for m in with_boundary_mets if m in mets_used_as_subs and m in mets_used_as_prods]
    else:
        # may add to with_boundary_mets based on strategy.
        # may update with_boundary_reacts according to newly added with_boundary_mets
        for j, r in enumerate(rxns):
            subs = []
            prods = []
            for i in range(len(S)):
                if S[i][j] > 0:
                    prods.append([S[i][j], i])
                elif S[i][j] < 0:
                    subs.append([-S[i][j], i])

            if not subs:
                raise Exception('reaction without substrate! should never happen!')
            elif not prods:
                if external_strategy[species] == 'BOUNDARY':
                    if len(subs) != 1:
                        raise Exception('EXCHANGE reaction with more than one substrate metabolite! unacceptable.')
                    new_boundary_met = mets[subs[0][1]] + '_BND'
                    if new_boundary_met in with_boundary_mets:
                        # print 'should happen only once for homo recon 1 and once for ecoli ijo 1366!'
                        prods = [[subs[0][0], with_boundary_mets.index(new_boundary_met)]]
                    else:
                        with_boundary_mets.append(new_boundary_met)
                        with_boundary_comparts.append('B')
                        prods = [[subs[0][0], next_boundary_met_idx]]
                        next_boundary_met_idx += 1
                else:
                    print 'skipping reaction without product: EXCHANGE!'
                    continue

            with_boundary_reacts[r] = [subs, prods]

        # with_boundary_S = [[] for i in range(len(with_boundary_mets))]
        # for ri, r in rxns:
        #     r_def = with_boundary_reacts[r]
        #     for e in r_def[0]:  # substrates
        #         if e[0] != 0:

        if external_strategy[species] == 'BOUNDARY':
            initial_external_metabolites_idxs = [i for i, m in enumerate(with_boundary_mets) if with_boundary_mets[i].endswith('_BND')]
        else:
            initial_external_metabolites_idxs = [i for i, m in enumerate(with_boundary_mets) if with_boundary_mets[i].endswith('_e')]
            # internal_metabolites = [m for m in with_boundary_mets if m not in initial_external_metabolites]

    if pseudo_species:
        my_util.mkdir_p(pseudo_species)
        inf = open('%s/verwoerd.tsv' % pseudo_species, 'w')
    else:
        inf = open('verwoerd_%s.tsv' % species, 'w')
    write_line(inf, '<Title>')
    write_line(inf, species.replace('/', ''))
    write_line(inf, '<Reactions>')
    write_line(inf, '\t'.join(rxns))
    write_line(inf, '<ReversibleReactions>')
    write_line(inf, '\t'.join([r for i, r in enumerate(rxns) if revs[i] == 1]))
    write_line(inf, '<InternalCompounds>')
    # write_line(inf, '\t'.join(['%s' % (with_boundary_mets[i]) for i in range(len(with_boundary_mets))]))
    write_line(inf, '\t'.join(['%s %s' % (with_boundary_mets[i], with_boundary_comparts[i]) for i in range(len(with_boundary_mets))]))
    write_line(inf, '<ExternalCompounds>')
    # write_line(inf, '\t'.join(['%s %s' % (with_boundary_mets[i], with_boundary_comparts[i]) for i in initial_external_metabolites_idxs]))
    write_line(inf, '<Stoichiometry>')
    write_line(inf, '%%MatrixMarket matrix coordinate real general')

    s_lines = []
    nnz_s = 0
    for ri, r in enumerate(rxns):
        if r not in with_boundary_reacts:   # exchange reactions
            continue

        rdef = with_boundary_reacts[r]
        for e in rdef[0]:
            s_lines.append('%d\t%d\t%s' % (e[1] + 1, ri + 1, -e[0]))
            nnz_s += 1
        for e in rdef[1]:
            s_lines.append('%d\t%d\t%s' % (e[1] + 1, ri + 1, e[0]))
            nnz_s += 1

    write_line(inf, '\t\t\t%d\t%d\t%d' % (len(with_boundary_mets), len(with_boundary_reacts), nnz_s))
    for sl in s_lines:
        write_line(inf, sl)

    inf.close()

    if pseudo_species:
        inf = open('%s/verwoerd-ext.txt' % pseudo_species, 'w')
    else:
        inf = open('verwoerd-ext_%s.txt' % species, 'w')
    for init_extern_idx in initial_external_metabolites_idxs:
        write_line(inf, with_boundary_mets[init_extern_idx])
    inf.close()


if __name__ == '__main__':
    for species in my_constants.species_sbml.keys():
        go(species)
