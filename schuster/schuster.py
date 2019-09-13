import shutil
import os

from util import my_constants
from util import my_util
from util import importer
from util.my_util import write_line


def methods_list_getter(species):
    method_dir = r'%s/method/schuster' % my_constants.basePath
    cts = eval(open('%s/thresholds.txt' % method_dir, 'r').read())
    return ['schuster%s' % l for l in cts[species]]


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


def go(species):
    method_dir = r'%s/method/schuster' % my_constants.basePath
    out_dir = r'%s/%s/schuster' % (my_constants.resultPath, species)
    my_util.mkdir_p(out_dir)

    source_file = '%s/dataset/networks/%s' % (my_constants.basePath, my_constants.species_sbml[species])

    S, mets, rxns, revs, met_names, rxn_names, biomass, met_comparts = importer.sbmlStoichiometricMatrix(source_file, True, read_species_compart=True, remove_biomass=True, normalize_stoich=True)

    external_strategy = eval(open('%s/external_strategy.txt' % method_dir, 'r').read())

    # both may either be with/without boundary reactions which is specified by the strategy
    # these will be used to create schuster.in file which is used by the method as the input network
    with_boundary_reacts = {}
    with_boundary_mets = list(mets)

    if external_strategy[species] == 'KEGG_PSEUDO_NETWORK':     # this implies networks are built from kegg by merging pathways
        subs_to_reacts = {}
        prods_to_reacts = {}
        for j, r in enumerate(rxns):
            subs = []
            prods = []
            for i in range(len(S)):
                if S[i][j] > 0:
                    prods.append([S[i][j], mets[i]])    # [stoichiometry, met_name]
                elif S[i][j] < 0:
                    subs.append([-S[i][j], mets[i]])    # [stoichiometry, met_name]
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
                the_react_subs = [ss[1] for ss in with_boundary_reacts[the_react][0]]
                the_react_prods = [ss[1] for ss in with_boundary_reacts[the_react][1]]
                if not (met in the_react_subs and met in the_react_prods):  # not polymeric
                    return False
            return True

        internal_metabolites = [m for m in with_boundary_mets if is_internal(m)]
        initial_external_metabolites = [m for m in with_boundary_mets if m not in internal_metabolites]
    else:
        # may add to with_boundary_mets based on strategy.
        # may update with_boundary_reacts according to newly added with_boundary_mets
        for j, r in enumerate(rxns):
            subs = []
            prods = []
            for i in range(len(S)):
                if S[i][j] > 0:
                    prods.append([S[i][j], mets[i]])    # [stoichiometry, met_name]
                elif S[i][j] < 0:
                    subs.append([-S[i][j], mets[i]])    # [stoichiometry, met_name]

            # this strategy means reactions without product are boundary reactions
            # note that: reversibility could potentially be problematic BUT I HAVE CHECKED BIGG NETWORKS, THESE KIND OF REACTIONS ARE REAL OUTSIDE BOUNDARY
            if not subs:
                raise Exception('reaction without substrate! should never happen!')
            elif not prods:
                if external_strategy[species] == 'BOUNDARY':
                    if len(subs) != 1:
                        raise Exception('EXCHANGE reaction with more than one substrate metabolite! unacceptable.')
                    new_boundary_met = subs[0][1] + '_BND'
                    if new_boundary_met not in with_boundary_mets:
                        # print 'the only exceptions are once for homo recon 1 and once for ecoli ijo 1366!'
                        with_boundary_mets.append(new_boundary_met)
                    prods = [[subs[0][0], new_boundary_met]]
                else:
                    print 'skipping reaction without product: EXCHANGE!'
                    continue

            with_boundary_reacts[r] = [subs, prods]

        if external_strategy[species] == 'BOUNDARY':
            initial_external_metabolites = [m for i, m in enumerate(with_boundary_mets) if with_boundary_mets[i].endswith('_BND')]
            internal_metabolites = [m for m in with_boundary_mets if m not in initial_external_metabolites]
        else:
            initial_external_metabolites = [m for i, m in enumerate(with_boundary_mets) if with_boundary_mets[i].endswith('_e')]
            internal_metabolites = [m for m in with_boundary_mets if m not in initial_external_metabolites]

    inf = open('schuster.in', 'w')
    write_line(inf, '-ENZREV')
    write_line(inf, ' '.join([r for i, r in enumerate(rxns) if revs[i] == 1]))
    write_line(inf, '')
    write_line(inf, '-ENZIRREV')
    write_line(inf, ' '.join([r for i, r in enumerate(rxns) if revs[i] != 1]))
    write_line(inf, '')
    write_line(inf, '-METINT')
    write_line(inf, ' '.join(internal_metabolites))
    write_line(inf, '')
    write_line(inf, '-METEXT')
    write_line(inf, ' '.join(initial_external_metabolites))
    write_line(inf, '')
    write_line(inf, '-CAT')
    for r_name, r in with_boundary_reacts.iteritems():
        react_str = '%s : %s = %s .' % (r_name, ' + '.join(['%d %s' % (p[0], p[1]) for p in r[0]]), ' + '.join(['%d %s' % (p[0], p[1]) for p in r[1]]))
        write_line(inf, react_str)
    inf.close()

    # for each species, this may be defined or simply be '' which is ignored
    thresholds = eval(open('%s/thresholds.txt' % method_dir, 'r').read())

    for thr in thresholds[species]:
        run_for_all_template_files('subsystem%d.out', os.remove)
        run_for_all_template_files(out_dir + '/subsystem%d.out', os.remove)

        # subnet file_name
        if my_constants.win:
            res = os.system('%s/src/subnet.exe schuster.in %s' % (method_dir, thr))
        else:
            res = os.system('subnet schuster.in %s' % thr)

        rmods = {}
        mmods = {}

        def read_module_and_move(resf_path, i):
            rmod, mmod = read_result_file(resf_path)
            rmods['%d' % i] = rmod
            mmods['%d' % i] = mmod

            int_out_dir_path = '%s/subsystem_%s_%d.out' % (out_dir, thr, i)
            shutil.copy(resf_path, int_out_dir_path)

        run_for_all_template_files('subsystem%d.out', read_module_and_move)

        outr = open('%s/react_modules_%s.txt' % (out_dir, thr), 'w')
        outrm = open('%s/metab_react_modules_%s.txt' % (out_dir, thr), 'w')
        outm = open('%s/metab_modules_%s.txt' % (out_dir, thr), 'w')
        for mname, rmod in rmods.iteritems():
            write_line(outr, ' '.join(rmod))
            write_line(outrm, ' '.join(mmods[mname]))
            if mmods[mname]:
                write_line(outm, ' '.join(mmods[mname]))
        outr.close()
        outrm.close()
        outm.close()

        shutil.copy('%s/metab_modules_%s.txt' % (out_dir, thr), '%s/final_modules_%s.txt' % (out_dir, thr))  # TODO: which file r/m should be selected as final modules?


if __name__ == '__main__':
    go('ecoli_core')
    # go('ecoli_iaf1260')
