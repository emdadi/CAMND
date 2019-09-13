import shutil
import os

from util import my_constants
from util import my_util
from util import importer
from util.my_util import write_line


def read_result_file(resf_path):
    resf = open(resf_path, 'r')

    parts = [[], [], [], [], []]
    part_idx = -1
    for l in resf:
        l = l.strip()

        if l.startswith('<Stoich'):
            break
        elif l.startswith('<Titl') or l.startswith('<Rea') or l.startswith('<Rev') or l.startswith('<Int') or l.startswith('<Ext'):
            part_idx += 1
        elif part_idx >= 0:
            parts[part_idx].extend([z.split(' ')[0] for z in l.split('\t') if len(z) > 0])  # metabolite names are 'metabname compart' when outputted from verwoerd!

    rmod = parts[1] + parts[2]
    mmod = parts[3]

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

    # for each species, this may be defined or simply be '' which is ignored
    thresholds = eval(open('%s/thresholds.txt' % method_dir, 'r').read())

    for thr in thresholds[species]:
        if pseudo_species:
            manual_results_dir = '%s/Subnetworks_%s' % (pseudo_species, thr)
        else:
            manual_results_dir = '%s_%s' % (species, thr)

        rmods = {}
        mmods = {}

        def read_module_and_move(resf_path, i):
            rmod, mmod = read_result_file(resf_path)
            rmods['%d' % i] = rmod
            mmods['%d' % i] = mmod

            int_out_dir_path = '%s/subsystem_%s_%d.out' % (out_dir, thr, i)
            shutil.copy(resf_path, int_out_dir_path)

        if pseudo_species:
            run_for_all_template_files(manual_results_dir + '/' + species.replace('/', '') + '_Block_%d.tsv', read_module_and_move)
        else:
            run_for_all_template_files(manual_results_dir + '/' + species + '_Block_%d.tsv', read_module_and_move)

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
    go('helico_iit341')
    go('saccaro_ind750')
    go('ecoli_iaf1260')
    go('ecoli_ijo1366')
    go('homo_recon1')
    go('mbarkeri_iaf692')
    go('arabidopsis_irs1597')
    go('mus_imm1415')
