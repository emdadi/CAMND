from util import my_constants
from util import importer
from util import my_util
import shutil
import correctedcmty

def go(species):
    out_dir = r'%s/%s/newman' % (my_constants.resultPath, species)
    my_util.mkdir_p(out_dir)

    S, mets, rxns, revs, met_names, rxn_names, biomass, met_comparts = importer.sbmlStoichiometricMatrix(r'%s/dataset/networks/%s' % (my_constants.basePath, my_constants.species_sbml[species]), True, read_species_compart=True, remove_biomass=True)

    # graph is represented as list of edges between metabs
    grph = my_util.graph_by_explode_reactions_to_complete_bipartite(S, mets)

    grphIdx = [(mets.index(m1), mets.index(m2)) for (m1, m2) in grph]
    grphIdx.sort()
    inf = open('newman.in', 'w')
    for eIdx in grphIdx:
        inf.write('%d,%d,1\n' % (eIdx[0], eIdx[1]))
    inf.close()

    correctedcmty.main(['DUMMY', 'newman.in', 'newman.out'])
    shutil.copy('newman.out', out_dir)

    resf = open('newman.out', 'r').read()
    stidx = resf.rindex('START_COMP')
    edidx = resf.rindex('END_COMP')
    res = resf[stidx + len('START_COMP'): edidx]
    res_lines = res.split('\n')

    out = open('%s/final_modules.txt' % out_dir, 'w')
    out.write('# each line one module!\n')
    for l in res_lines:
        l = l.strip()
        if l == '':
            continue
        mIdxs = eval(l)
        out.write(' '.join([mets[int(i)] for i in mIdxs]))
        out.write('\n')
    out.close()

if __name__ == '__main__':
    go('ecoli_core')
    # go('helico_iit341')
    # go('ecoli_ijr904')
    # go('ecoli_iaf1260')
    # go('ecoli_ijo1366')
    # go('saccaro_ind750')
    # go('homo_recon1')
    # go('toy_model')
