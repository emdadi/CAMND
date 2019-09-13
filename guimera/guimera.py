from util import my_constants
from util import importer
from util import my_util
import os
import shutil

def go(species):
    out_dir = r'%s/%s/guimera' % (my_constants.resultPath, species)
    my_util.mkdir_p(out_dir)

    S, mets, rxns, revs, met_names, rxn_names, biomass, met_comparts = importer.sbmlStoichiometricMatrix(r'%s/dataset/networks/%s' % (my_constants.basePath, my_constants.species_sbml[species]), True, read_species_compart=True, remove_biomass=True, normalize_stoich=True)

    # graph is represented as list of edges between metabs
    grph = my_util.graph_by_explode_reactions_to_complete_bipartite(S, mets)

    grphIdx = [(mets.index(m1), mets.index(m2)) for (m1, m2) in grph]
    grphIdx.sort()
    inf = open('guimera.in', 'w')
    for eIdx in grphIdx:
        inf.write('%d %d\n' % (eIdx[0], eIdx[1]))
    inf.close()

    # netcarto_cl net_file_name seed T_ini iteration_factor cooling_factor
    # T_ini, iteration_factor, and cooling_factor can be set to -1 to use the defaults (2/size_of_network, 1.0, and 0.995, respectively).
    res = os.system('netcarto_cl guimera.in %d -1 -1 -1 0')
    shutil.copy('modules.dat', out_dir)
    shutil.copy('roles.dat', out_dir)

    out = open('%s/final_modules.txt' % out_dir, 'w')
    out.write('# each line one module!\n')
    netcarto_outf = open('modules.dat', 'r')
    for l in netcarto_outf:
        didx = l.find('---')
        if didx != -1:
            l = l[didx + len('---'):].strip()
            mIdxs = l.split(' ')
            out.write(' '.join([mets[int(i)] for i in mIdxs]))
            out.write('\n')
    netcarto_outf.close()
    out.close()

if __name__ == '__main__':
    go('ecoli_core')
    #go('toy_model')
