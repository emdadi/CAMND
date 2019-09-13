import shutil

from util import my_constants
from util import my_util
from util import importer
from util.my_util import write_line


def go(species):
    method_dir = r'%s/method/muller' % my_constants.basePath
    out_dir = r'%s/%s/muller' % (my_constants.resultPath, species)
    my_util.mkdir_p(out_dir)

    source_file = '%s/dataset/networks/%s' % (my_constants.basePath, my_constants.species_sbml[species])

    S, mets, rxns, revs, met_names, rxn_names, biomass, met_comparts = importer.sbmlStoichiometricMatrix(source_file, True, read_species_compart=True, remove_biomass=False, normalize_stoich=False)

    f = open('muller.m', 'w')
    write_line(f, 'addpath %s/code' % method_dir)
    write_line(f, "model = readCbModel('%s')" % source_file)
    write_line(f, "model.c(%d) = 1;" % (rxns.index(biomass[0]) + 1))
    write_line(f, "changeCobraSolver('glpk');")
    write_line(f, '[modules, var, flux] = computeModulesOpt( model );')
    write_line(f, 'save mullerout.mat modules var flux')

    f.close()

    my_util.prepare_matlab_file_and_exec_and_wait_finish('muller', 'mullerout.mat', False)

    res_vars = my_util.try_load_matlab_result('mullerout.mat')
    raw_modules = res_vars['modules']  # if matlab has failed, this will throw exception!
    shutil.copy('mullerout.mat', out_dir)

    raw_modules = raw_modules.T.tolist()  # eacho row will be a module where reactions are marked
    modules = []
    for raw_module in raw_modules:
        modules.append([])
        for rIdx, in_module in enumerate(raw_module):
            if in_module == 1:
                modules[-1].append(rxns[rIdx])

    out = open('%s/final_modules.txt' % out_dir, 'w')
    out.write("#each row is a module of reactions. not all reactions are specified (nature of this method only select some reaction to be modules)\n")
    for m in modules:
        out.write(' '.join(m))
        out.write('\n')
    out.close()

if __name__ == '__main__':
     # go('ecoli_core')
    # go('helico_iit341')
    ## go('ecoli_ijr904')
    # go('ecoli_iaf1260')
    # go('ecoli_ijo1366')
    # go('saccaro_ind750')
    go( 'mus_imm1415')
    #   go('arabidopsis_irs1597')
    # dist = np.array([
    # [0,1,100,100],
    #     [1,0,100,100],
    #     [100,100,0,1],
    #     [100,100,1,0]
    # ], dtype=float)
    # print wpgma(dist)