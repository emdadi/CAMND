try
    addpath D:/project/comparison/method/muller2_new/code
    model = readCbModel('D:/project/comparison/dataset/networks/iAF692.xml')
    model.c(98) = 1;
    changeCobraSolver('glpk');
    [modules, var, flux] = computeModulesOpt( model );
    save muller2_newout.mat modules var flux
catch
    zzz = lasterror;
    save muller2_newout.mat zzz
end
quit
