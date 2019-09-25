import util.importer
from util.my_util import good_split


def read_fcdecomp(folders):
    parts = {}

    for folder in folders:
        f = open(folder + '/results.txt', 'r')
        for l in f:
            partstr = l[l.index('{'): l.index('}') + 1]
            a2part = eval(partstr)
            break
        f.close()

        sl_idx = folder.rfind('/')
        bsl_idx = folder.rfind('\\')
        if sl_idx > bsl_idx:
            base_name = folder[sl_idx + 1:]
        else:
            base_name = folder[bsl_idx + 1:]

        for n, p in a2part.iteritems():
            parts[base_name + '_' + n] = p

    return parts


def read_compart_by_label(sbml_path, to_remove_metabs):
    orig_s, orig_mets, orig_rxns, orig_revs, orig_met_names, orig_comparts = netpart.importer.sbmlStoichiometricMatrix(
        sbml_path, readSpeciesCompart=True)
    res = {}
    other = []
    for mi, met in enumerate(orig_mets):
        if met in to_remove_metabs:
            continue

        cidx = met.rfind('_')
        if cidx < 0:
            other.append(met)
        else:
            compart = met[cidx:]
            if compart not in res:
                res[compart] = []
            res[compart].append(met)
    if len(other) > 0:
        res['OTHER'] = other
    return res


def read_compart(sbml_path, ext_compart):
    orig_s, orig_mets, orig_rxns, orig_revs, orig_met_names, orig_comparts = netpart.importer.sbmlStoichiometricMatrix(
        sbml_path, readSpeciesCompart=True)
    res = {}
    ext = []
    for mi, met in enumerate(orig_mets):
        compart = orig_comparts[mi]
        if compart == ext_compart:
            ext.append(met)
        else:
            if compart not in res:
                res[compart] = list()
            res[compart].append(met)
    return res, ext


def read_subsys(sbml_path, to_change_subsys_map, null_subsys_name, reduced_metabs=None):
    orig_s, orig_mets, orig_rxns, orig_revs, orig_met_names, orig_comparts = netpart.importer.sbmlStoichiometricMatrix(
        sbml_path, readSpeciesCompart=True)

    orig_rom = dict((i, filter(lambda j: orig_s[i][j] != 0, range(len(orig_rxns)))) for i in range(len(orig_mets)))

    f = open(sbml_path, 'r')
    react = False
    react_to_subsys = {}
    for l in f:
        l = l.strip()
        if l.startswith('<reaction'):
            if react:
                react_to_subsys[react_name] = null_subsys_name
                react = False
            sti = l.index('id="')+len('id="')
            react_name = l[sti: l.index('"', sti)]
            react = True

        if l.find('SUBSYSTEM:') >= 0:
            sti = 'SUBSYSTEM:'
            subs = l[l.find(sti) + len(sti): l.find('</html', 1)].strip()
            if subs in to_change_subsys_map:
                subs = to_change_subsys_map[subs]
            react_to_subsys[react_name] = subs
            react = False
    f.close()

    res = {}

    for mi, ris in orig_rom.iteritems():
        m = orig_mets[mi]
        if reduced_metabs is not None and m not in reduced_metabs:
            continue
        rs = [orig_rxns[ri] for ri in ris]
        for r in rs:
            s = react_to_subsys[r]
            if s not in res:
                res[s] = []
            res[s].append(m)

    return res


def read_cmty(parts_file, metabs_file):
    f = open(metabs_file, 'r')
    for l in f:
        mets = eval(l.strip())
        break
    f.close()

    res = {}
    exts_and_virts = {}

    f = open(parts_file, 'r')
    i = 1
    for l in f:
        part_idxs = eval(l.strip())
        res[str(i)] = [mets[mi] for mi in part_idxs if mi < len(mets)]
        exts_and_virts[str(i)] = [mi for mi in part_idxs if mi >= len(mets)]
        i += 1
    f.close()

    return res, exts_and_virts


def read_fca2(parts_file):
    res = {}
    exts_and_virts = {}

    f = open(parts_file, 'r')
    i = 1
    for l in f:
        res[str(i)] = eval(l.strip())
        i += 1
    f.close()

    return res, exts_and_virts


def read_hepato(parts_file):
    parts = {}

    f = open(parts_file, 'r')
    i = 0
    for l in f:
        l = l.strip()
        ps = good_split(l, ',')

        if ps[0] == "'MY ID'":
            if i > 3:   # skip for first module
                parts[pid] = [ppid, prxns]
            pid = ps[1]
            prxns = []
        elif ps[0] == "'Parental ID'":
            ppid = ps[1]
        elif ps[0].startswith("'Reaction "):
            prxns.append('R%s' % (ps[0][1:-1].split(' ')[1]))

        i += 1

    res = {}
    exts_and_virts = {}

    return res, exts_and_virts

# z = read_hepato(r'C:\Users\Abolfazl\Desktop\thesis\results\journal.pcbi.1002262.s001.csv')


def read_method_metabolite_modules(species, method, conversion=None):
    # TODO: each method will have its own function and this will be a simple switch-case (when converison=None).
    # TODO: method function should take a folder as input (all special files for the method will be there)
    pass


def read_method_reaction_modules(species, method, conversion=None):
    pass


def read_method_metabolite_hierarchy(species, method, conversion=None):
    pass


def read_method_reaction_hierarchy(species, method, conversion=None):
    pass


def read_target_metabolite_modules(species, method):
    pass


def read_target_reaction_modules(species, method):
    pass