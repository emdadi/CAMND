import json
import os
import re

import libsbml

import my_constants
import eval_constants
import  my_util


"""
    metabmod: a map from module name to list of metabolites
    reactmod: a map from module name to list of reactions
"""


class MoreDataLoadedModel:
    def __init__(self, species):
        self.bigg = {}

        fname = my_constants.species_sbml[species].replace('.xml', '_moredata.json')
        inf = open('%s/dataset/networks/%s' % (my_constants.basePath, fname), 'r')
        elements = eval(inf.read().replace(': null', ': ""'))
        for e in elements:
            if e['bigg_id'] in self.bigg:
                print 'creating MoreDataLoadedModel for %s failed: more than one element with bigg_id %s' % (species, e['bigg_id'])
                # exit(1)
            self.bigg[e['bigg_id']] = e
        inf.close()


class GoLoadedModel:
    def gene_product_count(self):
        return self.gp_count

    def undefined_gene_product_count(self):
        return self.undefined_gp_count

    def __init__(self, sbml_model, species):
        self.gp_count = 0
        self.undefined_gp_count = 0
        self.sbml_model = sbml_model

        self.gene_product_definitions = {}
        for raw_gp in sbml_model.getPlugin('fbc').getListOfGeneProducts():
            self.gp_count += 1
            if not raw_gp.getAnnotation():
                self.undefined_gp_count += 1
                print 'WARNING: definition of gene product %s not found' % raw_gp.getId()
                continue

            description_tag = raw_gp.getAnnotation().getChild('RDF').getChild('Description')

            to_iterate_database_ids = []

            if description_tag.hasChild('is'):
                is_tag = description_tag.getChild('is')
                bag = is_tag.getChild('Bag')
                for li_idx in range(bag.getNumChildren()):
                    to_iterate_database_ids.append(bag.getChild(li_idx).getAttributes().getValue('resource'))

            if description_tag.hasChild('isEncodedBy'):
                is_encoded_by_tag = description_tag.getChild('isEncodedBy')
                bag = is_encoded_by_tag.getChild('Bag')
                for li_idx in range(bag.getNumChildren()):
                    to_iterate_database_ids.append(bag.getChild(li_idx).getAttributes().getValue('resource'))

            gp_db_ids = {}
            for db_id in to_iterate_database_ids:
                idx = db_id.find('http://identifiers.org/')
                if idx < 0:
                    print 'WARNING: 1 geneProduct is not of the structure http://identifiers.org/DB/ID'
                    continue

                db_id = db_id[len('http://identifiers.org/'):]
                sep_idx = db_id.find('/')
                if sep_idx < 0:
                    print 'WARNING: 2 geneProduct is not of the structure http://identifiers.org/DB/ID'
                    continue

                db_name, term_id = db_id[:sep_idx], db_id[sep_idx + 1:]

                if db_name not in gp_db_ids:
                    gp_db_ids[db_name] = []

                gp_db_ids[db_name].append(term_id)

            gp_db_ids = apply_gpd_corrections(gp_db_ids)

            self.gene_product_definitions[raw_gp.getId()] = gp_db_ids

        self.term_id_mapping = {}
        mappings_dir = '%s/evaluation/gossto/id_mappings' % my_constants.basePath
        second_keys = [f for f in os.listdir(mappings_dir) if os.path.isdir('%s/%s' % (mappings_dir, f))]
        for sk in second_keys:
            sk_mapping_dir = '%s/%s' % (mappings_dir, sk)
            species_plus_first_keys = [f for f in os.listdir(sk_mapping_dir) if os.path.isfile('%s/%s' % (sk_mapping_dir, f))]
            for spfk in species_plus_first_keys:
                if spfk.startswith(species):
                    fk = spfk[len(species) + 1:]
                    mapping_conf = eval_constants.go_idmapping_conf[sk]

                    src_id_col, dest_id_col = mapping_conf['src_id_col'], mapping_conf['dest_id_col']
                    if 'species_based_id_col_mapping' in mapping_conf and species in mapping_conf['species_based_id_col_mapping']:
                        overrider = mapping_conf['species_based_id_col_mapping'][species]
                        src_id_col, dest_id_col = overrider['src_id_col'], overrider['dest_id_col']

                    this_mapping = {}
                    self.term_id_mapping[(fk, sk)] = this_mapping
                    mapping_file = open('%s/%s' % (sk_mapping_dir, spfk))
                    for l in mapping_file:
                        l = l.strip()
                        if l == '' or any([l.startswith(ilm) for ilm in mapping_conf['ignore_line_marker']]):
                            continue
                        parts = l.split('\t')
                        # if parts[src_id_col] in this_mapping:
                        #     print 'WARNING: for %s, mapping %s->%s: duplicate id_mapping for %s (previous: %s, new: %s)' % (species, fk, sk, parts[src_id_col], this_mapping[parts[src_id_col]], parts[dest_id_col])
                        if 'dest_id_col_sep' in mapping_conf:
                            dest_value = parts[dest_id_col].split(mapping_conf['dest_id_col_sep'])[0].strip()
                        else:
                            dest_value = parts[dest_id_col]
                        this_mapping[parts[src_id_col]] = dest_value
                    mapping_file.close()


def apply_gpd_corrections(db_ids):
    res = {}

    for db_name, ids in db_ids.iteritems():
        if db_name == 'ncbigi':
            res['ncbigi'] = [id[id.index('GI:') + 3:] for id in ids]
        else:
            res[db_name] = ids

    return res


def get_metabolites_of_reaction_idx(S, mets, ri):
    return [mets[mi] for mi in range(len(S)) if S[mi][ri] != 0]


# reactions inside each mmodule are considered one rmodule. each reaction set between mmodules I and J is considered a new rmodule IJ
# for now: external metabolites are considered as one module, namely external module
def convert_metabmod_to_reactmod(metabmods, is_partial_module, S, mets, reacts):
    metab_to_module = {v: k for k in metabmods.keys() for v in metabmods[k]}

    if is_partial_module or len(mets) == 2766:  # TODO:BUG: this second part is for a bug in human network!!!
        done_metabs = set(metab_to_module.keys())
        for m in mets:
            if m not in done_metabs:
                metab_to_module[m] = eval_constants.EXTERNAL_MODULE_ID

    res = {}
    for ri, r in enumerate(reacts):
        r_metabs = get_metabolites_of_reaction_idx(S, mets, ri)

        # TODO:PARTIAL: r_mods is the signature of rmod which is a combination of mmods of terminals of
        # reaction r. but external modules may be present among these combination. they could be either
        # removed or not (in case of is_partial_module)
        r_mods = tuple(sorted(list(set([metab_to_module[m] for m in r_metabs]))))

        if eval_constants.USE_PUT_BETWEEN_REACTIONS_INTO_EXTERNAL_MODULE:
            if len(r_mods) == 1:  # reaction is completely inside an mmod (mmod may be the external module, ok no problem!)
                if r_mods not in res:
                    res[r_mods] = []
                res[r_mods].append(r)
            else:  # reaction is between more than one mmod, assign it to external module so it will be penalized
                if eval_constants.EXTERNAL_MODULE_ID not in res:
                    res[eval_constants.EXTERNAL_MODULE_ID] = []
                res[eval_constants.EXTERNAL_MODULE_ID].append(r)
        else:
            if r_mods not in res:
                res[r_mods] = []
            res[r_mods].append(r)

    return res


# besides mmodules, all reactions are marked as to be inside a mmodule or between two mmodules
# in metabmod_orgranized_reacts for each module a 3-tuple specifies reactions that are inner, wrongly-not-inner and outer accordingly
# wrongly inner reactions (with respect to a given mmmod) are those ...
def convert_reactmod_to_overlapping_metabmod(reactmods, is_partial_module, S, mets, reacts):
    react_to_module = {v: k for k in reactmods.keys() for v in reactmods[k]}

    if is_partial_module:
        done_metabs = set(react_to_module.keys())
        for r in reacts:
            if r not in done_metabs:
                react_to_module[r] = eval_constants.EXTERNAL_MODULE_ID

    metabmods = {}
    metabmod_orgranized_reacts = {}

    metabs_of_react = {}

    for rmod_name, rmod_reacts in reactmods.iteritems():
        metabmods[rmod_name] = set()  # mmod_name will be the same as rmod_name
        metabmod_orgranized_reacts[rmod_name] = ([], [], [])  # inner, wrongly-not-inner, outer

        for r in rmod_reacts:
            r_metabs = get_metabolites_of_reaction_idx(S, mets, reacts.index(r))
            metabmods[rmod_name] = metabmods[rmod_name].union(r_metabs)
            metabs_of_react[r] = r_metabs

    for r in reacts:
        if r in metabs_of_react:
            r_metabs = metabs_of_react[r]
        else:
            r_metabs = get_metabolites_of_reaction_idx(S, mets, reacts.index(r))

        for mmod_name, mmod in metabmods.iteritems():  # if all terminals of a reaction are in another module (rather than its main module) then it will be counted as OUTER EDGE for that module
            if is_partial_module and mmod_name == eval_constants.EXTERNAL_MODULE_ID:
                continue

            if mmod.issuperset(r_metabs):
                if mmod_name == react_to_module[r]:
                    metabmod_orgranized_reacts[mmod_name][0].append(r)  # inner
                else:
                    metabmod_orgranized_reacts[mmod_name][1].append(r)  # wrongly-not-inner
            elif mmod.intersection(r_metabs) != set():
                metabmod_orgranized_reacts[mmod_name][2].append(r)  # outer
            else:
                pass  # unrelated

    # TODO: wrongly-not-inner could be checked to filter them out from other modules those reactions (where they are currently marked as inner for those modules)
    return metabmods, metabmod_orgranized_reacts


def read_mmods(inf):
    f = open(inf, 'r')

    i = 0
    mods = {}
    for l in f:
        l = l.strip()
        if l == '' or l.startswith('#'):
            continue

        mods['%d' % i] = l.split(' ')
        i += 1

    f.close()

    return mods


def read_rmods(inf):
    f = open(inf, 'r')

    i = 0
    mods = {}
    for l in f:
        l = l.strip()
        if l == '' or l.startswith('#'):
            continue

        mods['%d' % i] = l.split(' ')
        i += 1

    f.close()

    return mods


# returns a list of enzyme_group where each entry is an iso-enzyme. structure of an enzyme_group is [group_tag, sub_enzyme_group1 [, ...]]
# where group_tag is either OR/AND and each sub_enzyme_group is of type enzyme_group
# NOTE: a reaction may be associated with no geneProduct. method will return None is that case
# IMPL: for each reaction, typically there is a fbc:geneProductAssociation which contains a set of fbc:geneProductRef grouped
# by nested fbc:or and fbc:and tags. (1) first fbc-enabled reaction is fetched, (2) then geneProductAssociation is extracted, (3) then
# using a recursive algorithm nested groups are traversed, (4) whilst any geneProductRef encountered is looked-up in listOfGeneProducts
# definition of sbml file, (5) then proper term_id is selected according to the species (each geneProduct is linked to several different
# gene/protein dbs where for each species one of them is the best candidate), (6) then the selected term_id is converted to the target
# term_id which can be used in GOssTo analysis, (7) and finally these converted term_ids are added to the return list. wow!
def get_annotation_term_ids_for_reaction(r, species, criteria, go_loaded_model, stats_mode=False):
    # 1
    fr = go_loaded_model.sbml_model.getReaction(r).getPlugin('fbc')

    # 2
    # TODO:MISSING: reactions without enzyme!
    z = fr.getGeneProductAssociation()
    if z is None:
        if stats_mode:
            return -1, None  # not enzymatic
        else:
            return None
    raw_gassoc_of_r = z.getAssociation()

    # 3
    if stats_mode:
        mapping_state, gene_product_associations_of_r = recursive_read_gene_product_association(raw_gassoc_of_r, species, criteria, go_loaded_model, stats_mode=True)
        return mapping_state, gene_product_associations_of_r
    else:
        gene_product_associations_of_r = recursive_read_gene_product_association(raw_gassoc_of_r, species, criteria, go_loaded_model, stats_mode=False)
        return gene_product_associations_of_r


def recursive_read_gene_product_association(assoc, species, criteria, go_loaded_model, stats_mode=False):
    assoc_type = assoc.getElementName()
    if assoc_type == 'or' or assoc_type == 'and':
        ret = [assoc_type]
        if stats_mode:
            any_missed = False
        for sub_assoc in assoc.getListOfAssociations():
            if stats_mode:
                ms, z = recursive_read_gene_product_association(sub_assoc, species, criteria, go_loaded_model, stats_mode=True)
            else:
                z = recursive_read_gene_product_association(sub_assoc, species, criteria, go_loaded_model, stats_mode=False)
            if z is not None:  # TODO:MISSING: z may be null because of missing data (improper geneProduct definition)
                ret.append(z)
            if stats_mode:
                if z is None or ms != 0:
                    any_missed = True
        if stats_mode:
            if len(ret) == 1:
                return -3, None  # missed: all parts of association are missed!
            elif any_missed:
                return -4, ret  # disrupted: some but not all missed
            else:
                return 0, ret  # ok
        else:
            if len(ret) == 1:
                return None
            else:
                return ret
    else:
        if stats_mode:
            mapping_state, res = convert_to_gossto_term_id(assoc.getGeneProduct(), species, criteria, go_loaded_model, stats_mode=True)
            return mapping_state, res
        else:
            return convert_to_gossto_term_id(assoc.getGeneProduct(), species, criteria, go_loaded_model, stats_mode=False)


def get_criteria_section(c):
    if c.startswith('go_distance'):
        return 'go_distance'
    else:
        return c


def convert_to_gossto_term_id(gene_product, species, criteria, go_loaded_model, stats_mode=False):
    term_id_mapping_for_species = eval_constants.term_id_mapping_for_species[get_criteria_section(criteria)][species]

    # TODO:MISSING: this is for the cases where geneProduct is not properly defined in SBML file (geneProduct is not associated with any genes!?)
    if gene_product not in go_loaded_model.gene_product_definitions:
        if stats_mode:
            return -2, None  # undefined geneProduct
        else:
            return None

    gpd = go_loaded_model.gene_product_definitions[gene_product]  # find the GeneProduct definition

    # TODO:IMPORTANT: each gene product may have multiple entries even in one database (iso-forms of a protein). the proper
    # strategy should consider which GOA file will be used for GO distance computation and ..., "[0]" is problematic!
    gene_product_in_sbml_id = gpd[term_id_mapping_for_species[0]][0]  # extract proper id for GeneProduct according to the species to be used against gene annotation database in GOssTo

    mapping_key = term_id_mapping_for_species
    if mapping_key[0] == mapping_key[1]:  # no real mapping needed, target is the id in sbml
        translated_gene_product = gene_product_in_sbml_id
    else:
        mapping = go_loaded_model.term_id_mapping[mapping_key]
        if gene_product_in_sbml_id not in mapping:  # TODO:MISSING: this is the case that SBML id is not found in UniprotKB Retrieve Id Mapping section (few)
            if stats_mode:
                return -3, None  # missed: mapping geneProduct source id to target id failed
            else:
                return None
        translated_gene_product = mapping[gene_product_in_sbml_id]
    if stats_mode:
        return 0, translated_gene_product  # everything ok!
    else:
        return translated_gene_product


def dump_all_ids_for_database_in_all_species(show_stats=False, species_filter=None):
    ret = {}
    ret_stats = {}
    for species in my_constants.species_sbml.keys():
        if species_filter and species not in species_filter:
            continue

        print 'DOING FOR SPECIES %s...' % species
        ret[species] = {}

        src_file = r'%s/dataset/networks/%s' % (my_constants.basePath, my_constants.species_sbml[species])

        reader = libsbml.SBMLReader()
        doc = reader.readSBML(src_file)
        m = doc.getModel()

        go_loaded_model = GoLoadedModel(m, species)
        ret_stats[species] = [go_loaded_model.gene_product_count(), go_loaded_model.undefined_gene_product_count(), {}, {}]

        for gpd in go_loaded_model.gene_product_definitions.itervalues():
            for db_name, term_ids in gpd.iteritems():
                # term_ids = apply_gpd_corrections({db_name: term_ids})[db_name]
                if db_name not in ret[species]:
                    ret[species][db_name] = set()
                    ret_stats[species][2][db_name] = 0
                    ret_stats[species][3][db_name] = 0
                ret[species][db_name].update(term_ids)
                ret_stats[species][2][db_name] += 1
                ret_stats[species][3][db_name] += len(term_ids) - 1

    out_dir = r'%s/gossto' % my_constants.evalResultPath
    my_util.mkdir_p(out_dir)

    if show_stats:
        print 'species\tgps\tundefined gps\tdbname\tgps in db\tids in db\tall one'
        for s, s_d_ids in ret.iteritems():
            for d, ids in s_d_ids.iteritems():
                print '%s\t%d\t%d\t%s\t%d\t%d\t%d' % (s, ret_stats[s][0], ret_stats[s][1], d, ret_stats[s][2][d], len(ids), ret_stats[s][3][d])
                f = open(r'%s/%s__%s__ids.txt' % (out_dir, s, d), 'w')
                f.write('\n'.join(ids))
                f.close()


def get_moredata_for_reaction(react, moredata_loaded_model):
    bigg_id = react[react.index('R_') + len('R_'):]
    if bigg_id in moredata_loaded_model.bigg:
        return moredata_loaded_model.bigg[bigg_id]
    else:
        if eval_constants.VERBOSE >= 2:
            print 'reaction bigg_id %s not found in moredata' % bigg_id
        return None


def get_moredata_for_metabolite(metab, moredata_loaded_model):
    bigg_id = metab[metab.index('M_') + len('M_'):]
    if bigg_id in moredata_loaded_model.bigg:
        return moredata_loaded_model.bigg[bigg_id]
    else:
        if eval_constants.VERBOSE >= 2:
            print 'metabolite bigg_id %s not found in moredata' % bigg_id
        return None


def get_all_template_files(dir, template):
    if not os.path.isdir(dir):
        return {}

    files = [f for f in os.listdir(dir) if os.path.isfile('%s/%s' % (dir, f))]

    res = {}
    template_re = re.compile(template)
    for f in files:
        m = template_re.match(f)
        if m:
            res[m.group(1)] = '%s/%s' % (dir, f)

    return res


def get_final_modules_files_for_method(method, species):
    if method.startswith('random'):
        method_dir = r'%s/method/random_mods' % my_constants.basePath
    else:
        method_dir = r'%s/method/%s' % (my_constants.basePath, method)

    out_dir = r'%s/%s/%s' % (my_constants.resultPath, species, method)

    if method == 'guimera':
        return {'': '%s/final_modules.txt' % out_dir}
    elif method == 'holme':
        return get_all_template_files(out_dir, 'final_modules_(.*)\.txt')
    elif method == 'muller':
        return {'': '%s/final_modules.txt' % out_dir}
    elif method == 'newman':
        return {'': '%s/final_modules.txt' % out_dir}
    elif method == 'poolman':
        return get_all_template_files(out_dir, 'final_modules_(.*)\.txt')
    elif method == 'random_mmods':
        return get_all_template_files(out_dir, 'final_modules_(.*)\.txt')
    elif method == 'random_rmods':
        return get_all_template_files(out_dir, 'final_modules_(.*)\.txt')
    elif method == 'schuster':
        return get_all_template_files(out_dir, 'final_modules_(.*)\.txt')
    elif method == 'sridharan':
        return get_all_template_files(out_dir, 'final_modules_(.*)\.txt')
    elif method == 'verwoerd':
        return get_all_template_files(out_dir, 'final_modules_(.*)\.txt')


def get_saved_final_results_for_species(species, out_dir=None, results_file_prefix=''):
    if out_dir is None:
        out_dir = '%s/%s' % (my_constants.evalResultPath, species)
    my_util.mkdir_p(out_dir)
    results_file = '%s/%sresults.txt' % (out_dir, results_file_prefix)

    if not os.path.isfile(results_file):
        drf = open(results_file, 'w')
        drf.write('{}')
        drf.close()

    already_saved_results_file = open(results_file, 'r')
    already_saved_final_results = eval(already_saved_results_file.read())
    already_saved_results_file.close()

    return already_saved_final_results


def update_saved_final_results_for_species(species, already_saved_final_results, eval_results, out_dir=None, results_file_prefix=''):
    if out_dir is None:
        out_dir = '%s/%s' % (my_constants.evalResultPath, species)
    my_util.mkdir_p(out_dir)
    results_file = '%s/%sresults.txt' % (out_dir, results_file_prefix)

    if already_saved_final_results is None:
        already_saved_final_results = {}

    criteria = sorted(list(eval_results.keys()))

    for c in criteria:
        if c not in already_saved_final_results:
            already_saved_final_results[c] = {}

        methods = sorted(list(eval_results[c].keys()))
        for m in methods:
            already_saved_final_results[c][m] = eval_results[c][m]

    i = 0
    bak_file = '%s/%sresults_%d.bak' % (out_dir, results_file_prefix, i)
    while os.path.isfile(bak_file):
        i += 1
        bak_file = '%s/%sresults_%d.bak' % (out_dir, results_file_prefix, i)
    try:
        os.rename(results_file, bak_file)
    except:
        print 'creating backup file failed! ignoring ...'

    already_saved_results_file = open(results_file, 'w')
    json.dump(already_saved_final_results, already_saved_results_file, indent=4, sort_keys=True)
    already_saved_results_file.close()


sort_order = ['id', 'name', 'namespace', 'alt_id', 'subset', 'def', 'is_a', 'relationship']


def good_term_sort(l):
    def sort_order_compare(x, y):
        try:
            xi = sort_order.index(x)
        except:
            xi = 100

        try:
            yi = sort_order.index(y)
        except:
            yi = 100

        return cmp(xi, yi)

    return sorted(l, cmp=sort_order_compare)


def write_term(term, out_file_handle, term2data):
    term_id = term['id']
    out_file_handle.write('[Term]\n')
    if term_id not in term2data:
        for k in good_term_sort(term.keys()):
            v = term[k]
            if isinstance(v, list):
                for vv in v:
                    out_file_handle.write('%s: %s\n' % (k, vv))
            else:
                out_file_handle.write('%s: %s\n' % (k, v))
    else:
        override_data = term2data[term_id]
        for k in good_term_sort(list(set(term.keys() + override_data.keys()))):
            if k in override_data:
                v = override_data[k]
            else:
                v = term[k]

            if isinstance(v, list):
                for vv in v:
                    out_file_handle.write('%s: %s\n' % (k, vv))
            else:
                out_file_handle.write('%s: %s\n' % (k, v))
    out_file_handle.write('\n')


def augment_obo_file(in_obo_file, out_obo_file, term2data):
    go_data = []
    go_f = open(in_obo_file, 'r')
    go_of = open(out_obo_file, 'w')
    in_term = False
    line_number = 0
    for l_raw in go_f:
        l = l_raw.strip()
        if in_term:
            sep_idx = l.find(':')
            if sep_idx >= 0:
                key, val = l[:sep_idx].strip(), l[sep_idx + 1:].strip()
                if key in go_data[-1]:
                    aslist = go_data[-1][key]
                    if not isinstance(aslist, list):
                        aslist = [aslist]
                    go_data[-1][key] = aslist + [val]
                else:
                    go_data[-1][key] = val

        if l.startswith('[') and l.endswith(']'):  # l.lower() == '[term]' or l.lower() == '[typedef]':
            if len(go_data) > 0 and 'id' not in go_data[-1]:
                print 'found bad term (no id) on line %d' % line_number
                exit(1)

            old_in_term = in_term
            in_term = l.lower() == '[term]'

            if old_in_term:
                write_term(go_data[-1], go_of, term2data)

            if in_term:
                go_data.append({})

        if not in_term:
            go_of.write(l_raw)

        line_number += 1

    if in_term:
        write_term(go_data[-1], go_of, term2data)

    go_of.close()
    go_f.close()

    go2data = {}
    for gd in go_data:
        go2data[gd['id']] = gd

    return go2data


def read_obo_file(obo_file):
    go_data = []
    go_f = open(obo_file, 'r')
    in_term = False
    line_number = 0
    for l in go_f:
        l = l.strip()
        if in_term:
            sep_idx = l.find(':')
            if sep_idx >= 0:
                key, val = l[:sep_idx].strip(), l[sep_idx + 1:].strip()
                if key in go_data[-1]:
                    aslist = go_data[-1][key]
                    if not isinstance(aslist, list):
                        aslist = [aslist]
                    go_data[-1][key] = aslist + [val]
                else:
                    go_data[-1][key] = val

        if l.startswith('[') and l.endswith(']'):  # l.lower() == '[term]' or l.lower() == '[typedef]':
            if len(go_data) > 0 and 'id' not in go_data[-1]:
                print 'found bad term (no id) on line %d' % line_number
                exit(1)

            in_term = l.lower() == '[term]'
            if in_term:
                go_data.append({})

        line_number += 1
    go_f.close()

    go2data = {}
    for gd in go_data:
        go2data[gd['id']] = gd

    return go2data


if __name__ == '__main__':
    dump_all_ids_for_database_in_all_species(show_stats=True, species_filter=None)

    if True:
        exit()

    S = [[0, 0, 0, 1, 0], [0, -1, 0, 0, 0], [0, 1, 0, 0, 0], [0, 1, 0, 0, -1], [-1, -1, 0, 0, 0], [1, 0, 0, 0, -1], [1, 0, 0, 0, 0], [0, 0, 1, -1, 1], [0, 0, 1, 0, 0], [0, 0, -1, 0, 0], [0, 0, -1, 0, 0], ]

    metabmods = {'1': ['x5', 'x4', 'x6'], '2': ['x1', 'x2', 'x3'], '3': ['x9'], '4': ['x0', 'x7', 'x8'], '5': ['x10']}

    mets = ['x0', 'x1', 'x2', 'x3', 'x4', 'x5', 'x6', 'x7', 'x8', 'x9', 'x10']

    reacts = ['A', 'B', 'C', 'D', 'E']

    print convert_metabmod_to_reactmod(metabmods, False, S, mets, reacts)

    S = [[-1, -1, 0, 0, 0, 0, 0, 0], [0, 1, 1, -1, 0, -1, 0, 0], [0, 0, 0, 0, 0, 1, 1, 0], [0, 0, 0, 0, 0, 0, -1, 1], [1, 0, -1, 0, -1, 0, 0, -1], [0, 0, 0, 1, 1, 0, 0, 0], ]

    reactmods = {'1': ['r1', 'r2', 'r4', 'r5'], '2': ['r3', 'r6', 'r7', 'r8']}

    mets = ['m1', 'm2', 'm3', 'm4', 'm5', 'm6']

    reacts = ['r1', 'r2', 'r3', 'r4', 'r5', 'r6', 'r7', 'r8']

    print convert_reactmod_to_overlapping_metabmod(reactmods, S, mets, reacts)

    # species = 'homo_recon1'
    # species = 'ecoli_core'
    species = 'helico_iit341'
    # species = 'ecoli_ijr904'
    # species = 'ecoli_iaf1260'
    # species = 'ecoli_ijo1366'
    # species = 'saccaro_ind750'

    out_dir = r'%s/%s/newman' % (my_constants.evalResultPath, species)

    srcFile = r'%s/dataset/networks/%s' % (my_constants.basePath, my_constants.species_sbml[species])
    reader = libsbml.SBMLReader()
    doc = reader.readSBML(srcFile)
    test_go_loaded_model = GoLoadedModel(doc.getModel())

    # print get_annotation_term_ids_for_reaction('R_2OXOADOXm', 'homo_recon1', test_go_loaded_model)
    # print get_annotation_term_ids_for_reaction('R_2HCO3_NAt', 'homo_recon1', test_go_loaded_model)
    # print get_annotation_term_ids_for_reaction('R_2AMADPTm', 'homo_recon1', test_go_loaded_model)

    # print get_annotation_term_ids_for_reaction('R_FRD7', 'ecoli_core', test_go_loaded_model)

    print get_annotation_term_ids_for_reaction('R_GAPD', 'helico_iit341', test_go_loaded_model)