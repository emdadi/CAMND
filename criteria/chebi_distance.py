from operator import itemgetter

import libsbml

from util import importer
from util import my_constants
from util.my_util import get_organism
from evaluation import eval_util, eval_constants


cached_metabolite_similarities = {}
cached_compound_similarities = {}

metabolite_pair_similarity_done = 0
metabolite_pair_similarity_missed_reason1 = 0
metabolite_pair_similarity_missed_reason2 = 0


class ChebiLoadedModel:
    def __init__(self, ordered_sbml_model_list):
        """
        :param ordered_sbml_model_list: dataset/forchebi/model and then dataset/model. this is because we need forchebi because of completeness in
            metabolite->chebi_accession mapping. also we need original dataset because it has more metabolites!!!!
        """
        self.met_count = 0
        self.undefined_met_count = 0
        self.alternate_to_canonical_id_map = {}
        self.canonical_chebi_ids_for_metabolite = {}

        a2cf = open(r'%s/evaluation/chebi/alternate_to_canonical_id_map.tsv' % my_constants.basePath, 'r')
        for l in a2cf:
            if l.startswith('!'):
                continue
            lps = l.strip().split('\t')
            self.alternate_to_canonical_id_map[lps[0]] = lps[1]
        a2cf.close()

        met_definitions = {}
        for sbml_model in ordered_sbml_model_list:
            for raw_met in sbml_model.getListOfSpecies():
                if raw_met.getId() in met_definitions:
                    continue

                self.met_count += 1

                if not raw_met.getAnnotation():
                    self.undefined_met_count += 1
                    print 'WARNING: definition of metabolite %s not found' % raw_met.getId()
                    continue

                description_tag = raw_met.getAnnotation().getChild('RDF').getChild('Description')

                to_iterate_database_ids = []

                if description_tag.hasChild('is'):
                    is_tag = description_tag.getChild('is')
                    bag = is_tag.getChild('Bag')
                    for li_idx in range(bag.getNumChildren()):
                        to_iterate_database_ids.append(bag.getChild(li_idx).getAttributes().getValue('resource'))

                met_db_ids = {}
                for db_id in to_iterate_database_ids:
                    idx = db_id.find('http://identifiers.org/')
                    if idx < 0:
                        print 'WARNING: 1 metabolite is not of the structure http://identifiers.org/DB/ID'
                        continue

                    db_id = db_id[len('http://identifiers.org/'):]
                    sep_idx = db_id.find('/')
                    if sep_idx < 0:
                        print 'WARNING: 2 metabolite is not of the structure http://identifiers.org/DB/ID'
                        continue

                    db_name, term_id = db_id[:sep_idx], db_id[sep_idx + 1:]

                    if db_name not in met_db_ids:
                        met_db_ids[db_name] = []

                    met_db_ids[db_name].append(term_id)

                met_definitions[raw_met.getId()] = met_db_ids

        self.all_chebi_accessions_for_species = set()
        for met_id, met_db_ids in met_definitions.iteritems():
            canonical_accessions = None
            if 'chebi' in met_db_ids:
                chebi_accessions = met_db_ids['chebi']
                canonical_accessions = set()
                for ca in chebi_accessions:
                    if ca in self.alternate_to_canonical_id_map:
                        cca = self.alternate_to_canonical_id_map[ca]
                    else:
                        cca = ca
                    canonical_accessions.add(cca)
            else:
                self.undefined_met_count += 1

            self.canonical_chebi_ids_for_metabolite[met_id] = canonical_accessions
            if canonical_accessions:
                self.all_chebi_accessions_for_species.update(canonical_accessions)


def get_cached_compound_similarities(species, type):
    global cached_compound_similarities

    if (species, type) in cached_compound_similarities:
        similarity_table, table_compound_order = cached_compound_similarities[(species, type)]
    else:
        similarity_table = {}
        table_compound_order = []

        similarity_file = open(r'%s/evaluation/chebi/compound_similarities/%s_ism_%s.txt' % (my_constants.basePath, get_organism(species), type.upper()), 'r')
        idx = 0
        first_row = True
        for sim_line in similarity_file:
            sim_line = sim_line.strip()
            if sim_line == '' or sim_line.startswith('!'):
                continue
            if first_row:
                first_row = False
                continue

            parts = sim_line.split('\t')
            similarity_table[parts[0]] = parts[1:]
            table_compound_order.append((parts[0], idx))

            idx += 1

        # convert similarity_table to map of gene to map of gene to sim value (map-of-list converted to map-of-map)
        for compound_name, compound_sim_row in similarity_table.iteritems():
            converted_list_of_sims = {}
            for g in table_compound_order:
                converted_list_of_sims[g[0]] = float(compound_sim_row[g[1]])
            similarity_table[compound_name] = converted_list_of_sims

        table_compound_order = sorted(table_compound_order, key=itemgetter(0))

        cached_compound_similarities[(species, type)] = (similarity_table, table_compound_order)

        # print 'simtable done'

    return similarity_table, table_compound_order


def do_compute_distance(cached_metabolite_similarities, compound_similarity_table, table_compound_order, inf, species, factor1, is_rmod, is_partial_module, S, mets, reacts, revs, chebi_loaded_model):
    if is_rmod:
        mods = eval_util.read_mmods(inf)
        mods, metabmod_orgranized_reacts = eval_util.convert_reactmod_to_overlapping_metabmod(mods, is_partial_module, S, mets, reacts)
    else:
        mods = eval_util.read_mmods(inf)
        if is_partial_module:
            done_mets = set()
            for mod_n, mod_mets in mods.iteritems():
                done_mets.update(mod_mets)

            not_done_metabs = set(mets) - set(done_mets)
            if len(not_done_metabs) > 0:
                mods[eval_constants.EXTERNAL_MODULE_ID] = not_done_metabs

    # mods = filter_out_reactions_by_type(mods, [(0, 'R_EX_')])

    mods_similarity_table = {}
    for mod1_name, mod1 in mods.iteritems():
        mods_similarity_table[mod1_name] = {}
        for mod2_name, mod2 in mods.iteritems():
            sim = 0
            for m1 in mod1:
                for m2 in mod2:
                    if mod1_name != mod2_name or m1 != m2:  # designed to work for overlapping modules. only filters out a gene-similarity-with-itself for intra-module
                        sim += compute_similarity_of_metabolite_pair(cached_metabolite_similarities, m1, m2, compound_similarity_table, species, factor1, chebi_loaded_model)
            mods_similarity_table[mod1_name][mod2_name] = sim

    if eval_constants.VERBOSE >= 2:
        print 'compute_similarity_of_metabolite_pair stat for species %s: done: %d, missed_for_no_compound: %d, missed_for_not_in_chebi: %d' % (species, metabolite_pair_similarity_done, metabolite_pair_similarity_missed_reason1, metabolite_pair_similarity_missed_reason2)
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


def get_compound_ids_for_metabolite(m, chebi_loaded_model):
    compound_grp = chebi_loaded_model.canonical_chebi_ids_for_metabolite[m]
    if compound_grp == []:
        compound_grp = None
    return compound_grp


def remove_unmapped_compounds(compound_grp, simtbl):
    res = []

    for c in compound_grp:
        if c is not None and c in simtbl:
            res.append(c)

    if res == []:
        return None
    else:
        return res


def compute_similarity_of_metabolite_pair(cached_metabolite_similarities, m1, m2, sim_table, species, factor1, chebi_loaded_model):
    global metabolite_pair_similarity_done, metabolite_pair_similarity_missed_reason1, metabolite_pair_similarity_missed_reason2

    if (m1, m2, species, factor1) in cached_metabolite_similarities:  # and not stats_mode:
        return cached_metabolite_similarities[(m1, m2, species, factor1)]
    else:
        # global z
        # if z % 100000 == 0:
        #     print 'compute_similarity_of_metabolite_pair %s %s %d %d' % (species, factor1, z, len(cached_reaction_similarities))
        # z+= 1
        compound_grp1 = get_compound_ids_for_metabolite(m1, chebi_loaded_model)
        compound_grp2 = get_compound_ids_for_metabolite(m2, chebi_loaded_model)

        # TODO:MISSING:
        # reason1: metabolite has no chebi in model file!
        if compound_grp1 is None or compound_grp2 is None:
            metabolite_pair_similarity_missed_reason1 += 1
            cached_metabolite_similarities[(m1, m2, species, factor1)] = 0
            return 0

        compound_grp1 = remove_unmapped_compounds(compound_grp1, sim_table)
        compound_grp2 = remove_unmapped_compounds(compound_grp2, sim_table)

        # TODO:MISSING:
        # reason2: compounds of metabolite are not present in COA (because COA is not completely present in OBO)
        if compound_grp1 is None or compound_grp2 is None:
            metabolite_pair_similarity_missed_reason2 += 1
            cached_metabolite_similarities[(m1, m2, species, factor1)] = 0
            return 0

        # apply "best-match average" approach for the pair of metabolites
        max_similars_to_c1s = [max([sim_table[c1][c2] for c2 in compound_grp2]) for c1 in compound_grp1]
        max_similars_to_c2s = [max([sim_table[c2][c1] for c1 in compound_grp1]) for c2 in compound_grp2]
        sim = sum(max_similars_to_c1s) * 1.0 / len(max_similars_to_c1s)
        sim += sum(max_similars_to_c2s) * 1.0 / len(max_similars_to_c2s)
        sim /= 2

        cached_metabolite_similarities[(m1, m2, species, factor1)] = sim

        metabolite_pair_similarity_done += 1
        return sim


def compute_chebi_distance(inf, species, type, is_rmod, is_partial_module, S, mets, reacts, revs, chebi_loaded_model):  # type either of cc/bp/bf
    similarity_table, table_compound_order = get_cached_compound_similarities(species, type)
    return do_compute_distance(cached_metabolite_similarities, similarity_table, table_compound_order, inf, species, type, is_rmod, is_partial_module, S, mets, reacts, revs, chebi_loaded_model)


if __name__ == '__main__':
    # compute_stats_for_all_species('%s/evaluation/gossto/stats_final.txt' % my_constants.basePath, '%s/evaluation/gossto/stats_reaction_types.txt' % my_constants.basePath)
    #
    # if True:
    #     exit()

    species = 'ecoli_core'

    out_dir = r'%s/%s/newman' % (my_constants.resultPath, species)

    src_file = r'%s/dataset/networks/%s' % (my_constants.basePath, my_constants.species_sbml[species])

    S, mets, rxns, revs, met_names, rxn_names, biomass, met_comparts = importer.sbmlStoichiometricMatrix(src_file, True, read_species_compart=True)

    reader = libsbml.SBMLReader()
    doc = reader.readSBML(src_file)
    chebi_loaded_model = ChebiLoadedModel([doc.getModel()])

    print compute_chebi_distance('%s/final_modules.txt' % out_dir, species, 'mf', False, False, S, mets, rxns, revs, chebi_loaded_model)