from operator import itemgetter

import traceback
import re
import pprint as pp
import operator

import libsbml

from util.my_util import get_organism
from util import my_constants
from util import importer
from cc_go_util import *


cached_reaction_similarities = {}
cached_similarity_tables = {}


def get_cached_similarity_tables(species, type, restric_to_reliable_genes):
    global cached_similarity_tables

    if (species, type, restric_to_reliable_genes) in cached_similarity_tables:
        similarity_table, table_gene_order = cached_similarity_tables[(species, type, restric_to_reliable_genes)]
    else:
        similarity_table = {}
        table_gene_order = []

        similarity_file = open(r'%s/evaluation/gossto/gene_similarities/%s_%s_ism_%s.txt' % (my_constants.basePath, 'GOOD' if restric_to_reliable_genes else 'FULL', get_organism(species), type.upper()), 'r')

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
            table_gene_order.append((parts[0], idx))

            idx += 1

        # convert similarity_table to map of gene to map of gene to sim value (map-of-list converted to map-of-map)
        for gene_name, gene_sim_row in similarity_table.iteritems():
            converted_list_of_sims = {}
            for g in table_gene_order:
                converted_list_of_sims[g[0]] = float(gene_sim_row[g[1]])
            similarity_table[gene_name] = converted_list_of_sims

        table_gene_order = sorted(table_gene_order, key=itemgetter(0))

        cached_similarity_tables[(species, type, restric_to_reliable_genes)] = (similarity_table, table_gene_order)

        # print 'simtable done'

    return similarity_table, table_gene_order


def compute_go_distance(inf, species, type, restric_to_reliable_genes, is_mmod, is_partial_module, S, mets, reacts, revs, go_loaded_model, moredata_loaded_model):  # type either of cc/bp/bf
    global cached_reaction_similarities

    similarity_table, table_gene_order = get_cached_similarity_tables(species, type, restric_to_reliable_genes)
    return general_compute_distance(cached_reaction_similarities, similarity_table, table_gene_order, inf, species, 'go_distance', type, restric_to_reliable_genes, is_mmod, is_partial_module, S, mets, reacts, revs, go_loaded_model, moredata_loaded_model)


def compute_stats_for_all_species(res_path, reaction_types_res_path, species_filter=None):
    res = open(res_path, 'w')
    res.write('species\ttype\tall\tnot enzymatic\tundefined enzymes\tmissed: from sbml to simtbl\tdisrupted: from sbml to simtbl\trescued by EC Number\n')

    reaction_types_res = open(reaction_types_res_path, 'w')

    for species, species_file in my_constants.species_sbml.iteritems():
        if species_filter and species not in species_filter:
            continue

        reaction_types = {}
        species_row_counter = 0

        for type in ['mf', 'bp', 'cc']:
            for restric_to_reliable_genes in [False, True]:
                by_reaction_type_classifier = []

                try:
                    similarity_table, table_gene_order = get_cached_similarity_tables(species, type, restric_to_reliable_genes)
                except:
                    traceback.print_exc()
                    continue

                S, mets, rxns, revs, met_names, rxn_names, biomass, met_comparts = importer.sbmlStoichiometricMatrix(r'%s/dataset/networks/%s' % (my_constants.basePath, species_file), True, read_species_compart=True)

                src_file = r'%s/dataset/networks/%s' % (my_constants.basePath, species_file)
                reader = libsbml.SBMLReader()
                doc = reader.readSBML(src_file)
                go_loaded_model = eval_util.GoLoadedModel(doc.getModel(), species)
                moredata_loaded_model = eval_util.MoreDataLoadedModel(species)

                all_rxns, not_enzymatic, undefined_gpd, disrupted, missed, ok, rescued = 0, 0, 0, 0, 0, 0, 0
                for ri, r in enumerate(rxns):
                    all_rxns += 1
                    enz_mapping_state, enz_grp = get_even_ecnumber_annotation_term_ids_for_reaction(r, species, 'go_distance', go_loaded_model, moredata_loaded_model, stats_mode=True)
                    # enz_mapping_state, enz_grp = eval_util.get_annotation_term_ids_for_reaction(r, species, 'go_distance', go_loaded_model, stats_mode=True)

                    # enz_all = accumulate_all_enzymes_altogether(enz_grp)
                    if enz_mapping_state == -1:
                        decided_reaction_type = -1
                        not_enzymatic += 1
                    elif enz_mapping_state == -2:
                        decided_reaction_type = -2
                        undefined_gpd += 1
                    elif enz_mapping_state == -3:
                        decided_reaction_type = -3
                        missed += 1
                    elif enz_mapping_state == -5:
                        decided_reaction_type = -5
                        rescued += 1
                    else:
                        enz_all = accumulate_all_enzymes_altogether(enz_grp)
                        if all([e not in similarity_table for e in enz_all]):
                            decided_reaction_type = -3
                            missed += 1
                        elif any([e not in similarity_table for e in enz_all]):     # enz_mapping_state == -4:
                            decided_reaction_type = -4
                            disrupted += 1
                        else:
                            decided_reaction_type = 0
                            ok += 1

                    reaction_type = re.findall('[a-z]+$|^R_EX_|^R_DM_|t2$', r)
                    if reaction_type:
                        rt = reaction_type[0]
                    else:
                        rt = 'unknown'
                    if rt not in reaction_types:
                        reaction_types[rt] = set()
                    reaction_types[rt].add((r, decided_reaction_type, tuple([m for mi, m in enumerate(mets) if S[mi][ri] < 0]), tuple([m for mi, m in enumerate(mets) if S[mi][ri] > 0])))
                    by_reaction_type_classifier.append((decided_reaction_type, rt))

                # for r1 in rxns:
                #     for r2 in rxns:
                #         all += 1
                #         code, desc = compute_similarity_of_gene_pair(r1, r2, similarity_table, species, type, restric_to_reliable_genes, go_loaded_model, stats_mode=True)
                #         if code == 0:
                #             ok += 1
                #         elif code == -1:
                #             not_enzymatic += 1
                #         elif code == -2:
                #             disrupted += 1
                #         elif code == -3:
                #             missed += 1
                #         else:
                #             print 'ERROR: unknown go_distance:compute_similarity_of_gene_pair status code!'
                #             exit(1)
                # res.write('%s\t%s\t%d\t%d\t%d\t%d\n' % (species, '%s_%s' % (type.lower(), 'g' if restric_to_reliable_genes else 'f'), all_rxns, not_enzymatic, disrupted, missed))

                ne_distribution = {}
                for ft in reaction_types.keys():
                    ne_distribution[ft] = 0
                for rtc in by_reaction_type_classifier:
                    if rtc[0] == -1:
                        ne_distribution[rtc[1]] += 1

                found_types0 = sorted(list(reaction_types.keys()), key=lambda x: ne_distribution[x], reverse=True)
                found_types = []
                try:
                    rexi = found_types0.index('R_EX_')
                    found_types.append(found_types0[rexi])
                except:
                    pass

                try:
                    ti = found_types0.index('t')
                    found_types.append(found_types0[ti])
                except:
                    pass

                found_types.append('SUM OTHERS')
                
                sum_others = 0
                for ft in found_types0:
                    if ft not in {'R_EX_', 't'}:
                        found_types.append(ft)
                        sum_others += ne_distribution[ft]

                ne_distribution['SUM OTHERS'] = sum_others

                ne_distribution_line = ''
                if species_row_counter == 0:
                    ne_distribution_line = '\t'.join(found_types)
                elif species_row_counter == 1:
                    ne_distribution_line = '\t'.join([str(ne_distribution[k]) for k in found_types])

                species_row_counter += 1

                res.write('%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t\t\t\t\t%s\n' % (species, '%s_%s' % (type.lower(), 'g' if restric_to_reliable_genes else 'f'), all_rxns, not_enzymatic, undefined_gpd, missed, disrupted, rescued, ne_distribution_line))

        reaction_types_res.write('\n--------------------------------------------FOR %s:\n' % species)
        reaction_types_res.write(' '.join(sorted(list(reaction_types.keys()))) + '\n')
        for k, v in reaction_types.iteritems():
            stat = reduce(lambda x, y: tuple(map(operator.add, x, y)), [(1, 0, 0, 0, 0, 0) if z[1] == -1 else (0, 1, 0, 0, 0, 0) if z[1] == -2 else (0, 0, 1, 0, 0, 0) if z[1] == -3 else (0, 0, 0, 1, 0, 0) if z[1] == -4 else (0, 0, 0, 0, 1, 0) if z[1] == -5 else (0, 0, 0, 0, 0, 1) for z in v])
            reaction_types_res.write(k + ': ne: %d, undef: %d, mis: %d, dis: %d, rescued: %d, ok: %d\n' % (stat[0], stat[1], stat[2], stat[3], stat[4], stat[5]))
            reaction_types_res.write('     ')
            pp.pprint(v, reaction_types_res)

    res.close()
    reaction_types_res.close()


if __name__ == '__main__':
    species_filt = None
    species_filt = ['mus_imm1415']
    compute_stats_for_all_species('%s/evaluation/gossto/stats_final.txt' % my_constants.basePath, '%s/evaluation/gossto/stats_reaction_types.txt' % my_constants.basePath, species_filt)

    if True:
        exit()

    species = 'ecoli_core'

    out_dir = r'%s/%s/newman' % (my_constants.resultPath, species)

    S, mets, rxns, revs, met_names, rxn_names, biomass, met_comparts = importer.sbmlStoichiometricMatrix(r'%s/dataset/networks/%s' % (my_constants.basePath, my_constants.species_sbml[species]), True, read_species_compart=True)

    print compute_go_distance('%s/final_modules.txt' % out_dir, 'mf', True, False, S, mets, rxns, revs, GO_LOADED_MODEL, MOREDATA_LOADED_MODEL)