#__author__ = 'fatemeh'
from operator import itemgetter
import os

import libsbml

import my_constants
import importer
from cc_go_util import *
from my_util import get_organism


cached_reaction_similarities = {}
cached_similarity_tables = {}

def preprocess(inf,species,is_mmod):
    src_file = r'%s' % (my_constants.species_sbml[species])
    reader = libsbml.SBMLReader()
    doc = reader.readSBML(src_file)
    go_loaded_model = eval_util.GoLoadedModel(doc.getModel(), species)
    moredata_loaded_model = eval_util.MoreDataLoadedModel(species)
    sbml_file = r'%s/dataset/networks/%s' % (my_constants.basePath, my_constants.species_sbml[species])
    S, mets, rxns, revs, met_names, rxn_names, biomass, met_comparts \
        = importer.sbmlStoichiometricMatrix(sbml_file, True, read_species_compart=True)# remove_biomass=method_conf['remove_biomass'], normalize_stoich=method_conf['normalize_stoich'])
    go1=compute_coexpression_of_enzymes_aggregate(inf,species,is_mmod,True,S,mets,rxns,revs,go_loaded_model, moredata_loaded_model)
    return go1
def get_cached_similarity_tables(species, condition):
    global cached_similarity_tables

    if (species, condition) in cached_similarity_tables:
        similarity_table, table_gene_order = cached_similarity_tables[(species, condition)]
    else:
        similarity_table = {}
        table_gene_order = []

        similarity_file = open(r'%s/evaluation/gene_expr/coexp_%s_%s.tsv' % (my_constants.basePath, get_organism(species), condition), 'r')
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
                #if g[0]>len(converted_list_of_sims):
                 #   print('g[0]',g[0])
                  #  print('len_convert',len(converted_list_of_sims))
                #if g[1]> len(gene_sim_row):
                 #   print('g[1]',g[1])
                  #  print('len_gene',len(gene_sim_row))
                converted_list_of_sims[g[0]] = float(gene_sim_row[g[1]])
            similarity_table[gene_name] = converted_list_of_sims

        table_gene_order = sorted(table_gene_order, key=itemgetter(0))

        cached_similarity_tables[(species, condition)] = (similarity_table, table_gene_order)

        # print 'simtable done'

    return similarity_table, table_gene_order


def compute_coexpression_of_enzymes_in_condition(inf, species, condition, is_mmod, is_partial_module, S, mets, reacts, revs, go_loaded_model, moredata_loaded_model):
    similarity_table, table_gene_order = get_cached_similarity_tables(species, condition)
    return general_compute_distance(cached_reaction_similarities, similarity_table, table_gene_order, inf, species, 'coexpression_of_enzymes', condition, "DUMMY", is_mmod, is_partial_module, S, mets, reacts, revs, go_loaded_model, moredata_loaded_model)


def compute_coexpression_of_enzymes_aggregate(inf, species, is_mmod, is_partial_module, S, mets, reacts, revs, go_loaded_model, moredata_loaded_model):
    all_conditions_dir = '%s/evaluation/gene_expr/%s' % (my_constants.basePath, get_organism(species))
    all_conditions = [c for c in os.listdir(all_conditions_dir) if os.path.isdir('%s/%s' % (all_conditions_dir, c))]
    coexp_results = {}
    for condition in all_conditions:
        coexp_results[condition] = compute_coexpression_of_enzymes_in_condition(inf, species, condition, is_mmod, is_partial_module, S, mets, reacts, revs, go_loaded_model, moredata_loaded_model)
    return sum(coexp_results.values()) * 1.0 / len(coexp_results)


if __name__ == '__main__':
    species = 'ecoli_core'

    out_dir = r'%s/%s/newman' % (my_constants.resultPath, species)

    src_file = r'%s/dataset/networks/%s' % (my_constants.basePath, my_constants.species_sbml[species])

    S, mets, rxns, revs, met_names, rxn_names, biomass, met_comparts = importer.sbmlStoichiometricMatrix(src_file, True, read_species_compart=True)

    reader = libsbml.SBMLReader()
    doc = reader.readSBML(src_file)
    go_loaded_model = eval_util.GoLoadedModel(doc.getModel(), species)
    moredata_loaded_model = eval_util.MoreDataLoadedModel(species)
    print preprocess('/Users/fatemeh/Documents/comparison_1/RESULTS/ecoli_core/poolman/final_modules_0.5.txt','ecoli_core',True)
    #print compute_coexpression_of_enzymes_aggregate('%s/final_modules.txt' % out_dir, 'ecoli_core', True, False, S, mets, rxns, revs, go_loaded_model, moredata_loaded_model)