
criteria_conf = {
    'modularity': {'is_mmod': True, 'is_rmod': False, 'short_name': 'mod', 'is_pvalue': False},
    'go_distance_cc_G': {'is_mmod': False, 'is_rmod': True, 'short_name': 'g_c_G', 'is_pvalue': True},
    'go_distance_bp_G': {'is_mmod': False, 'is_rmod': True, 'short_name': 'g_p_G', 'is_pvalue': True},
    'go_distance_mf_G': {'is_mmod': False, 'is_rmod': True, 'short_name': 'g_f_G', 'is_pvalue': True},
    'go_distance_cc_F': {'is_mmod': False, 'is_rmod': True, 'short_name': 'g_c_F', 'is_pvalue': True},
    'go_distance_bp_F': {'is_mmod': False, 'is_rmod': True, 'short_name': 'g_p_F', 'is_pvalue': True},
    'go_distance_mf_F': {'is_mmod': False, 'is_rmod': True, 'short_name': 'g_f_F', 'is_pvalue': True},
    'coexpression_of_enzymes': {'is_mmod': False, 'is_rmod': True, 'short_name': 'coexp', 'is_pvalue': True},
    'cohesion_coupling': {'is_mmod': False, 'is_rmod': True, 'short_name': 'cc', 'is_pvalue': True},
    'module_count': {'is_mmod': True, 'is_rmod': True, 'short_name': 'cnt', 'is_pvalue': False},
    'size_distribution': {'is_mmod': True, 'is_rmod': True, 'short_name': 'sz_dst', 'is_pvalue': False},
    'efficacy': {'is_mmod': True, 'is_rmod': True, 'short_name': 'fcy', 'is_pvalue': False},
    'chebi_distance_mf': {'is_mmod': True, 'is_rmod': False, 'short_name': 'chebi', 'is_pvalue': True},
    'kegg_pathway_agreement': {'is_mmod': True, 'is_rmod': True, 'short_name': 'kegg', 'is_pvalue': False},
}

methods_to_evaluate = {
    # 'geryk': {},
    'guimera': {'is_mmod': True, 'remove_biomass': True, 'normalize_stoich': True, 'modularity_do_bipartite_explode': True, 'is_partial_module': False},
    'holme': {'is_mmod': True, 'remove_biomass': True, 'normalize_stoich': True, 'modularity_do_bipartite_explode': False, 'is_partial_module': False, 'methods_list_getter': 'holme.methods_list_getter(species)'},
    # 'kelk': {'is_mmod': False, 'remove_biomass': False, 'normalize_stoich': False, 'modularity_do_bipartite_explode': False, 'is_partial_module': True},
    'muller': {'is_mmod': False, 'remove_biomass': False, 'normalize_stoich': False, 'modularity_do_bipartite_explode': False, 'is_partial_module': True},
    'newman': {'is_mmod': True, 'remove_biomass': True, 'normalize_stoich': True, 'modularity_do_bipartite_explode': True, 'is_partial_module': False},
    'poolman': {'is_mmod': False, 'remove_biomass': True, 'normalize_stoich': True, 'modularity_do_bipartite_explode': False, 'is_partial_module': False, 'methods_list_getter': 'poolman.methods_list_getter(species)'},
    # 'random_rmods': {'is_mmod': False, 'remove_biomass': True, 'normalize_stoich': True, 'modularity_do_bipartite_explode': False, 'is_partial_module': False},
    # 'random_mmods': {'is_mmod': True, 'remove_biomass': True, 'normalize_stoich': True, 'modularity_do_bipartite_explode': False, 'is_partial_module': False},
    'schuster': {'is_mmod': True, 'remove_biomass': True, 'normalize_stoich': True, 'modularity_do_bipartite_explode': False, 'is_partial_module': True, 'methods_list_getter': 'schuster.methods_list_getter(species)'},
    'sridharan': {'is_mmod': False, 'remove_biomass': False, 'normalize_stoich': False, 'modularity_do_bipartite_explode': False, 'is_partial_module': False, 'methods_list_getter': 'sridharan.methods_list_getter(species)'},
    'verwoerd': {'is_mmod': True, 'remove_biomass': True, 'normalize_stoich': True, 'modularity_do_bipartite_explode': False, 'is_partial_module': True, 'methods_list_getter': 'verwoerd_prepare.methods_list_getter(species)'},
    # 'yoon': {}
    'Ding': {'is_mmod': True, 'remove_biomass': True, 'normalize_stoich': True, 'modularity_do_bipartite_explode': False, 'is_partial_module': True},
}

# go_annotation_term_id_for_species = {
#     'toy_model': ('?', '?'),
#     'ecoli_core': ('uniprot', 'uniprot'),
#     'helico_iit341': ('ncbigene', 'ncbigene'),
#     # 'ecoli_ijr904': 'uniprot',
#     'ecoli_iaf1260': ('ncbigi', 'uniprot'),
#     'ecoli_ijo1366': ('ncbigi', 'uniprot'),
#     'saccaro_ind750': ('sgd', 'sgd'),
#     'homo_recon1': '?',
# }

gossto_term_id_mapping_for_species = {
    'toy_model': ('?', '?'),
    'ecoli_core': ('uniprot', 'uniprot'),
    'helico_iit341': ('ncbigene', 'uniprot'),
    'ecoli_iaf1260': ('ncbigene', 'uniprot'),
    'ecoli_ijo1366': ('ncbigene', 'uniprot'),
    'saccaro_ind750': ('ncbigene', 'uniprot'),
    'homo_recon1': ('ncbigene', 'genesymbol'),
    'mbarkeri_iaf692': ('ncbigi', 'uniprot'),
    'arabidopsis_irs1597': ('ncbigene', 'ncbigene'),  # dummy!
    'mus_imm1415': ('ncbigene', 'genesymbol'),
}

coexp_term_id_mapping_for_species = {
    'toy_model': ('?', '?'),
    'ecoli_core': ('ncbigene', 'alias1'),
    'helico_iit341': ('ncbigene', 'alias1'),
    'ecoli_iaf1260': ('ncbigene', 'alias1'),
    'ecoli_ijo1366': ('ncbigene', 'alias1'),
    'saccaro_ind750': ('ncbigene', 'alias1'),
    'homo_recon1': ('ncbigene', 'alias1'),
    'mbarkeri_iaf692': ('ncbigene', 'alias1'),
    'arabidopsis_irs1597': ('ncbigene', 'ncbigene'),    # dummy
    'mus_imm1415': ('ncbigene', 'genesymbol'),
}

term_id_mapping_for_species = {
    'go_distance': gossto_term_id_mapping_for_species,
    'coexpression_of_enzymes': coexp_term_id_mapping_for_species,
}

go_idmapping_conf = {
    'uniprot': {'src_id_col': 0, 'dest_id_col': 1, 'ignore_line_marker': ['yourlist', '!'], 'species_based_id_col_mapping': {'mbarkeri_iaf692': {'src_id_col': 0, 'dest_id_col': 2}}},
    'genesymbol': {'src_id_col': 2, 'dest_id_col': 5, 'ignore_line_marker': ['tax_id']},
    'alias1': {'src_id_col': 2, 'dest_id_col': 6, 'dest_id_col_sep': ',', 'ignore_line_marker': ['tax_id']},
}

poor_go_mappings = {
    ('homo_recon1', 'go_distance_mf_G'), ('homo_recon1', 'go_distance_bp_G'),
    ('ecoli_core', 'go_distance_bp_G'),
    ('ecoli_ijo1366', 'go_distance_mf_G'), ('ecoli_ijo1366', 'go_distance_bp_G'), ('ecoli_ijo1366', 'go_distance_cc_F'), ('ecoli_ijo1366', 'go_distance_cc_G'),
    ('helico_iit341', 'go_distance_mf_G'), ('helico_iit341', 'go_distance_bp_G'), ('helico_iit341', 'go_distance_cc_F'), ('helico_iit341', 'go_distance_cc_G'),
    ('ecoli_iaf1260', 'go_distance_mf_G'), ('ecoli_iaf1260', 'go_distance_bp_G'), ('ecoli_iaf1260', 'go_distance_cc_F'), ('ecoli_iaf1260', 'go_distance_cc_G'),
    ('saccaro_ind750', 'go_distance_cc_G'), ('saccaro_ind750', 'go_distance_mf_G'), ('saccaro_ind750', 'go_distance_bp_G'),
    ('mbarkeri_iaf692', 'go_distance_cc_F'), ('mbarkeri_iaf692', 'go_distance_mf_G'), ('mbarkeri_iaf692', 'go_distance_bp_G'), ('mbarkeri_iaf692', 'go_distance_cc_G'),
    ('arabidopsis_irs1597', 'go_distance_bp_F'), ('arabidopsis_irs1597', 'go_distance_cc_F'), ('arabidopsis_irs1597', 'go_distance_mf_G'), ('arabidopsis_irs1597', 'go_distance_bp_G'), ('arabidopsis_irs1597', 'go_distance_cc_G'),
    ('mus_imm1415', 'go_distance_mf_G'), ('mus_imm1415', 'go_distance_bp_G'), ('mus_imm1415', 'go_distance_cc_G'),
}

EXTERNAL_MODULE_ID = 1000000  # for partial methods, metabolites outside all modules are marked with this module
USE_PUT_BETWEEN_REACTIONS_INTO_EXTERNAL_MODULE = True
ZSCORE_RANDOM_SAMPLE_COUNT = 100
PVALUE_STABILITY_FOLD_SIZE = 10
VERBOSE = 2
ZSCORE_PROCESS_COUNT = 10
