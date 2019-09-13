import eval_util
from evaluation.eval_util import MoreDataLoadedModel
from util import my_constants
from util import importer


class KeggLoadedModel(MoreDataLoadedModel):
    def __init__(self, species):
        MoreDataLoadedModel.__init__(self, species)  # contains mapping of metabolites/reactions to kegg_ids

        self.kegg_reactions = {}
        self.kegg_compounds = {}
        self.load_kegg_moredata()  # contains kegg_ids mapping to kegg data (esp. modules/pathways for each compound/reaction)

        self.cached_overlapping_clusterings = {}  # map from (element_type, indirect_kegg_property) to a clustering

    def load_kegg_moredata(self):
        inf = open('%s/evaluation/kegg/all_kegg_moredata.json' % my_constants.basePath, 'r')
        elements = eval(inf.read().replace(': null', ': ""'))
        for ke in elements:
            if ke['type'] == 'reaction':
                if ke['kegg_id'] in self.kegg_reactions:
                    print 'creating KeggLoadedModel for %s failed: more than one reaction with kegg_id %s' % (species, ke['kegg_id'])
                    # exit(1)
                self.kegg_reactions[ke['kegg_id']] = ke
            elif ke['type'] == 'metabolite':
                if ke['kegg_id'] in self.kegg_compounds:
                    print 'creating KeggLoadedModel for %s failed: more than one metabolites with kegg_id %s' % (species, ke['kegg_id'])
                    # exit(1)
                self.kegg_compounds[ke['kegg_id']] = ke

        inf.close()

    def get_overlapping_clustering_for_induced_by(self, element_type, indirect_kegg_property):
        if (element_type, indirect_kegg_property) in self.cached_overlapping_clusterings:
            return self.cached_overlapping_clusterings[(element_type, indirect_kegg_property)]
        else:
            clust = {}

            if element_type == 'reaction':
                is_metab = False
                bigg_id_prefix = 'R_'
                kegg_prop_in_bigg = 'KEGG Reaction'
            elif element_type == 'metabolite':
                is_metab = True
                bigg_id_prefix = 'M_'
                kegg_prop_in_bigg = 'KEGG Compound'
            else:
                print 'unknown element_type in get_base_clustering %s' % element_type
                exit(1)

            for eid, e in self.bigg.iteritems():
                if e['type'] == element_type:
                    values_of_indirect_kegg_property, state = self.__general_map_element_to_aggr_prop(is_metab, eid, None, kegg_prop_in_bigg, indirect_kegg_property)
                    for v in values_of_indirect_kegg_property:
                        if v not in clust:
                            clust[v] = set()
                        clust[v].add(bigg_id_prefix + eid)

            deduplicated_clust = {}
            while clust:
                cn, c = clust.popitem()
                deduplicated_clust[cn] = c

                cn_dup_names = []
                for cn2, c2 in clust.iteritems():
                    if c2 == c:
                        cn_dup_names.append(cn2)
                for cdn in cn_dup_names:
                    clust.pop(cdn)

            self.cached_overlapping_clusterings[(element_type, indirect_kegg_property)] = deduplicated_clust
            return deduplicated_clust

    # def __map_metabolite_to_module(self, m):
    #     return self.__general_map_element_to_aggr_prop(True, m, 'M_', 'KEGG Compound', 'modules')
    #
    # def __map_reaction_to_module(self, m):
    #     return self.__general_map_element_to_aggr_prop(False, m, 'R_', 'KEGG Reaction', 'modules')
    #
    # def __map_metabolite_to_pathway(self, m):
    #     return self.__general_map_element_to_aggr_prop(True, m, 'M_', 'KEGG Compound', 'pathways')
    #
    # def __map_reaction_to_pathway(self, m):
    #     return self.__general_map_element_to_aggr_prop(False, m, 'R_', 'KEGG Reaction', 'pathways')

    def __general_map_element_to_aggr_prop(self, is_metab, m, prefix, kegg_prop_in_bigg, aggr_prop):
        module = []

        if prefix:
            bid = m[m.index(prefix) + len(prefix):]
        else:
            bid = m

        if bid not in self.bigg:
            return [], 1  # metabolite/reaction not found in bigg_moredata.json
        moredata = self.bigg[bid]

        if kegg_prop_in_bigg not in moredata or not moredata[kegg_prop_in_bigg]:
            return [], 2  # moredata for the metabolite/reaction does not refer to any kegg compound/reaction
        partial_miss = 0
        kegg_data_part = self.kegg_compounds if is_metab else self.kegg_reactions
        for keggelem in moredata[kegg_prop_in_bigg]:
            if keggelem in kegg_data_part:
                if aggr_prop in kegg_data_part[keggelem]:
                    module += kegg_data_part[keggelem][aggr_prop]
                else:
                    partial_miss += 1
            else:
                partial_miss += 1

        if not module:
            return [], 3  # non of compounds/reactions have associated 'modules' in kegg

        return module, -partial_miss  # DOES return a list of modules. but there may be some modules missed during mapping, if so return the number of such misses (reason1: kegg_id not found in proper section, reason2: found kegg element does not have any modules)


# def pair_of_sets_distance_by_entropy(s1, s2):  # NOTE: ref. Measuring similarity between sets of overlapping clusters
#     pass  # so hard! skip!


def pair_of_sets_distance_by_d1(ss1, ss2):  # NOTE: ref. Measuring similarity between sets of overlapping clusters
    return len(ss1.union(ss2)) - len(ss1.intersection(ss2))


def pair_of_sets_distance_by_d2(ss1, ss2):  # NOTE: ref. Measuring similarity between sets of overlapping clusters
    return 1 - (len(ss1.intersection(ss2)) * 1.0 / len(ss1.union(ss2)))


def pair_of_overlapping_clusterings_distance_best_match(c1, c2, set_distance_computer):
    sum1 = 0
    for sn1, s1 in c1.iteritems():
        sum1 += min([set_distance_computer(s1, s2) for s2 in c2.itervalues()])

    sum2 = 0
    for sn2, s2 in c2.iteritems():
        sum2 += min([set_distance_computer(s2, s1) for s1 in c1.itervalues()])

    return (sum1 + sum2) * 1.0 / (len(c1) + len(c2))  # TODO: needs normalization


def compute_kegg_pathway_agreement(inf, is_rmod, is_partial_module, S, mets, reacts, revs, kegg_loaded_model):
    if is_rmod:
        mods = eval_util.read_mmods(inf)
    else:
        mods = eval_util.read_rmods(inf)

    # TODO: PARTIAL: how partial modules should be handled in KEGG agreement computations?

    if is_rmod:
        base_clustering = kegg_loaded_model.get_overlapping_clustering_for_induced_by('reaction', 'pathways')
    else:
        base_clustering = kegg_loaded_model.get_overlapping_clustering_for_induced_by('metabolite', 'modules')

    mods_converted_to_set = {}
    for mn, mod in mods.iteritems():
        mods_converted_to_set[mn] = set(mod)

    return 1 - pair_of_overlapping_clusterings_distance_best_match(mods_converted_to_set, base_clustering, pair_of_sets_distance_by_d2)


if __name__ == '__main__':
    species = 'ecoli_core'

    kegg_loaded_model = KeggLoadedModel(species)

    out_dir = r'%s/%s/newman' % (my_constants.resultPath, species)
    S, mets, rxns, revs, met_names, rxn_names, biomass, met_comparts = importer.sbmlStoichiometricMatrix(r'%s/dataset/networks/%s' % (my_constants.basePath, my_constants.species_sbml[species]), True, read_species_compart=True)
    print 'newman:', compute_kegg_pathway_agreement('%s/final_modules.txt' % out_dir, False, False, S, mets, rxns, revs, kegg_loaded_model)

    out_dir = r'%s/%s/poolman' % (my_constants.resultPath, species)
    S, mets, rxns, revs, met_names, rxn_names, biomass, met_comparts = importer.sbmlStoichiometricMatrix(r'%s/dataset/networks/%s' % (my_constants.basePath, my_constants.species_sbml[species]), True, read_species_compart=True)
    print 'poolman0.5', compute_kegg_pathway_agreement('%s/final_modules_0.5.txt' % out_dir, True, False, S, mets, rxns, revs, kegg_loaded_model)
