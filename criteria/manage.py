#__author__ = 'fatemeh'
import modularity1
import go_distance_1
import efficacy_1
import coexpression_1
import cohesion_coupling_1
import module_count
import chebi_distance_1
def process(inf,is_mmod,species,criteria):
    f1=open('result.txt','w')
    r=0.0
    if 'modularity' in criteria:
        r=modularity1.preprocess(inf,is_mmod,species)
        f1.write('modularity\t')
        f1.write(str(r))
        f1.write('\n')
    if 'module_count' in criteria:
        r=module_count.compute_module_count(inf,is_mmod)
        f1.write('module_count\t')
        f1.write(str(r))
        f1.write('\n')
    if 'go_distance_bp_F' in criteria:
        if species!='arabidopsis_irs1597':
            r=go_distance_1.preprocess(inf,species,'go_distance_bp_F',is_mmod)
            f1.write('go_distance_bp_F\t')
            f1.write(str(r))
            f1.write('\n')
        else:
            f1.write('go_distance_bp_F\t')
            f1.write('---')
            f1.write('\n')
    if 'go_distance_bp_G' in criteria:
        if species not in ['arabidopsis_irs1597','homo_recon1','mbarkeri_iaf692','ecoli_iaf1260','ecoli_ijo1366']:
            r=go_distance_1.preprocess(inf,species,'go_distance_bp_G',is_mmod)
            f1.write('go_distance_bp_G\t')
            f1.write(str(r))
            f1.write('\n')
        else:
            f1.write('go_distance_bp_G\t')
            f1.write('---')
            f1.write('\n')
    if 'go_distance_cc_G' in criteria:
        if species not in ['arabidopsis_irs1597','homo_recon1','mbarkeri_iaf692','ecoli_iaf1260','ecoli_ijo1366']:
            r=go_distance_1.preprocess(inf,species,'go_distance_cc_G',is_mmod)
            f1.write('go_distance_cc_G\t')
            f1.write(str(r))
            f1.write('\n')
        else:
            f1.write('go_distance_cc_G\t')
            f1.write('---')
            f1.write('\n')
    if 'go_distance_mf_G' in criteria:
        if species not in ['arabidopsis_irs1597','homo_recon1','mbarkeri_iaf692','ecoli_iaf1260','ecoli_ijo1366']:
            r=go_distance_1.preprocess(inf,species,'go_distance_mf_G',is_mmod)
            f1.write('go_distance_mf_G\t')
            f1.write(str(r))
            f1.write('\n')
        else:
            f1.write('go_distance_mf_G\t')
            f1.write('---')
            f1.write('\n')
    if 'go_distance_cc_F' in criteria:
        if species not in ['arabidopsis_irs1597','ecoli_iaf1260','ecoli_ijo1366']:
            r=go_distance_1.preprocess(inf,species,'go_distance_cc_F',is_mmod)
            f1.write('go_distance_cc_F\t')
            f1.write(str(r))
            f1.write('\n')
        else:
            f1.write('go_distance_cc_F\t')
            f1.write('---')
            f1.write('\n')
    if 'go_distance_mf_F' in criteria:
        r=go_distance_1.preprocess(inf,species,'go_distance_mf_F',is_mmod)
        f1.write('go_distance_mf_F\t')
        f1.write(str(r))
        f1.write('\n')
    if 'efficacy' in criteria:
        r=efficacy_1.preprocess(inf,is_mmod)
        f1.write('efficacy\t')
        f1.write(str(r))
        f1.write('\n')
    if 'cohesion_coupling' in criteria:
        if True:
            f1.write('cohesion_coupling\t')
            f1.write('---')
            f1.write('\n')
        else:
            r=cohesion_coupling_1.preprocess(inf,species,is_mmod)
            f1.write('cohesion_coupling\t')
            f1.write(str(r))
            f1.write('\n')
    if  'coexpression_of_enzymes' in criteria:
        if species!='arabidopsis_irs1597':
            r=coexpression_1.preprocess(inf,species,is_mmod)
            f1.write('coexpression_of_enzymes\t')
            f1.write(str(r))
            f1.write('\n')
        else:
            f1.write('coexpression_of_enzymes\t')
            f1.write('---')
            f1.write('\n')
    if 'chebi_distance' in criteria:
        if species!='arabidopsis_irs1597':
            r=chebi_distance_1.preprocess(inf,species,is_mmod)
            f1.write('chebi_distance\t')
            f1.write(str(r))
            f1.write('\n')
        else:
            f1.write('chebi_distance\t')
            f1.write('--- ')
            f1.write('\n')
    f1.close()
if __name__ == '__main__':
    #process('/Users/fatemeh/Documents/comparison_1/RESULTS/ecoli_core/guimera/final_modules.txt',True,'ecoli_core',['chebi_distance','efficacy','modularity','module_count'])
    process('/Users/fatemeh/Documents/comparison_1/RESULTS/helico_iit341/poolman/final_modules_0.5.txt',False,'helico_iit341',['go_distance_bp_F','go_distance_bp_G'])
