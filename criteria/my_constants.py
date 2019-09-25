import sys

if sys.platform.startswith('win'):
    basePath = r'D:/project/comparison'
    resultPath = r'D:/project/comparison/RESULTS'
    evalResultPath = r'D:/project/comparison/EVAL_RESULTS'
    rScriptPath = r'"C:/Program Files/R/R-3.4.1/bin/x64/Rscript.exe"'
    win = True
else:
    basePath = '/Users/fatemeh/Documents/comparison_1'
    resultPath = '/Users/fatemeh/Documents/comparison_1/RESULTS'
    evalResultPath = '/Users/fatemeh/Documents/comparison_1/EVAL_RESULTS'
    rScriptPath = r'/Users/fatemeh/bin/Rscript'
    # resultPath = r'/home/a.rezvan/fca3%s/results' % ('' if len(sys.argv) != 2 else ('/' + sys.argv[1]))
    win = False

species_sbml = {
                # "toy_model": "toy-model-kelk.xml",
                "ecoli_core": "e_coli_core.xml",
                "helico_iit341": "iIT341.xml",
                # "ecoli_ijr904": "iJR904.xml",
                "ecoli_iaf1260": "iAF1260.xml",
                "ecoli_ijo1366": "iJO1366.xml",
                "saccaro_ind750": "iND750.xml",
                "homo_recon1": "RECON1.xml",
                "mbarkeri_iaf692": "iAF692.xml",
                "arabidopsis_irs1597": "iRS1597.xml",
                "mus_imm1415": "iMM1415.xml",
}

species_artificial_biomass = {
    "toy_model": False,
    "ecoli_core": True,
    "helico_iit341": True,
    "ecoli_iaf1260": True,
    "ecoli_ijo1366": True,
    "saccaro_ind750": True,
    "homo_recon1": False,
    "mbarkeri_iaf692": True,
    "arabidopsis_irs1597": True,
    "mus_imm1415": True,
}

species_sort_order = [
    'toy_model',
    'ecoli_core',
    'helico_iit341',
    'mbarkeri_iaf692',
    'saccaro_ind750',
    'arabidopsis_irs1597',
    'ecoli_iaf1260',
    'ecoli_ijo1366',
    'mus_imm1415',
    'homo_recon1',
]
