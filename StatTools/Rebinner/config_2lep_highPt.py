## Rebinner Configuration File

# The input pattern used to glob for sample files.
input_pattern = '../SkimTreeBuilder/samples/%s.root'

# The rebinner output directory.
destination = 'rebinning_2lep_highPt_scan1'

# The name of the tree from which to load events.
treename = 'train'

# The name and boundary values of the distribution to rebin.
target = 'CMS_vhbb_BDT_Zll_HighPT_13TeV', -1, 1

# Aliases for the branches required by the autocategorizer
# if they have different names than their expected defaults
# i.e. "bin", "is_signal", and "weight".
aliases = {
    'bin': 'bin_index_Zll_highPt',
    'is_signal': 'is_signal_Zll',
}

# Rebinning scheme search settings.
settings = {
    'n_bins': range(3, 9),
    'metrics': ['asimov', 'poisson'],
    'n_minbkg': [10, 30, 50, 100, 300, 500, 1000],
    'smooth_bkg': [True, False],
    'unc_tol': [0.35],
    'width_tol': [0.01, 0.025, 0.05],
}

# The samples from which to load events for rebinning.
samples = [
    'ZH125_ZLL_powheg', 
    'ggZH125_ZLL_powheg',

    'TT_powheg',

    'TToLeptons_s',
    'TBarToLeptons_t_powheg',
    #'T_tW',
    'Tbar_tW',

    'DYToLL_HT100to200',
    'DYToLL_HT200to400',
    'DYToLL_HT400to600',
    'DYToLL_HT600to800',
    'DYToLL_HT800to1200',
    'DYToLL_HT1200to2500',
    'DYToLL_HT2500toInf',

    'WJets_madgraph',
    'WJets-HT100To200',
    'WJets-HT200To400',
    'WJets-HT400To600',
    'WJets-HT600To800',
    'WJets-HT800To1200',
    'WJets-HT1200To2500',
    'WJets-HT2500ToInf',
    'WBJets-Pt100To200',
    'WBJets-Pt200ToInf',
    'WJets_BGenFilter-Pt100To200',
    'WJets_BGenFilter-Pt200ToInf',

    'WW',
    'WZ',
    'ZZ',
]

