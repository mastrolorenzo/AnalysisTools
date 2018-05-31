## Rebinner Configuration File

# The input pattern used to glob for sample files.
input_pattern = '../SkimTreeBuilder/samples_2016/%s.root'

# The rebinner output directory.
destination = 'rebinning_0lep_scan1'

# The name of the tree from which to load events.
treename = 'train'

# The name and boundary values of the distribution to rebin.
target = 'CMS_vhbb_BDT_Znn_13TeV', -1, 1

# Aliases for the branches required by the autocategorizer
# if they have different names than their expected defaults
# i.e. "bin", "is_signal", and "weight".
aliases = {
    'bin': 'bin_index_Znn',
    'is_signal': 'is_signal_Znn',
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
    'ZH125_ZNuNu_powheg',
    'ggZH125_ZNuNu_powheg',
    'WminusH125_powheg',
    'WplusH125_powheg',

    'TT_powheg',

    'TToLeptons_s',
    'TBarToLeptons_t_powheg',
    #'T_tW',
    'Tbar_tW',

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

    'ZJetsToNuNu_HT100To200',
    'ZJetsToNuNu_HT200To400',
    'ZJetsToNuNu_HT400To600',
    'ZJetsToNuNu_HT600To800',
    'ZJetsToNuNu_HT800To1200',
    'ZJetsToNuNu_HT1200To2500',
    'ZJetsToNuNu_HT2500ToInf',

    'WW',
    'WZ',
    'ZZ',
]

