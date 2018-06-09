## Rebinner Configuration File

# The input pattern used to glob for sample files.
input_pattern = '../SkimTreeBuilder/samples_2016/%s.root'

# The rebinner output directory.
destination = 'rebinning_1lep_BDTG_2016'

# The name of the tree from which to load events.
treename = 'train'

# The name and boundary values of the distribution to rebin.
target = 'CMS_vhbb_BDTG_Wln_13TeV', -1, 1

# Aliases for the branches required by the autocategorizer
# if they have different names than their expected defaults
# i.e. "bin", "is_signal", and "weight".
aliases = {
    'bin': 'bin_index_Wln',
    'is_signal': 'is_signal_Wln',
}

# Rebinning scheme search settings.
settings = {
    'n_bins': range(3, 9),
    'metrics': ['asimov', 'poisson'],
    'n_minbkg': [10, 20, 30, 50, 100, 200, 300],
    'smooth_bkg': [True, False],
    'unc_tol': [0.35],
}

# The samples from which to load events for rebinning.
samples = [
    'WminusH125_powheg',
    'WplusH125_powheg',

    'TT_powheg',

    'TToLeptons_s',
    'TBarToLeptons_t_powheg',
    'TToLeptons_t_powheg',
    'Tbar_tW',
    'T_tW',

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

