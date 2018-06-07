## Rebinner Configuration File

# The input pattern used to glob for sample files.
input_pattern = '../SkimTreeBuilder/samples_2017/%s.root'

# The rebinner output directory.
destination = 'rebinning_0lep2017_BDTG'

# The name of the tree from which to load events.
treename = 'train'

# The name and boundary values of the distribution to rebin.
target = 'CMS_vhbb_BDTG_Znn_HighPT_13TeV', -1, 1

# Aliases for the branches required by the autocategorizer
# if they have different names than their expected defaults
# i.e. "bin", "is_signal", and "weight".
aliases = {
    'bin': 'bin_index_Znn',
    'is_signal': 'is_signal_Znn',
}

# Rebinning scheme search settings.
settings = {
    'n_bins': range(3, 7),
    'metrics': ['asimov', 'poisson'],
    'n_minbkg': [10, 20, 30, 50, 100, 200, 300],
    'smooth_bkg': [True, False],
    'unc_tol': [0.35],
}

# The samples from which to load events for rebinning.
samples = [
    'ZH125_ZNuNu_powheg',
    'ggZH125_ZNuNu_powheg',
    'WminusH125_powheg',
    'WplusH125_powheg',

    'TT_AllHadronic',
    'TT_DiLep',
    'TT_SingleLep',

    'ST_s-c_4f_lep_PSw',
    'ST_t-c_top_4f_inc',
    'ST_tW_antitop_5f_inc',
    'ST_tW_top_5f_inc_PSw',

    'WJets_madgraph',
    'WJets-HT100To200',
    'WJets-HT200To400',
    'WJets-HT400To600',
    'WJets-HT600To800',
    'WJets-HT800To1200',
    'WJets-HT1200To2500',

    'ZJetsToNuNu_HT100To200',
    'ZJetsToNuNu_HT200To400',
    'ZJetsToNuNu_HT400To600',
    'ZJetsToNuNu_HT600To800',
    'ZJetsToNuNu_HT800To1200',

    'WW',
    'WZ',
    'ZZ',

    'QCD_HT300to500',
    'QCD_HT500to700',
    'QCD_HT700to1000',
    'QCD_HT1000to1500',
    'QCD_HT1500to2000',
    'QCD_HT2000toInf',
]
