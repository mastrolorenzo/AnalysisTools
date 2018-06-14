## Rebinner Configuration File

# The input pattern used to glob for sample files.
input_pattern = '../SkimTreeBuilder/samples_2017/%s.root'

# The rebinner output directory.
destination = 'rebin_2lep_highPt_BDTG_2017'

# The name of the tree from which to load events.
treename = 'train'

# The name and boundary values of the distribution to rebin.
target = 'CMS_vhbb_BDTG_Zll_HighPT_13TeV', -1, 1

# Aliases for the branches required by the autocategorizer
# if they have different names than their expected defaults
# i.e. "bin", "is_signal", and "weight".
aliases = {
    'bin': 'bdt_index_Zll_highPt',
    'is_signal': 'is_signal_Zll',
}

# Rebinning scheme search settings.
settings = {
    'n_bins': range(3, 10),
    'metrics': ['asimov', 'poisson'],
    'n_minbkg': [10, 20, 30, 50, 100, 200, 300],
    'smooth_bkg': [True, False],
    'unc_tol': [0.35],
}

# The samples from which to load events for rebinning.
samples = [
    'ZH125_ZLL_powheg',
    'ggZH125_ZLL_powheg',

    'TT_AllHadronic',
    'TT_DiLep',
    'TT_SingleLep',

    'ST_s-c_4f_lep_PSw',
    'ST_t-c_antitop_4f_inc',
    'ST_t-c_top_4f_inc',
    'ST_tW_antitop_5f_inc',
    'ST_tW_top_5f_inc_PSw',

    'DYToLL_madgraph',
    'DYToLL_HT100to200_madgraph',
    'DYToLL_HT200to400_madgraph',
    'DYToLL_HT400to600_madgraph',
    'DYToLL_HT600to800_madgraph',
    'DYToLL_HT800to1200_madgraph',
    'DYToLL_HT1200to2500_madgraph',
    'DYToLL_HT2500toInf_madgraph',

    'DYToLL_M4to50_HT70to100_madgraph',
    'DYToLL_M4to50_HT100to200_madgraph',
    'DYToLL_M4to50_HT200to400_madgraph',
    'DYToLL_M4to50_HT400to600_madgraph',
    'DYToLL_M4to50_HT600toInf_madgraph',

    'WW',
    'WZ',
    'ZZ',
]

