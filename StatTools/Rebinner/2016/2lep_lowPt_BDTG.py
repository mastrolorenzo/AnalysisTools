## Rebinner Configuration File

# The input pattern used to glob for sample files.
input_pattern = '../SkimTreeBuilder/samples_2016/%s.root'

# The rebinner output directory.
destination = 'rebin_2lep_lowPt_BDTG_2016'

# The name of the tree from which to load events.
treename = 'train'

# The name and boundary values of the distribution to rebin.
target = 'CMS_vhbb_BDTG_Zll_LowPT_13TeV', -1, 1

# Aliases for the branches required by the autocategorizer
# if they have different names than their expected defaults
# i.e. "bin", "is_signal", and "weight".
aliases = {
    'bin': 'bdt_index_Zll_lowPt',
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

    'TT_powheg',

    'TToLeptons_s',
    'TBarToLeptons_t_powheg',
    'TToLeptons_t_powheg',
    'Tbar_tW',
    'T_tW',

    'DYToLL_madgraph',
    'DYToLL_HT100to200',
    'DYToLL_HT200to400',
    'DYToLL_HT400to600',
    'DYToLL_HT600to800',
    'DYToLL_HT800to1200',
    'DYToLL_HT1200to2500',
    'DYToLL_HT2500toInf',
    'DYBJets-Pt100To200',
    'DYBJets-Pt200ToInf',
    'DYJets-BGenFilter-Pt100To200',
    'DYJets-BGenFilter-Pt200ToInf',

    'WW',
    'WZ',
    'ZZ',
]

