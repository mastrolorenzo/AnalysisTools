## SkimTreeBuilder Configuration File
from new_branches import *

# The input pattern used to glob for sample output files.
input_pattern = '/eos/uscms/store/group/lpchbb/VHbbAnalysisNtuples/2016V4_June4_incJES/%s/output*.root'

# The output directory for skimmed trees.
destination = 'skims_for_rebinning_2016'

# The list of existing branches to be kept in the skimmed
# tree in addition to the default subset of branches.
keep_branches = [
    'hJets_btagged_1',
    'H_mass',
    'V_pt',
    '*_fit_fallback',
    'CMS_vhbb_BDTG_Znn_HighPT_13TeV',
    'CMS_vhbb_BDTG_Wln_13TeV',
    'CMS_vhbb_BDTG_Zll_LowPT_13TeV',
    'CMS_vhbb_BDTG_Zll_HighPT_13TeV',
    'CMS_vhbb_DNN_Znn_13TeV',
    'CMS_vhbb_DNN_Wen_13TeV',
    'CMS_vhbb_DNN_Wmn_13TeV',
    'CMS_vhbb_DNN_Zll_LowPT_13TeV',
    'CMS_vhbb_DNN_Zll_HighPT_13TeV',
]

# The list of new branches to be created for the skimmed tree.
new_branches = [
    # Convert DNN to DNN^3
    CMS_vhbb_DNN_Znn_13TeV_F,
    CMS_vhbb_DNN_Wen_13TeV_F,
    CMS_vhbb_DNN_Wmn_13TeV_F,
    CMS_vhbb_DNN_Zll_LowPT_13TeV_F,
    CMS_vhbb_DNN_Zll_HighPT_13TeV_F,
    # 0-lepton rebinner branches
    is_signal_Znn_I,
    bdt_index_Znn_I,
    dnn_index_Znn_I,
    # 1-lepton rebinner branches
    is_signal_Wln_I,
    bdt_index_Wln_I,
    dnn_index_Wen_I,
    dnn_index_Wmn_I,
    # 2-lepton rebinner branches
    is_signal_Zll_I,
    bdt_index_Zll_lowPt_I,
    bdt_index_Zll_highPt_I,
    dnn_index_Zll_lowPt_I,
    dnn_index_Zll_highPt_I,
]

# The samples for which to build skimmed trees.
# name : [scale_factor, [input_tokens]]
samples = {
    'ZH125_ZNuNu_powheg':              [None, ['ZH125_ZNuNu_powheg']],
    'ggZH125_ZNuNu_powheg':            [None, ['ggZH125_ZNuNu_powheg']],
    'WminusH125_powheg':               [None, ['WminusH125_powheg']],
    'WplusH125_powheg':                [None, ['WplusH125_powheg']],
    'ZH125_ZLL_powheg':                [None, ['ZH125_ZLL_powheg']],
    'ggZH125_ZLL_powheg':              [None, ['ggZH125_ZLL_powheg']],

    'TT_powheg':                       [None, ['TT_powheg']],

    'TToLeptons_s':                    [None, ['TToLeptons_s']],
    'TBarToLeptons_t_powheg':          [None, ['TBarToLeptons_t_powheg']],
    'TToLeptons_t_powheg':             [None, ['TToLeptons_t_powheg']],
    'T_tW':                            [None, ['T_tW']],
    'Tbar_tW':                         [None, ['Tbar_tW']],

    'DYToLL_madgraph':                 [None, ['DYToLL_madgraph']],
    'DYToLL_HT100to200':               [None, ['DYToLL_HT100to200']],
    'DYToLL_HT200to400':               [None, ['DYToLL_HT200to400']],
    'DYToLL_HT400to600':               [None, ['DYToLL_HT400to600']],
    'DYToLL_HT600to800':               [None, ['DYToLL_HT600to800']],
    'DYToLL_HT800to1200':              [None, ['DYToLL_HT800to1200']],
    'DYToLL_HT1200to2500':             [None, ['DYToLL_HT1200to2500']],
    'DYToLL_HT2500toInf':              [None, ['DYToLL_HT2500toInf']],

    'WJets_madgraph':                  [None, ['WJets_madgraph']],
    'WJets-HT100To200':                [None, ['WJets-HT100To200']],
    'WJets-HT200To400':                [None, ['WJets-HT200To400']],
    'WJets-HT400To600':                [None, ['WJets-HT400To600']],
    'WJets-HT600To800':                [None, ['WJets-HT600To800']],
    'WJets-HT800To1200':               [None, ['WJets-HT800To1200']],
    'WJets-HT1200To2500':              [None, ['WJets-HT1200To2500']],
    'WJets-HT2500ToInf':               [None, ['WJets-HT2500ToInf']],
    'WBJets-Pt100To200':               [None, ['WBJets-Pt100To200']],
    'WBJets-Pt200ToInf':               [None, ['WBJets-Pt200ToInf']],
    'WJets_BGenFilter-Pt100To200':     [None, ['WJets_BGenFilter-Pt100To200']],
    'WJets_BGenFilter-Pt200ToInf':     [None, ['WJets_BGenFilter-Pt200ToInf']],

    'ZJetsToNuNu_HT100To200':          [None, ['ZJetsToNuNu_HT100To200']],
    'ZJetsToNuNu_HT200To400':          [None, ['ZJetsToNuNu_HT200To400']],
    'ZJetsToNuNu_HT400To600':          [None, ['ZJetsToNuNu_HT400To600']],
    'ZJetsToNuNu_HT600To800':          [None, ['ZJetsToNuNu_HT600To800']],
    'ZJetsToNuNu_HT800To1200':         [None, ['ZJetsToNuNu_HT800To1200']],
    'ZJetsToNuNu_HT1200To2500':        [None, ['ZJetsToNuNu_HT1200To2500']],
    'ZJetsToNuNu_HT2500ToInf':         [None, ['ZJetsToNuNu_HT2500ToInf']],

    'WW':                              [None, ['WW']],
    'WZ':                              [None, ['WZ']],
    'ZZ':                              [None, ['ZZ']],
}

