## SkimTreeBuilder Configuration File
from new_branches import *
from scale_factors_HIG16044 import *

# The input pattern used to glob for sample output files.
input_pattern = '/eos/uscms/store/group/lpchbb/VHbbNanoPostProc/2016/NanoNtuples/Analyzed/2016V4_RunALL/%s/output*.root'

# The output directory for skimmed trees.
destination = 'samples_2016'

# The list of existing branches to be kept in the skimmed
# tree in addition to the default subset of branches.
keep_branches = [
    'hJets_pt_0',
    'hJets_btagged_0',
    'hJets_pt_1',
    'hJets_btagged_1',
    'HJ1_HJ2_dEta',
    'H_mass',
    'V_pt',
    'V_mass',
    'HVdPhi',
    'MET_pt',
    'MET_Pt',
    'nAddJets_2lep',
    'SA5',
    'n_recoil_jets_fit',
    '*_fit_fallback',
    'H_mass_sigma_fit',
]

# The list of new branches to be created for the skimmed tree.
new_branches = [
    # 0-lepton rebinner branches
    is_signal_Znn_I,
    bin_index_Znn_I,
    # 1-lepton rebinner branches
    is_signal_Wln_I,
    bin_index_Wln_I,
    # 2-lepton rebinner branches
    is_signal_Zll_I,
    bin_index_Zll_lowPt_I,
    bin_index_Zll_highPt_I,
]

# The samples for which to build skimmed trees.
# name : [scale_factor, [input_tokens]]
samples = {
    'Run2016BToG_DoubleEle_ReMiniAOD': [None, ['Run2016BToG_DoubleEle_ReMiniAOD']],
    'Run2016BToG_DoubleMu_ReMiniAOD':  [None, ['Run2016BToG_DoubleMu_ReMiniAOD']],
    'Run2016BToG_Ele_ReMiniAOD':       [None, ['Run2016BToG_Ele_ReMiniAOD']],
    'Run2016BToG_MET_ReMiniAOD':       [None, ['Run2016BToG_MET_ReMiniAOD']],
    'Run2016BToG_Mu_ReMiniAOD':        [None, ['Run2016BToG_Mu_ReMiniAOD']],

    'ZH125_ZNuNu_powheg':              [None, ['ZH125_ZNuNu_powheg']],
    'ggZH125_ZNuNu_powheg':            [None, ['ggZH125_ZNuNu_powheg']],
    'WminusH125_powheg':               [None, ['WminusH125_powheg']],
    'WplusH125_powheg':                [None, ['WplusH125_powheg']],
    'ZH125_ZLL_powheg':                [None, ['ZH125_ZLL_powheg']],
    'ggZH125_ZLL_powheg':              [None, ['ggZH125_ZLL_powheg']],

    'TT_powheg':                       [sf_tt, ['TT_powheg']],

    'TToLeptons_s':                    [None, ['TToLeptons_s']],
    'TBarToLeptons_t_powheg':          [None, ['TBarToLeptons_t_powheg']],
    #'T_tW':                            [None, ['T_tW']],
    'Tbar_tW':                         [None, ['Tbar_tW']],

    'DYToLL_HT100to200':               [sf_zjets, ['DYToLL_HT100to200']],
    'DYToLL_HT200to400':               [sf_zjets, ['DYToLL_HT200to400']],
    'DYToLL_HT400to600':               [sf_zjets, ['DYToLL_HT400to600']],
    'DYToLL_HT600to800':               [sf_zjets, ['DYToLL_HT600to800']],
    'DYToLL_HT800to1200':              [sf_zjets, ['DYToLL_HT800to1200']],
    'DYToLL_HT1200to2500':             [sf_zjets, ['DYToLL_HT1200to2500']],
    'DYToLL_HT2500toInf':              [sf_zjets, ['DYToLL_HT2500toInf']],
    'DYJets_1J':                       [sf_zjets, ['DYJets_1J']],
    'DYJets_2J':                       [sf_zjets, ['DYJets_2J']],
    'DYJets_3J':                       [sf_zjets, ['DYJets_3J']],
    'DYBJets-Pt100To200':              [sf_zjets, ['DYBJets-Pt100To200']],
    'DYBJets-Pt200ToInf':              [sf_zjets, ['DYBJets-Pt200ToInf']],

    'WJets_madgraph':                  [sf_wjets, ['WJets_madgraph']],
    'WJets-HT100To200':                [sf_wjets, ['WJets-HT100To200']],
    'WJets-HT200To400':                [sf_wjets, ['WJets-HT200To400']],
    'WJets-HT400To600':                [sf_wjets, ['WJets-HT400To600']],
    'WJets-HT600To800':                [sf_wjets, ['WJets-HT600To800']],
    'WJets-HT800To1200':               [sf_wjets, ['WJets-HT800To1200']],
    'WJets-HT1200To2500':              [sf_wjets, ['WJets-HT1200To2500']],
    'WJets-HT2500ToInf':               [sf_wjets, ['WJets-HT2500ToInf']],
    'WBJets-Pt100To200':               [sf_wjets, ['WBJets-Pt100To200']],
    'WBJets-Pt200ToInf':               [sf_wjets, ['WBJets-Pt200ToInf']],
    'WJets_BGenFilter-Pt100To200':     [sf_wjets, ['WJets_BGenFilter-Pt100To200']],
    'WJets_BGenFilter-Pt200ToInf':     [sf_wjets, ['WJets_BGenFilter-Pt200ToInf']],

    'ZJetsToNuNu_HT100To200':          [sf_zjets, ['ZJetsToNuNu_HT100To200']],
    'ZJetsToNuNu_HT200To400':          [sf_zjets, ['ZJetsToNuNu_HT200To400']],
    'ZJetsToNuNu_HT400To600':          [sf_zjets, ['ZJetsToNuNu_HT400To600']],
    'ZJetsToNuNu_HT600To800':          [sf_zjets, ['ZJetsToNuNu_HT600To800']],
    'ZJetsToNuNu_HT800To1200':         [sf_zjets, ['ZJetsToNuNu_HT800To1200']],
    'ZJetsToNuNu_HT1200To2500':        [sf_zjets, ['ZJetsToNuNu_HT1200To2500']],
    'ZJetsToNuNu_HT2500ToInf':         [sf_zjets, ['ZJetsToNuNu_HT2500ToInf']],

    'WW':                              [None, ['WW']],
    'WZ':                              [None, ['WZ']],
    'WZ_fil':                          [None, ['WZ_fil']],
    'ZZ':                              [None, ['ZZ']],
    'ZZTo2L2Q':                        [None, ['ZZTo2L2Q']],
    'ZZTo2Q2Nu':                       [None, ['ZZTo2Q2Nu']],
}

