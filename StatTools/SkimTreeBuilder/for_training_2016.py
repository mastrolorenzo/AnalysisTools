## SkimTreeBuilder Configuration File
from new_branches import *

# The input pattern used to glob for sample output files.
input_pattern = '/eos/uscms/store/group/lpchbb/VHbbAnalysisNtuples/CMSConnect_June12_2016V4_FixedBtagInput/haddjobs/sum_%s*root'

# The output directory for skimmed trees.
destination = 'skim_for_training_2016'

# The list of existing branches to be kept in the skimmed
# tree in addition to the default subset of branches.
keep_branches = [
    'hJetInd1',
    'hJetInd2',
    'Jet_btagCMVA',
    'hJets_leadingPt',
    'hJets_subleadingPt',
    'hJets_btagged_0',
    'hJets_btagged_1',
    'HJ1_HJ2_dEta',
    'HJ1_HJ2_dPhi',
    'HJ1_HJ2_dR',
    'H_pt',
    'H_mass',
    'V_pt',
    'V_mass',
    'V_mt',
    'MET_pt',
    'MET_Pt',
    'HVdPhi',
    'jjVPtRatio',
    'SA5',
    'otherJetsBestBtag',
    'otherJetsHighestPt',
    'minDPhiFromOtherJets',
    'Top1_mass_fromLepton_regPT_w4MET',
    'nAddJets302p5_puid',
    'lepMetDPhi',
    '*_fit_fallback',
    'nAddJets_2lep',
    'n_recoil_jets_fit',
    'H_mass_sigma_fit',
]

# The list of new branches to be created for the skimmed tree.
new_branches = [
    hJets_DeepCSV_0_F,
    hJets_DeepCSV_1_F,
    hJets_CMVA_0_F,
    hJets_CMVA_1_F,
    abs_HJ1_HJ2_dEta_F,
    abs_HVdPhi_F,
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
    'DYBJets-Pt100To200':              [None,['DYBJets-Pt100To200']],
    'DYBJets-Pt200ToInf':              [None,['DYBJets-Pt200ToInf']],
    'DYJets-BGenFilter-Pt100To200':    [None,['DYJets-BGenFilter-Pt100To200']],
    'DYJets-BGenFilter-Pt200ToInf':    [None,['DYJets-BGenFilter-Pt200ToInf']],

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

    'ZBJetsToNuNu_Pt-100to200':        [None,['ZBJetsToNuNu_Pt-100to200']],
    'ZBJetsToNuNu_Pt-200toInf':        [None,['ZBJetsToNuNu_Pt-200toInf']],
    'ZJetsToNuNu_BGenFilter_Pt-100to200':   [None,['ZJetsToNuNu_BGenFilter_Pt-100to200']],
    'ZJetsToNuNu_BGenFilter_Pt-200toInf':   [None,['ZJetsToNuNu_BGenFilter_Pt-200toInf']],

    'WW':                              [None, ['WW']],
    'WZ':                              [None, ['WZ']],
    'ZZ':                              [None, ['ZZ']],
}

