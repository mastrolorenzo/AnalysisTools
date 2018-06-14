## SkimTreeBuilder Configuration File
from new_branches import *

# The input pattern used to glob for sample output files.
input_pattern = '%s'

# The output directory for skimmed trees.
destination = 'skims_for_training_2017'

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
    'nAddJets252p5_puid',
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
    'ZH125_ZNuNu_powheg':                [None, ['/eos/uscms/store/group/lpchbb/VHbbAnalysisNtuples/2017V5_June2/ZH125_ZNuNu_powheg/*.root']],
    'ggZH125_ZNuNu_powheg':              [None, ['/eos/uscms/store/group/lpchbb/VHbbAnalysisNtuples/2017V5_June2/ggZH125_ZNuNu_powheg/*.root']],
    'WminusH125_powheg':                 [None, ['/eos/uscms/store/group/lpchbb/VHbbAnalysisNtuples/2017V5_June2/WminusH125_powheg/*.root']],
    'WplusH125_powheg':                  [None, ['/eos/uscms/store/group/lpchbb/VHbbAnalysisNtuples/2017V5_June2/WplusH125_powheg/*.root']],
    'ZH125_ZLL_powheg':                  [None, ['/eos/uscms/store/group/lpchbb/VHbbAnalysisNtuples/2017V5_June2/ZH125_ZLL_powheg/*.root']],
    'ggZH125_ZLL_powheg':                [None, ['/eos/uscms/store/group/lpchbb/VHbbAnalysisNtuples/2017V5_June2/ggZH125_ZLL_powheg/*.root']],

    'TT_AllHadronic':                    [None, ['/eos/uscms/store/group/lpchbb/VHbbAnalysisNtuples/2017V5_June1_incJES/TT_AllHadronic/*.root']],
    'TT_DiLep':                          [None, ['/eos/uscms/store/group/lpchbb/VHbbAnalysisNtuples/CMSConnect_June2_2017V5/TT_DiLep/*.root']],
    'TT_SingleLep':                      [None, ['/eos/uscms/store/group/lpchbb/VHbbAnalysisNtuples/2017V5_June1_incJES/TT_SingleLep/*.root']],

    'ST_s-c_4f_lep_PSw':                 [None, ['/eos/uscms/store/group/lpchbb/VHbbAnalysisNtuples/2017V5_June1_incJES/ST_s-c_4f_lep_PSw/*.root']],
    'ST_t-c_antitop_4f_inc':             [None, ['/eos/uscms/store/group/lpchbb/VHbbAnalysisNtuples/2017V5_June2/ST_t-c_antitop_4f_inc/*.root']],
    'ST_t-c_top_4f_inc':                 [None, ['/eos/uscms/store/group/lpchbb/VHbbAnalysisNtuples/2017V5_June1_incJES/ST_t-c_top_4f_inc/*.root']],
    'ST_tW_antitop_5f_inc':              [None, ['/eos/uscms/store/group/lpchbb/VHbbAnalysisNtuples/2017V5_June1_incJES/ST_tW_antitop_5f_inc/*.root']],
    'ST_tW_top_5f_inc_PSw':              [None, ['/eos/uscms/store/group/lpchbb/VHbbAnalysisNtuples/2017V5_June1_incJES/ST_tW_top_5f_inc_PSw/*.root']],

    'DYToLL_madgraph':                   [None, ['/eos/uscms/store/group/lpchbb/VHbbAnalysisNtuples/2017V5_June1_incJES/DYToLL_madgraph/*.root']],
    'DYToLL_HT100to200_madgraph':        [None, ['/eos/uscms/store/group/lpchbb/VHbbAnalysisNtuples/2017V5_June1_incJES/DYToLL_HT100to200_madgraph/*.root']],
    'DYToLL_HT200to400_madgraph':        [None, ['/eos/uscms/store/group/lpchbb/VHbbAnalysisNtuples/2017V5_June1_incJES/DYToLL_HT200to400_madgraph/*.root']],
    'DYToLL_HT400to600_madgraph':        [None, ['/eos/uscms/store/group/lpchbb/VHbbAnalysisNtuples/2017V5_June1_incJES/DYToLL_HT400to600_madgraph/*.root']],
    'DYToLL_HT600to800_madgraph':        [None, ['/eos/uscms/store/group/lpchbb/VHbbAnalysisNtuples/2017V5_June1_incJES/DYToLL_HT600to800_madgraph/*.root']],
    'DYToLL_HT800to1200_madgraph':       [None, ['/eos/uscms/store/group/lpchbb/VHbbAnalysisNtuples/2017V5_June1_incJES/DYToLL_HT800to1200_madgraph/*.root']],
    'DYToLL_HT1200to2500_madgraph':      [None, ['/eos/uscms/store/group/lpchbb/VHbbAnalysisNtuples/2017V5_June1_incJES/DYToLL_HT1200to2500_madgraph/*.root']],
    'DYToLL_HT2500toInf_madgraph':       [None, ['/eos/uscms/store/group/lpchbb/VHbbAnalysisNtuples/2017V5_June1_incJES/DYToLL_HT2500toInf_madgraph/*.root']],
    'DYToLL_M4to50_HT70to100_madgraph':  [None, ['/eos/uscms/store/group/lpchbb/VHbbAnalysisNtuples/2017V5_June1_incJES/DYToLL_M4to50_HT70to100_madgraph/*.root']],
    'DYToLL_M4to50_HT100to200_madgraph': [None, ['/eos/uscms/store/group/lpchbb/VHbbAnalysisNtuples/2017V5_June1_incJES/DYToLL_M4to50_HT100to200_madgraph/*.root']],
    'DYToLL_M4to50_HT200to400_madgraph': [None, ['/eos/uscms/store/group/lpchbb/VHbbAnalysisNtuples/2017V5_June1_incJES/DYToLL_M4to50_HT200to400_madgraph/*.root']],
    'DYToLL_M4to50_HT400to600_madgraph': [None, ['/eos/uscms/store/group/lpchbb/VHbbAnalysisNtuples/2017V5_June1_incJES/DYToLL_M4to50_HT400to600_madgraph/*.root']],
    'DYToLL_M4to50_HT600toInf_madgraph': [None, ['/eos/uscms/store/group/lpchbb/VHbbAnalysisNtuples/2017V5_June1_incJES/DYToLL_M4to50_HT600toInf_madgraph/*.root']],

    'WJets_madgraph':                    [None, ['/eos/uscms/store/group/lpchbb/VHbbAnalysisNtuples/2017V5_June1_incJES/WJets_madgraph/*.root']],
    'WJets-HT100To200':                  [None, ['/eos/uscms/store/group/lpchbb/VHbbAnalysisNtuples/2017V5_June1_incJES/WJets-HT100To200/*.root']],
    'WJets-HT200To400':                  [None, ['/eos/uscms/store/group/lpchbb/VHbbAnalysisNtuples/2017V5_June1_incJES/WJets-HT200To400/*.root']],
    'WJets-HT400To600':                  [None, ['/eos/uscms/store/group/lpchbb/VHbbAnalysisNtuples/2017V5_June1_incJES/WJets-HT400To600/*.root']],
    'WJets-HT600To800':                  [None, ['/eos/uscms/store/group/lpchbb/VHbbAnalysisNtuples/2017V5_June1_incJES/WJets-HT600To800/*.root']],
    'WJets-HT800To1200':                 [None, ['/eos/uscms/store/group/lpchbb/VHbbAnalysisNtuples/2017V5_June1_incJES/WJets-HT800To1200/*.root']],
    'WJets-HT1200To2500':                [None, ['/eos/uscms/store/group/lpchbb/VHbbAnalysisNtuples/2017V5_June1_incJES/WJets-HT1200To2500/*.root']],

    'ZJetsToNuNu_HT100To200':            [None, ['/eos/uscms/store/group/lpchbb/VHbbAnalysisNtuples/2017V5_June1_incJES/ZJetsToNuNu_HT100To200/*.root']],
    'ZJetsToNuNu_HT200To400':            [None, ['/eos/uscms/store/group/lpchbb/VHbbAnalysisNtuples/2017V5_June1_incJES/ZJetsToNuNu_HT200To400/*.root']],
    'ZJetsToNuNu_HT400To600':            [None, ['/eos/uscms/store/group/lpchbb/VHbbAnalysisNtuples/2017V5_June1_incJES/ZJetsToNuNu_HT400To600/*.root']],
    'ZJetsToNuNu_HT600To800':            [None, ['/eos/uscms/store/group/lpchbb/VHbbAnalysisNtuples/2017V5_June1_incJES/ZJetsToNuNu_HT600To800/*.root']],
    'ZJetsToNuNu_HT800To1200':           [None, ['/eos/uscms/store/group/lpchbb/VHbbAnalysisNtuples/2017V5_June1_incJES/ZJetsToNuNu_HT800To1200/*.root']],

    'WW':                                [None, ['/eos/uscms/store/group/lpchbb/VHbbAnalysisNtuples/2017V5_June2/WW/*.root']],
    'WZ':                                [None, ['/eos/uscms/store/group/lpchbb/VHbbAnalysisNtuples/2017V5_June2/WZ/*.root']],
    'ZZ':                                [None, ['/eos/uscms/store/group/lpchbb/VHbbAnalysisNtuples/2017V5_June2/ZZ/*.root']],
}

