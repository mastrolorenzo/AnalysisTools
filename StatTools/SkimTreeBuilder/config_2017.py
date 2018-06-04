## SkimTreeBuilder Configuration File
from new_branches import *

# The input pattern used to glob for sample output files.
input_pattern = '%s'

# The output directory for skimmed trees.
destination = 'samples_2017'

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
    'ZH125_ZNuNu_powheg':                [None, ['/eos/uscms/store/group/lpchbb/VHbbAnalysisNtuples/2017V5_June2/ZH125_ZNuNu_powheg/*.root']],
    'ggZH125_ZNuNu_powheg':              [None, ['/eos/uscms/store/group/lpchbb/VHbbAnalysisNtuples/2017V5_June2/ggZH125_ZNuNu_powheg/*.root']],
    'WminusH125_powheg':                 [None, ['/eos/uscms/store/group/lpchbb/VHbbAnalysisNtuples/2017V5_June2/WminusH125_powheg/*.root']],
    'WplusH125_herwig':                  [None, ['/eos/uscms/store/group/lpchbb/VHbbAnalysisNtuples/2017V5_June2/WplusH125_herwig/*.root']],
    'WplusH125_powheg':                  [None, ['/eos/uscms/store/group/lpchbb/VHbbAnalysisNtuples/2017V5_June2/WplusH125_powheg/*.root']],
    'ZH125_ZLL_powheg':                  [None, ['/eos/uscms/store/group/lpchbb/VHbbAnalysisNtuples/2017V5_June2/ZH125_ZLL_powheg/*.root']],
    'ggZH125_ZLL_powheg':                [None, ['/eos/uscms/store/group/lpchbb/VHbbAnalysisNtuples/2017V5_June2/ggZH125_ZLL_powheg/*.root']],

    'TT_AllHadronic':                    [None, ['/eos/uscms/store/group/lpchbb/VHbbAnalysisNtuples/2017V5_June1_incJES/TT_AllHadronic/*.root']],
    'TT_DiLep':                          [None, ['/eos/uscms/store/group/lpchbb/VHbbAnalysisNtuples/CMSConnect_June2_2017V5/TT_DiLep/*.root']],
    'TT_SingleLep':                      [None, ['/eos/uscms/store/group/lpchbb/VHbbAnalysisNtuples/2017V5_June1_incJES/TT_SingleLep/*.root']],

    'ST_s-c_4f_lep_PSw':                 [None, ['/eos/uscms/store/group/lpchbb/VHbbAnalysisNtuples/2017V5_June1_incJES/ST_s-c_4f_lep_PSw/*.root']],
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
    'WW_1L1Nu2Q':                        [None, ['/eos/uscms/store/group/lpchbb/VHbbAnalysisNtuples/2017V5_June2/WW_1L1Nu2Q/*.root']],
    'WW_LNu2Q_nlo':                      [None, ['/eos/uscms/store/group/lpchbb/VHbbAnalysisNtuples/2017V5_June2/WW_LNu2Q_nlo/*.root']],
    'WZ':                                [None, ['/eos/uscms/store/group/lpchbb/VHbbAnalysisNtuples/2017V5_June2/WZ/*.root']],
    'ZZ':                                [None, ['/eos/uscms/store/group/lpchbb/VHbbAnalysisNtuples/2017V5_June2/ZZ/*.root']],
    'ZZ_4L':                             [None, ['/eos/uscms/store/group/lpchbb/VHbbAnalysisNtuples/2017V5_June2/ZZ_4L/*.root']],

    ###'QCD_HT200to300':                    [None, ['/eos/uscms/store/group/lpchbb/VHbbAnalysisNtuples/2017V5_June1_incJES/QCD_HT200to300/*.root']],
    'QCD_HT300to500':                    [None, ['/eos/uscms/store/group/lpchbb/VHbbAnalysisNtuples/2017V5_June1_incJES/QCD_HT300to500/*.root']],
    'QCD_HT500to700':                    [None, ['/eos/uscms/store/group/lpchbb/VHbbAnalysisNtuples/2017V5_June1_incJES/QCD_HT500to700/*.root']],
    'QCD_HT700to1000':                   [None, ['/eos/uscms/store/group/lpchbb/VHbbAnalysisNtuples/2017V5_June1_incJES/QCD_HT700to1000/*.root']],
    'QCD_HT1000to1500':                  [None, ['/eos/uscms/store/group/lpchbb/VHbbAnalysisNtuples/2017V5_June1_incJES/QCD_HT1000to1500/*.root']],
    'QCD_HT1500to2000':                  [None, ['/eos/uscms/store/group/lpchbb/VHbbAnalysisNtuples/2017V5_June1_incJES/QCD_HT1500to2000/*.root']],
    'QCD_HT2000toInf':                   [None, ['/eos/uscms/store/group/lpchbb/VHbbAnalysisNtuples/2017V5_June1_incJES/QCD_HT2000toInf/*.root']],
}

