the_samples_dict = {
#   name                              index   scale     legend              input tokens (sample name in the previous step)
    'data_obs':                       [0,       1.,     'data_obs',     ['Run2017_MET_MiniAOD', 'Run2017_Ele_ReMiniAOD','Run2017_Mu_ReMiniAOD','Run2017_DoubleEle_ReMiniAOD','Run2017_DoubleMu_ReMiniAOD']],
    'ZH125_Znn_powheg':               [-12504,  1.,     'ZH_hbb',       ['ZH125_ZNuNu_powheg']],
    'ggZH125_Znn_powheg':             [-12505,  1.,     'ggZH_hbb',     ['ggZH125_ZNuNu_powheg']],
    'WH125_WminusLNu_powheg':         [-12501,  1.,     'WH_hbb',       ['WminusH125_powheg']],
    'WH125_WplusLNu_powheg':          [-12500,  1.,     'WH_hbb',       ['WplusH125_powheg']],
    'TT_DiLep':                       [202,     1.,     'TT',           ['TT_DiLep']],
    'TT_SingleLep':                   [203,     1.,     'TT',           ['TT_SingleLep']],
    'TT_AllHadronic':                 [212,     1.,     'TT',           ['TT_AllHadronic']],
    'TToLeptons_s':                   [16,      1.,     's_Top',        ['ST_s-c_4f_lep_PSw']],
    'TToLeptons_t':                   [18,      1.,     's_Top',        ['ST_t-c_top_4f_inc']],
    'TBarToLeptons_t':                [19,      1.,     's_Top',        ['ST_t-c_antitop_4f_inc']],
    'T_tW':                           [20,      1.,     's_Top',        ['ST_tW_top_5f_inc_PSw']],
    'Tbar_tW':                        [21,      1.,     's_Top',        ['ST_tW_antitop_5f_inc']],
    ##'QCDHT100To200':                  [1 ,    1.,       'QCD',          ['QCD_HT100to200']],
    'QCDHT200To300':                  [2 ,      1.,     'QCD',          ['QCD_HT200to300']],
    'QCDHT300To500':                  [3 ,      1.,     'QCD',          ['QCD_HT300to500']],
    'QCDHT500To700':                  [4 ,      1.,     'QCD',          ['QCD_HT500to700']],
    'QCDHT700To1000':                 [5 ,      1.,     'QCD',          ['QCD_HT700to1000']],
    'QCDHT1000To1500':                [6 ,      1.,     'QCD',          ['QCD_HT1000to1500']],
    'QCDHT1500To2000':                [7 ,      1.,     'QCD',          ['QCD_HT1500to2000']],
    'QCDHT2000ToInf':                 [8 ,      1.,     'QCD',          ['QCD_HT2000toInf']],
    'Z_udcsg':                        [15000,   1.,     'Zj0b',         ['ZJetsToNuNu_HT100To200']],
    'Z_b':                            [15001,   1.,     'Zj1b',         ['ZJetsToNuNu_HT100To200']],
    'Z_bb':                           [15002,   1.,     'Zj2b',         ['ZJetsToNuNu_HT100To200']],
    'Z_udcsg_HT200to400':             [15100,   1.,     'Zj0b',         ['ZJetsToNuNu_HT200To400']],
    'Z_b_HT200to400':                 [15101,   1.,     'Zj1b',         ['ZJetsToNuNu_HT200To400']],
    'Z_bb_HT200to400':                [15102,   1.,     'Zj2b',         ['ZJetsToNuNu_HT200To400']],
    'Z_udcsg_HT400to600':             [15200,   1.,     'Zj0b',         ['ZJetsToNuNu_HT400To600']],
    'Z_b_HT400to600':                 [15201,   1.,     'Zj1b',         ['ZJetsToNuNu_HT400To600']],
    'Z_bb_HT400to600':                [15202,   1.,     'Zj2b',         ['ZJetsToNuNu_HT400To600']],
    'Z_udcsg_HT600to800':             [15300,   1.,     'Zj0b',         ['ZJetsToNuNu_HT600To800']],
    'Z_b_HT600to800':                 [15301,   1.,     'Zj1b',         ['ZJetsToNuNu_HT600To800']],
    'Z_bb_HT600to800':                [15302,   1.,     'Zj2b',         ['ZJetsToNuNu_HT600To800']],
    'Z_udcsg_HT800to1200':            [15400,   1.,     'Zj0b',         ['ZJetsToNuNu_HT800To1200']],
    'Z_b_HT800to1200':                [15401,   1.,     'Zj1b',         ['ZJetsToNuNu_HT800To1200']],
    'Z_bb_HT800to1200':               [15402,   1.,     'Zj2b',         ['ZJetsToNuNu_HT800To1200']],
    'Z_udcsg_HT2500toInf':            [15600,   1.,     'Zj0b',         ['ZJetsToNuNu_HT2500ToInf']],
    'Z_b_HT2500toInf':                [15601,   1.,     'Zj1b',         ['ZJetsToNuNu_HT2500ToInf']],
    'Z_bb_HT2500toInf':               [15602,   1.,     'Zj2b',         ['ZJetsToNuNu_HT2500ToInf']],
    'ZBJetsToNuNu_Pt-100to200_udcsg': [16000,   1.,     'Zj0b',         ['ZBJetsToNuNu_Pt-100to200']],
    'ZBJetsToNuNu_Pt-100to200_b':     [16001,   1.,     'Zj1b',         ['ZBJetsToNuNu_Pt-100to200']],
    'ZBJetsToNuNu_Pt-100to200_bb':    [16002,   1.,     'Zj2b',         ['ZBJetsToNuNu_Pt-100to200']],
    'ZBJetsToNuNu_Pt-200toInf_udcsg': [16100,   1.,     'Zj0b',         ['ZBJetsToNuNu_Pt-200toInf']],
    'ZBJetsToNuNu_Pt-200toInf_b':     [16101,   1.,     'Zj1b',         ['ZBJetsToNuNu_Pt-200toInf']],
    'ZBJetsToNuNu_Pt-200toInf_bb':    [16102,   1.,     'Zj2b',         ['ZBJetsToNuNu_Pt-200toInf']],
    'ZJetsToNuNu_BGenFilter_Pt-100to200_udcsg': [16200, 1., 'Zj0b',     ['ZJetsToNuNu_BGenFilter_Pt-100to200']],
    'ZJetsToNuNu_BGenFilter_Pt-100to200_b':     [16201, 1., 'Zj1b',     ['ZJetsToNuNu_BGenFilter_Pt-100to200']],
    'ZJetsToNuNu_BGenFilter_Pt-100to200_bb':    [16202, 1., 'Zj2b',     ['ZJetsToNuNu_BGenFilter_Pt-100to200']],
    'ZJetsToNuNu_BGenFilter_Pt-200toInf_udcsg': [16300, 1., 'Zj0b',     ['ZJetsToNuNu_BGenFilter_Pt-200toInf']],
    'ZJetsToNuNu_BGenFilter_Pt-200toInf_b':     [16301, 1., 'Zj1b',     ['ZJetsToNuNu_BGenFilter_Pt-200toInf']],
    'ZJetsToNuNu_BGenFilter_Pt-200toInf_bb':    [16302, 1., 'Zj2b',     ['ZJetsToNuNu_BGenFilter_Pt-200toInf']],
    #'Z_udcsg_HT1200to2500':          [15500,   1.,     'Zj0b',         ['ZJetsToNuNu_HT1200To2500']],
    #'Z_b_HT1200to2500':              [15501,   1.,     'Zj1b',         ['ZJetsToNuNu_HT1200To2500']],
    #'Z_bb_HT1200to2500':             [15502,   1.,     'Zj2b',         ['ZJetsToNuNu_HT1200To2500']],
    'W_udcsg':                        [4000,    1.,     'Wj0b',         ['WJets_madgraph']],
    'W_b':                            [4001,    1.,     'Wj1b',         ['WJets_madgraph']],
    'W_bb':                           [4002,    1.,     'Wj2b',         ['WJets_madgraph']],
    'W_udcsg_HT100To200':             [4100,    1.,     'Wj0b',         ['WJets-HT100To200']],
    'W_b_HT100To200':                 [4101,    1.,     'Wj1b',         ['WJets-HT100To200']],
    'W_bb_HT100To200':                [4102,    1.,     'Wj2b',         ['WJets-HT100To200']],
    'W_udcsg_HT200To400':             [4200,    1.,     'Wj0b',         ['WJets-HT200To400']],
    'W_b_HT200To400':                 [4201,    1.,     'Wj1b',         ['WJets-HT200To400']],
    'W_bb_HT200To400':                [4202,    1.,     'Wj2b',         ['WJets-HT200To400']],
    'W_udcsg_HT400To600':             [4300,    1.,     'Wj0b',         ['WJets-HT400To600']],
    'W_b_HT400To600':                 [4301,    1.,     'Wj1b',         ['WJets-HT400To600']],
    'W_bb_HT400To600':                [4302,    1.,     'Wj2b',         ['WJets-HT400To600']],
    'W_udcsg_HT600To800':             [4400,    1.,     'Wj0b',         ['WJets-HT600To800']],
    'W_b_HT600To800':                 [4401,    1.,     'Wj1b',         ['WJets-HT600To800']],
    'W_bb_HT600To800':                [4402,    1.,     'Wj2b',         ['WJets-HT600To800']],
    'W_udcsg_HT800To1200':            [4500,    1.,     'Wj0b',         ['WJets-HT800To1200']],
    'W_b_HT800To1200':                [4501,    1.,     'Wj1b',         ['WJets-HT800To1200']],
    'W_bb_HT800To1200':               [4502,    1.,     'Wj2b',         ['WJets-HT800To1200']],
    'W_udcsg_HT1200To2500':           [4600,    1.,     'Wj0b',         ['WJets-HT1200To2500']],
    'W_b_HT1200To2500':               [4601,    1.,     'Wj1b',         ['WJets-HT1200To2500']],
    'W_bb_HT1200To2500':              [4602,    1.,     'Wj2b',         ['WJets-HT1200To2500']],
    'W_udcsg_HT2500ToInf':            [4700,    1.,     'Wj0b',         ['WJets-HT2500ToInf']],
    'W_b_HT2500ToInf':                [4701,    1.,     'Wj1b',         ['WJets-HT2500ToInf']],
    'W_bb_HT2500ToInf':               [4702,    1.,     'Wj2b',         ['WJets-HT2500ToInf']],
    'W_udcsg_b100to200':              [5000,    1.,     'Wj0b',         ['WBJets-Pt100To200']],
    'W_b_b100to200':                  [5001,    1.,     'Wj1b',         ['WBJets-Pt100To200']],
    'W_bb_b100to200':                 [5002,    1.,     'Wj2b',         ['WBJets-Pt100To200']],
    'W_udcsg_b200toInf':              [5100,    1.,     'Wj0b',         ['WBJets-Pt200ToInf']],
    'W_b_b200toInf':                  [5101,    1.,     'Wj1b',         ['WBJets-Pt200ToInf']],
    'W_bb_b200toInf':                 [5102,    1.,     'Wj2b',         ['WBJets-Pt200ToInf']],
    'W_udcsg_bgen100to200':           [5300,    1.,     'Wj0b',         ['WJets_BGenFilter-Pt100To200']],
    'W_b_bgen100to200':               [5301,    1.,     'Wj1b',         ['WJets_BGenFilter-Pt100To200']],
    'W_bb_bgen100to200':              [5302,    1.,     'Wj2b',         ['WJets_BGenFilter-Pt100To200']],
    'W_udcsg_bgen200toInf':           [5400,    1.,     'Wj0b',         ['WJets_BGenFilter-Pt200ToInf']],
    'W_b_bgen200toInf':               [5401,    1.,     'Wj1b',         ['WJets_BGenFilter-Pt200ToInf']],
    'W_bb_bgen200toInf':              [5402,    1.,     'Wj2b',         ['WJets_BGenFilter-Pt200ToInf']],
    'ZZ_lf':                          [3500,    1.,     'VVLF',         ['ZZ']],
    'ZZ_b':                           [3501,    1.,     'VVHF',         ['ZZ']],
    'ZZ_bb':                          [3502,    1.,     'VVHF',         ['ZZ']],
    'WZ_1L1Nu2Q_lf':                  [3200,    1.,     'VVLF',         ['WZ_1L1Nu2Q']],
    'WZ_1L1Nu2Q_b':                   [3201,    1.,     'VVHF',         ['WZ_1L1Nu2Q']],
    'WZ_1L1Nu2Q_bb':                  [3202,    1.,     'VVHF',         ['WZ_1L1Nu2Q']],
    'WW_1L1Nu2Q_lf':                  [3400,    1.,     'VVLF',         ['WW_1L1Nu2Q']],
    'WW_2L1Nu2Q_b':                   [3401,    1.,     'VVHF',         ['WW_1L1Nu2Q']],
    'WW_2L1Nu2Q_bb':                  [3402,    1.,     'VVHF',         ['WW_1L1Nu2Q']],
}
