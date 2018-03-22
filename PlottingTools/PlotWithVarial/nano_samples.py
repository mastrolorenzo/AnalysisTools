the_samples_dict = {
#   name                              index   scale     legend              input tokens (sample name in the previous step)
    'Data':                          [0,      1.,       'Data',             ['Run2017_Ele_ReMiniAOD','Run2017_Mu_ReMiniAOD','Run2017_MET_MiniAOD','Run2017_DoubleEle_ReMiniAOD','Run2017_DoubleMu_ReMiniAOD']],
    'TT_DiLep':                      [50,     1.,       't#bar{t}@DL',      ['TT_DiLep']],
    'TT_SingleLep':                  [51,     1.,       't#bar{t}@SL',      ['TT_SingleLep']],
    'TT_AllHadronic':                [52,     1.,       't#bar{t}@Had',     ['TT_AllHadronic']],
    'TToLeptons_s':                  [16,     1.,       'Top',              ['ST_s-c_4f_lep']],
    'TToLeptons_t':                  [17,     1.,       'Top',              ['ST_t-c_top_4f_inc']],
    'TBarToLeptons_t':               [18,     1.,       'Top',              ['ST_t-c_antitop_4f_inc']],
    'T_tW':                          [20,     1.,       'Top',              ['ST_tW_top_5f_noFullHad']],
    'Tbar_tW':                       [21,     1.,       'Top',              ['ST_tW_antitop_5f_noFullHad']],
    'Z_udcsg':                       [2300,   1./191860448,       'Z+udcsg',          ['DYToLL_madgraph']],
    'Z_b':                           [2301,   1./191860448,       'Z+b',              ['DYToLL_madgraph']],
    'Z_bb':                          [2302,   1./191860448,       'Z+bb',             ['DYToLL_madgraph']],
    'Z_udcsg_HT100To200':            [6100,   1.,       'Z+udcsg',          ['DYToLL_HT100to200_madgraph']],
    'Z_b_HT100To200':                [6101,   1.,       'Z+b',              ['DYToLL_HT100to200_madgraph']],
    'Z_bb_HT100To200':               [6102,   1.,       'Z+bb',             ['DYToLL_HT100to200_madgraph']],
    'Z_udcsg_HT200To400':            [6200,   1.,       'Z+udcsg',          ['DYToLL_HT200to400_madgraph']],
    'Z_b_HT200To400':                [6201,   1.,       'Z+b',              ['DYToLL_HT200to400_madgraph']],
    'Z_bb_HT200To400':               [6202,   1.,       'Z+bb',             ['DYToLL_HT200to400_madgraph']],
    'Z_udcsg_HT400To600':            [6300,   1.,       'Z+udcsg',          ['DYToLL_HT400to600_madgraph']],
    'Z_b_HT400To600':                [6301,   1.,       'Z+b',              ['DYToLL_HT400to600_madgraph']],
    'Z_bb_HT400To600':               [6302,   1.,       'Z+bb',             ['DYToLL_HT400to600_madgraph']],
    'Z_udcsg_HT600To800':            [6400,   1.,       'Z+udcsg',          ['DYToLL_HT600to800_madgraph']],
    'Z_b_HT600To800':                [6401,   1.,       'Z+b',              ['DYToLL_HT600to800_madgraph']],
    'Z_bb_HT600To800':               [6402,   1.,       'Z+bb',             ['DYToLL_HT600to800_madgraph']],
    'Z_udcsg_HT800To1200':           [6500,   1.,       'Z+udcsg',          ['DYToLL_HT800to1200_madgraph']],
    'Z_b_HT800To1200':               [6501,   1.,       'Z+b',              ['DYToLL_HT800to1200_madgraph']],
    'Z_bb_HT800To1200':              [6502,   1.,       'Z+bb',             ['DYToLL_HT800to1200_madgraph']],
    'Z_udcsg_HT1200To2500':          [6600,   1.,       'Z+udcsg',          ['DYToLL_HT1200to2500_madgraph']],
    'Z_b_HT1200To2500':              [6601,   1.,       'Z+b',              ['DYToLL_HT1200to2500_madgraph']],
    'Z_bb_HT1200To2500':             [6602,   1.,       'Z+bb',             ['DYToLL_HT1200to2500_madgraph']],
    'Z_udcsg_HT2500ToInf':           [6700,   1.,       'Z+udcsg',          ['DYToLL_HT2500toInf_madgraph']],
    'Z_b_HT2500ToInf':               [6701,   1.,       'Z+b',              ['DYToLL_HT2500toInf_madgraph']],
    'Z_bb_HT2500ToInf':              [6702,   1.,       'Z+bb',             ['DYToLL_HT2500toInf_madgraph']],
    'Z_udcsg_LowMHT100To200':        [6900,   1.,       'Z+udcsg',          ['DYToLL_M4to50_HT100to200_madgraph']],
    'Z_b_LowMHT100To200':            [6901,   1.,       'Z+b',              ['DYToLL_M4to50_HT100to200_madgraph']],
    'Z_bb_LowMHT100To200':           [6902,   1.,       'Z+bb',             ['DYToLL_M4to50_HT100to200_madgraph']],
}

sample_colors = {
    'WH125':            880,
    'ZH125':            808,
    't#bar{t}@POWHEG':  600,
    'Top':              432,
    'VV+LF':            922,
    'VV(bb)':           632,
    'W+udcsg':          418,
    'W+b':              410,
    'W+bb':             416,
    'Z+udcsg':          394,
    'Z+b':              398,
    'Z+bb':             400,
}
