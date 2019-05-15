# Configuration file for 2018 samples

# The name of the TTree in the ntuple.
treename = 'Events'

# The stacking order of physics processes in the plots.
# The list is ordered from top to bottom of the stack.
stacking_order = [
    'ZH(b#bar{b})',
    'WH(b#bar{b})',
    'ggZH(b#bar{b})',
    'QCD',
    'VV+LF',
    'VV+HF',
    'Z+udcsg',
    'Z+b',
    'Z+b#bar{b}',
    'W+udcsg',
    'W+b',
    'W+b#bar{b}',
    't#bar{t}',
    'Single top',
    #'2017',
]

# The fill and line colors assigned to each physics process
# in the plots. The current palette copies HIG-16-044.
sample_colors = {
    #'2017': 1,
    'ZH(b#bar{b})':   2,
    'ggZH(b#bar{b})': 625,
    'WH(b#bar{b})':   634,
    'QCD':            613,
    'VV+LF':          13,
    'VV+HF':          17,
    'Z+udcsg':        402,
    'Z+b':            397,
    'Z+b#bar{b}':     5,
    'W+udcsg':        814,
    'W+b':            815,
    'W+b#bar{b}':     81,
    't#bar{t}':       4,
    'Single top':     70,
}

# The input pattern used to glob for sample files given their input token.
#input_pattern = '/nfs/dust/cms/user/hkaveh/VHbbAnalysisNtuples/vhbb_2018_1204/haddjobs/pca/sum_%s*_pca.root'
input_pattern = '/nfs/dust/cms/user/hkaveh/VHbbAnalysisNtuples/vhbb_2018_1204/haddjobs/sum_%s*.root'

#input_pattern = '/nfs/dust/cms/user/hkaveh/VHbbAnalysisNtuples/vhbb_2018_1203/*/output_%s*.root'

def get_samples(channel, signal_overlay=False, **kwargs):
    """Return a dicionary of samples for a given channel.

    Parameters
    ----------
    channel : string
        Either "Znn", "Wln", or "Zll".
    signal_overlay : bool, optional
        Flag to add signal overlays. The default is False.
    **kwargs
        The following additional keyword arguments are used to
        assign scale factors to specific processes:
          - sf_tt
          - sf_wj0b
          - sf_wj1b
          - sf_wj2b
          - sf_zj0b
          - sf_zj1b
          - sf_zj2b
    """
    sf_tt = kwargs.get('sf_tt', 1.)
    sf_wj0b = kwargs.get('sf_wj0b', 1.)
    sf_wj1b = kwargs.get('sf_wj1b', 1.)
    sf_wj2b = kwargs.get('sf_wj2b', 1.)
    sf_zj0b = kwargs.get('sf_zj0b', 1.)
    sf_zj1b = kwargs.get('sf_zj1b', 1.)
    sf_zj2b = kwargs.get('sf_zj2b', 1.)
    
    #sf_zj0b = 1
    #sf_zj1b = 1
    #sf_zj2b = 1       
    #sf_lumi = 38588.8/58830.
    
    sf_lumi = 1.
    # Samples common to all channels.

    samples = {
        'ST_s':           [16, 1.*sf_lumi, 'Single top', ['ST_s-c_4f_lep']],
        'ST_t_antitop':   [19, 1.*sf_lumi, 'Single top', ['ST_t-c_antitop_4f_inc']],
        'ST_t_top':       [18, 1.*sf_lumi, 'Single top', ['ST_t-c_top_4f_inc']],
        'ST_tW_antitop':  [21, 1.*sf_lumi, 'Single top', ['ST_tW_antitop_5f_inc']],
        'ST_tW_top':      [20, 1.*sf_lumi, 'Single top', ['ST_tW_top_5f_inc']],

        'TT_AllHadronic': [212, sf_tt*sf_lumi, 't#bar{t}', ['TT_AllHadronic']],
        'TT_DiLep':       [202, sf_tt*sf_lumi, 't#bar{t}', ['TT_DiLep_NoPSWeights']],
        'TT_SingleLep':   [203, sf_tt*sf_lumi, 't#bar{t}', ['TT_SingleLep']],

        'WWLF': ['sampleIndex%100!=2', 1.*sf_lumi, 'VV+LF', ['WW']],
        'WWHF': ['sampleIndex%100==2', 1.*sf_lumi, 'VV+HF', ['WW']],
        'WZLF': ['sampleIndex%100!=2', 1.*sf_lumi, 'VV+LF', ['WZ']],
        'WZHF': ['sampleIndex%100==2', 1.*sf_lumi, 'VV+HF', ['WZ']],
        'ZZLF': ['sampleIndex%100!=2', 1.*sf_lumi, 'VV+LF', ['ZZ']],
        'ZZHF': ['sampleIndex%100==2', 1.*sf_lumi, 'VV+HF', ['ZZ']],
    }

    # Samples specific to the Znn channel.
    if channel == 'Znn':
        samples.update({
            'Data': [0, 1., 'Data', ['Run2017_Ele_ReMiniAOD',
                                     #'Run2017_MET_MiniAOD', FIXME
                                     'Run2017_Mu_ReMiniAOD']],

            'ZnnH_powheg':    [-12504, 1., 'ZH(b#bar{b})',   ['ZH125_ZNuNu_powheg']],
            'ggZnnH_powheg':  [-12505, 1., 'ggZH(b#bar{b})', ['ggZH125_ZNuNu_powheg']],
            'WminusH_powheg': [-12501, 1., 'WH(b#bar{b})',   ['WminusH125_powheg']],
            'WplusH_powheg':  [-12500, 1., 'WH(b#bar{b})',   ['WplusH125_powheg']],

            'ZJetsToNuNu_HT100To200_udcsg':  ['sampleIndex%100==0', sf_zj0b, 'Z+udcsg',    ['ZJetsToNuNu_HT100To200']],
            'ZJetsToNuNu_HT100To200_b':      ['sampleIndex%100==1', sf_zj1b, 'Z+b',        ['ZJetsToNuNu_HT100To200']],
            'ZJetsToNuNu_HT100To200_bb':     ['sampleIndex%100==2', sf_zj2b, 'Z+b#bar{b}', ['ZJetsToNuNu_HT100To200']],
            'ZJetsToNuNu_HT200To400_udcsg':  ['sampleIndex%100==0', sf_zj0b, 'Z+udcsg',    ['ZJetsToNuNu_HT200To400']],
            'ZJetsToNuNu_HT200To400_b':      ['sampleIndex%100==1', sf_zj1b, 'Z+b',        ['ZJetsToNuNu_HT200To400']],
            'ZJetsToNuNu_HT200To400_bb':     ['sampleIndex%100==2', sf_zj2b, 'Z+b#bar{b}', ['ZJetsToNuNu_HT200To400']],
            'ZJetsToNuNu_HT400To600_udcsg':  ['sampleIndex%100==0', sf_zj0b, 'Z+udcsg',    ['ZJetsToNuNu_HT400To600']],
            'ZJetsToNuNu_HT400To600_b':      ['sampleIndex%100==1', sf_zj1b, 'Z+b',        ['ZJetsToNuNu_HT400To600']],
            'ZJetsToNuNu_HT400To600_bb':     ['sampleIndex%100==2', sf_zj2b, 'Z+b#bar{b}', ['ZJetsToNuNu_HT400To600']],
            'ZJetsToNuNu_HT600To800_udcsg':  ['sampleIndex%100==0', sf_zj0b, 'Z+udcsg',    ['ZJetsToNuNu_HT600To800']],
            'ZJetsToNuNu_HT600To800_b':      ['sampleIndex%100==1', sf_zj1b, 'Z+b',        ['ZJetsToNuNu_HT600To800']],
            'ZJetsToNuNu_HT600To800_bb':     ['sampleIndex%100==2', sf_zj2b, 'Z+b#bar{b}', ['ZJetsToNuNu_HT600To800']],
            'ZJetsToNuNu_HT800To1200_udcsg': ['sampleIndex%100==0', sf_zj0b, 'Z+udcsg',    ['ZJetsToNuNu_HT800To1200']],
            'ZJetsToNuNu_HT800To1200_b':     ['sampleIndex%100==1', sf_zj1b, 'Z+b',        ['ZJetsToNuNu_HT800To1200']],
            'ZJetsToNuNu_HT800To1200_bb':    ['sampleIndex%100==2', sf_zj2b, 'Z+b#bar{b}', ['ZJetsToNuNu_HT800To1200']],
        })

        if signal_overlay:
            samples.update({
                'ZnnH125_powheg':    [-12504, 1., 'ZH(b#bar{b})',   ['ZH125_ZNuNu_powheg']],
                'ggZnnH125_powheg':  [-12505, 1., 'ggZH(b#bar{b})', ['ggZH125_ZNuNu_powheg']],
                'WminusH125_powheg': [-12501, 1., 'WH(b#bar{b})',   ['WminusH125_powheg']],
                'WplusH125_powheg':  [-12500, 1., 'WH(b#bar{b})',   ['WplusH125_powheg']],
            })

    # Samples specific to the Wln channel.
    if channel == 'Wln':
        samples.update({
            'Data': [0, 1., 'Data', ['Run2018_Ele',
                                     'Run2018_Mu']],

            'WminusH_powheg':   [-12501, 1., 'WH(b#bar{b})', ['WminusH125_powheg']],
            'WplusH_powheg':    [-12500, 1., 'WH(b#bar{b})', ['WplusH125_powheg']],
        })

        if signal_overlay:
            samples.update({
                'WminusH125_powheg': [-12501, 1., 'WH(b#bar{b})', ['WminusH125_powheg']],
                'WplusH125_powheg':  [-12500, 1., 'WH(b#bar{b})', ['WplusH125_powheg']],
            })

    # Samples shared by the Znn and Wln channels.
    if channel == 'Znn' or channel == 'Wln':
        samples.update({
            'WJets_madgraph_udcsg':     ['sampleIndex%100==0', sf_wj0b*sf_lumi, 'W+udcsg',    ['WJets_madgraph']],
            'WJets_madgraph_b':         ['sampleIndex%100==1', sf_wj1b*sf_lumi, 'W+b',        ['WJets_madgraph']],
            'WJets_madgraph_bb':        ['sampleIndex%100==2', sf_wj2b*sf_lumi, 'W+b#bar{b}', ['WJets_madgraph']],
            'WJets-HT100To200_udcsg':   ['sampleIndex%100==0', sf_wj0b*sf_lumi, 'W+udcsg',    ['WJets-HT100To200']],
            'WJets-HT100To200_b':       ['sampleIndex%100==1', sf_wj1b*sf_lumi, 'W+b',        ['WJets-HT100To200']],
            'WJets-HT100To200_bb':      ['sampleIndex%100==2', sf_wj2b*sf_lumi, 'W+b#bar{b}', ['WJets-HT100To200']],
            'WJets-HT200To400_udcsg':   ['sampleIndex%100==0', sf_wj0b*sf_lumi, 'W+udcsg',    ['WJets-HT200To400']],
            'WJets-HT200To400_b':       ['sampleIndex%100==1', sf_wj1b*sf_lumi, 'W+b',        ['WJets-HT200To400']],
            'WJets-HT200To400_bb':      ['sampleIndex%100==2', sf_wj2b*sf_lumi, 'W+b#bar{b}', ['WJets-HT200To400']],
            'WJets-HT400To600_udcsg':   ['sampleIndex%100==0', sf_wj0b*sf_lumi, 'W+udcsg',    ['WJets-HT400To600']],
            'WJets-HT400To600_b':       ['sampleIndex%100==1', sf_wj1b*sf_lumi, 'W+b',        ['WJets-HT400To600']],
            'WJets-HT400To600_bb':      ['sampleIndex%100==2', sf_wj2b*sf_lumi, 'W+b#bar{b}', ['WJets-HT400To600']],
            'WJets-HT600To800_udcsg':   ['sampleIndex%100==0', sf_wj0b*sf_lumi, 'W+udcsg',    ['WJets-HT600To800']],
            'WJets-HT600To800_b':       ['sampleIndex%100==1', sf_wj1b*sf_lumi, 'W+b',        ['WJets-HT600To800']],
            'WJets-HT600To800_bb':      ['sampleIndex%100==2', sf_wj2b*sf_lumi, 'W+b#bar{b}', ['WJets-HT600To800']],
            'WJets-HT800To1200_udcsg':  ['sampleIndex%100==0', sf_wj0b*sf_lumi, 'W+udcsg',    ['WJets-HT800To1200']],
            'WJets-HT800To1200_b':      ['sampleIndex%100==1', sf_wj1b*sf_lumi, 'W+b',        ['WJets-HT800To1200']],
            'WJets-HT800To1200_bb':     ['sampleIndex%100==2', sf_wj2b*sf_lumi, 'W+b#bar{b}', ['WJets-HT800To1200']],
            'WJets-HT1200To2500_udcsg': ['sampleIndex%100==0', sf_wj0b*sf_lumi, 'W+udcsg',    ['WJets-HT1200To2500']],
            'WJets-HT1200To2500_b':     ['sampleIndex%100==1', sf_wj1b*sf_lumi, 'W+b',        ['WJets-HT1200To2500']],
            'WJets-HT1200To2500_bb':    ['sampleIndex%100==2', sf_wj2b*sf_lumi, 'W+b#bar{b}', ['WJets-HT1200To2500']],

            #'WJets-HT100To200_udcsg':   ['sampleIndex%100==0', sf_wj0b, 'W+udcsg',    ['WJets_BGenFilter-Pt100To200']],
            #'WJets-HT100To200_b':       ['sampleIndex%100==1', sf_wj1b, 'W+b',        ['WJets_BGenFilter-Pt100To200']],
            #'WJets-HT100To200_bb':      ['sampleIndex%100==2', sf_wj2b, 'W+b#bar{b}', ['WJets_BGenFilter-Pt100To200']],
            #'WJets-HT200To400_udcsg':   ['sampleIndex%100==0', sf_wj0b, 'W+udcsg',    ['WJets_BGenFilter-Pt200To400']],
            #'WJets-HT200To400_b':       ['sampleIndex%100==1', sf_wj1b, 'W+b',        ['WJets_BGenFilter-Pt200To400']],
            #'WJets-HT200To400_bb':      ['sampleIndex%100==2', sf_wj2b, 'W+b#bar{b}', ['WJets_BGenFilter-Pt200To400']],
            
            'QCD_HT200to300':   [2, 1., 'QCD', ['QCD_HT200to300']],
            'QCD_HT300to500':   [3, 1., 'QCD', ['QCD_HT300to500']],
            'QCD_HT500to700':   [4, 1., 'QCD', ['QCD_HT500to700']],
            'QCD_HT700to1000':  [5, 1., 'QCD', ['QCD_HT700to1000']],
            'QCD_HT1000to1500': [6, 1., 'QCD', ['QCD_HT1000to1500']],
            'QCD_HT1500to2000': [7, 1., 'QCD', ['QCD_HT1500to2000']],
            'QCD_HT2000toInf':  [8, 1., 'QCD', ['QCD_HT2000toInf']],
            
            'DYToLL_madgraph_udcsg':              ['sampleIndex%100==0', sf_zj0b*sf_lumi, 'Z+udcsg',    ['DYToLL_madgraph']],
            'DYToLL_madgraph_b':                  ['sampleIndex%100==1', sf_zj1b*sf_lumi, 'Z+b',        ['DYToLL_madgraph']],
            'DYToLL_madgraph_bb':                 ['sampleIndex%100==2', sf_zj2b*sf_lumi, 'Z+b#bar{b}', ['DYToLL_madgraph']],
            'DYToLL_HT100to200_madgraph_udcsg':   ['sampleIndex%100==0', sf_zj0b*sf_lumi, 'Z+udcsg',    ['DYToLL_HT100to200_madgraph']],
            'DYToLL_HT100to200_madgraph_b':       ['sampleIndex%100==1', sf_zj1b*sf_lumi, 'Z+b',        ['DYToLL_HT100to200_madgraph']],
            'DYToLL_HT100to200_madgraph_bb':      ['sampleIndex%100==2', sf_zj2b*sf_lumi, 'Z+b#bar{b}', ['DYToLL_HT100to200_madgraph']],
           # 'DYToLL_HT200to400_madgraph_udcsg':   ['sampleIndex%100==0', sf_zj0b, 'Z+udcsg',    ['DYToLL_HT200to400_madgraph']],
           # 'DYToLL_HT200to400_madgraph_b':       ['sampleIndex%100==1', sf_zj1b, 'Z+b',        ['DYToLL_HT200to400_madgraph']],
           # 'DYToLL_HT200to400_madgraph_bb':      ['sampleIndex%100==2', sf_zj2b, 'Z+b#bar{b}', ['DYToLL_HT200to400_madgraph']],
            'DYToLL_HT400to600_madgraph_udcsg':   ['sampleIndex%100==0', sf_zj0b*sf_lumi, 'Z+udcsg',    ['DYToLL_HT400to600_madgraph']],
            'DYToLL_HT400to600_madgraph_b':       ['sampleIndex%100==1', sf_zj1b*sf_lumi, 'Z+b',        ['DYToLL_HT400to600_madgraph']],
            'DYToLL_HT400to600_madgraph_bb':      ['sampleIndex%100==2', sf_zj2b*sf_lumi, 'Z+b#bar{b}', ['DYToLL_HT400to600_madgraph']],
            'DYToLL_HT600to800_madgraph_udcsg':   ['sampleIndex%100==0', sf_zj0b*sf_lumi, 'Z+udcsg',    ['DYToLL_HT600to800_madgraph']],
            'DYToLL_HT600to800_madgraph_b':       ['sampleIndex%100==1', sf_zj1b*sf_lumi, 'Z+b',        ['DYToLL_HT600to800_madgraph']],
            'DYToLL_HT600to800_madgraph_bb':      ['sampleIndex%100==2', sf_zj2b*sf_lumi, 'Z+b#bar{b}', ['DYToLL_HT600to800_madgraph']],
            'DYToLL_HT800to1200_madgraph_udcsg':  ['sampleIndex%100==0', sf_zj0b*sf_lumi, 'Z+udcsg',    ['DYToLL_HT800to1200_madgraph']],
            'DYToLL_HT800to1200_madgraph_b':      ['sampleIndex%100==1', sf_zj1b*sf_lumi, 'Z+b',        ['DYToLL_HT800to1200_madgraph']],
            'DYToLL_HT800to1200_madgraph_bb':     ['sampleIndex%100==2', sf_zj2b*sf_lumi, 'Z+b#bar{b}', ['DYToLL_HT800to1200_madgraph']],
          
        })

    # Samples specific to the Zll channel.

    if channel == 'Zll':
        samples.update({
            'Data': [0, 1., 'Data', ['Run2018_DoubleEle',
                                     'Run2018_DoubleMu']],
            
      #      'Data': [0, (58830./41530.), 'Data', ['Run2017_DoubleEle_ReMiniAOD',
      #                                            'Run2017_DoubleMu_ReMiniAOD']],
            'ZH_ZLL_powheg':   [-12502, 1., 'ZH(b#bar{b})',   ['ZH125_ZLL_powheg']],
            'ggZH_ZLL_powheg': [-12503, 1., 'ggZH(b#bar{b})', ['ggZH125_ZLL_powheg']],
            'DYToLL_madgraph_udcsg':              ['sampleIndex%100==0', sf_zj0b*sf_lumi, 'Z+udcsg',    ['DYToLL_madgraph']],
            'DYToLL_madgraph_b':                  ['sampleIndex%100==1', sf_zj1b*sf_lumi, 'Z+b',        ['DYToLL_madgraph']],
            'DYToLL_madgraph_bb':                 ['sampleIndex%100==2', sf_zj2b*sf_lumi, 'Z+b#bar{b}', ['DYToLL_madgraph']],
            'DYToLL_HT100to200_madgraph_udcsg':   ['sampleIndex%100==0', sf_zj0b*sf_lumi, 'Z+udcsg',    ['DYToLL_HT100to200_madgraph']],
            'DYToLL_HT100to200_madgraph_b':       ['sampleIndex%100==1', sf_zj1b*sf_lumi, 'Z+b',        ['DYToLL_HT100to200_madgraph']],
            'DYToLL_HT100to200_madgraph_bb':      ['sampleIndex%100==2', sf_zj2b*sf_lumi, 'Z+b#bar{b}', ['DYToLL_HT100to200_madgraph']],
           # 'DYToLL_HT200to400_madgraph_udcsg':   ['sampleIndex%100==0', sf_zj0b, 'Z+udcsg',    ['DYToLL_HT200to400_madgraph']],
           # 'DYToLL_HT200to400_madgraph_b':       ['sampleIndex%100==1', sf_zj1b, 'Z+b',        ['DYToLL_HT200to400_madgraph']],
           # 'DYToLL_HT200to400_madgraph_bb':      ['sampleIndex%100==2', sf_zj2b, 'Z+b#bar{b}', ['DYToLL_HT200to400_madgraph']],
            'DYToLL_HT400to600_madgraph_udcsg':   ['sampleIndex%100==0', sf_zj0b*sf_lumi, 'Z+udcsg',    ['DYToLL_HT400to600_madgraph']],
            'DYToLL_HT400to600_madgraph_b':       ['sampleIndex%100==1', sf_zj1b*sf_lumi, 'Z+b',        ['DYToLL_HT400to600_madgraph']],
            'DYToLL_HT400to600_madgraph_bb':      ['sampleIndex%100==2', sf_zj2b*sf_lumi, 'Z+b#bar{b}', ['DYToLL_HT400to600_madgraph']],
            'DYToLL_HT600to800_madgraph_udcsg':   ['sampleIndex%100==0', sf_zj0b*sf_lumi, 'Z+udcsg',    ['DYToLL_HT600to800_madgraph']],
            'DYToLL_HT600to800_madgraph_b':       ['sampleIndex%100==1', sf_zj1b*sf_lumi, 'Z+b',        ['DYToLL_HT600to800_madgraph']],
            'DYToLL_HT600to800_madgraph_bb':      ['sampleIndex%100==2', sf_zj2b*sf_lumi, 'Z+b#bar{b}', ['DYToLL_HT600to800_madgraph']],
            'DYToLL_HT800to1200_madgraph_udcsg':  ['sampleIndex%100==0', sf_zj0b*sf_lumi, 'Z+udcsg',    ['DYToLL_HT800to1200_madgraph']],
            'DYToLL_HT800to1200_madgraph_b':      ['sampleIndex%100==1', sf_zj1b*sf_lumi, 'Z+b',        ['DYToLL_HT800to1200_madgraph']],
            'DYToLL_HT800to1200_madgraph_bb':     ['sampleIndex%100==2', sf_zj2b*sf_lumi, 'Z+b#bar{b}', ['DYToLL_HT800to1200_madgraph']],
           # 'DYToLL_HT1200to2500_madgraph_udcsg': ['sampleIndex%100==0', sf_zj0b, 'Z+udcsg',    ['DYToLL_HT1200to2500_madgraph']],
           # 'DYToLL_HT1200to2500_madgraph_b':     ['sampleIndex%100==1', sf_zj1b, 'Z+b',        ['DYToLL_HT1200to2500_madgraph']],
           # 'DYToLL_HT1200to2500_madgraph_bb':    ['sampleIndex%100==2', sf_zj2b, 'Z+b#bar{b}', ['DYToLL_HT1200to2500_madgraph']],
           # 'DYToLL_HT2500toInf_madgraph_udcsg':  ['sampleIndex%100==0', sf_zj0b, 'Z+udcsg',    ['DYToLL_HT2500toInf_madgraph']],
           # 'DYToLL_HT2500toInf_madgraph_b':      ['sampleIndex%100==1', sf_zj1b, 'Z+b',        ['DYToLL_HT2500toInf_madgraph']],
           # 'DYToLL_HT2500toInf_madgraph_bb':     ['sampleIndex%100==2', sf_zj2b, 'Z+b#bar{b}', ['DYToLL_HT2500toInf_madgraph']],

           # 'DYToLL_M4to50_HT70to100_madgraph':  ['sampleIndex%100==0', sf_zj0b, 'Z+udcsg',    ['DYToLL_M4to50_HT70to100_madgraph']],
           # 'DYToLL_M4to50_HT70to100_madgraph':  ['sampleIndex%100==1', sf_zj1b, 'Z+b',        ['DYToLL_M4to50_HT70to100_madgraph']],
           # 'DYToLL_M4to50_HT70to100_madgraph':  ['sampleIndex%100==2', sf_zj2b, 'Z+b#bar{b}', ['DYToLL_M4to50_HT70to100_madgraph']],
           # 'DYToLL_M4to50_HT100to200_madgraph': ['sampleIndex%100==0', sf_zj0b, 'Z+udcsg',    ['DYToLL_M4to50_HT100to200_madgraph']],
           # 'DYToLL_M4to50_HT100to200_madgraph': ['sampleIndex%100==1', sf_zj1b, 'Z+b',        ['DYToLL_M4to50_HT100to200_madgraph']],
           # 'DYToLL_M4to50_HT100to200_madgraph': ['sampleIndex%100==2', sf_zj2b, 'Z+b#bar{b}', ['DYToLL_M4to50_HT100to200_madgraph']],
           # 'DYToLL_M4to50_HT200to400_madgraph': ['sampleIndex%100==0', sf_zj0b, 'Z+udcsg',    ['DYToLL_M4to50_HT200to400_madgraph']],
           # 'DYToLL_M4to50_HT200to400_madgraph': ['sampleIndex%100==1', sf_zj1b, 'Z+b',        ['DYToLL_M4to50_HT200to400_madgraph']],
           # 'DYToLL_M4to50_HT200to400_madgraph': ['sampleIndex%100==2', sf_zj2b, 'Z+b#bar{b}', ['DYToLL_M4to50_HT200to400_madgraph']],
           # 'DYToLL_M4to50_HT400to600_madgraph': ['sampleIndex%100==0', sf_zj0b, 'Z+udcsg',    ['DYToLL_M4to50_HT400to600_madgraph']],
           # 'DYToLL_M4to50_HT400to600_madgraph': ['sampleIndex%100==1', sf_zj1b, 'Z+b',        ['DYToLL_M4to50_HT400to600_madgraph']],
           # 'DYToLL_M4to50_HT400to600_madgraph': ['sampleIndex%100==2', sf_zj2b, 'Z+b#bar{b}', ['DYToLL_M4to50_HT400to600_madgraph']],
           # 'DYToLL_M4to50_HT600toInf_madgraph': ['sampleIndex%100==0', sf_zj0b, 'Z+udcsg',    ['DYToLL_M4to50_HT600toInf_madgraph']],
           # 'DYToLL_M4to50_HT600toInf_madgraph': ['sampleIndex%100==1', sf_zj1b, 'Z+b',        ['DYToLL_M4to50_HT600toInf_madgraph']],
           # 'DYToLL_M4to50_HT600toInf_madgraph': ['sampleIndex%100==2', sf_zj2b, 'Z+b#bar{b}', ['DYToLL_M4to50_HT600toInf_madgraph']],
        })

        if signal_overlay:
            samples.update({
                #'ZH125_ZLL_powheg':   [-12502, 1.*sf_lumi, 'ZH(b#bar{b})',   ['ZH125_ZLL_powheg']],
                
                #'d_2017': [0, (58830./41530.), '2017', ['Run2017_DoubleMu_ReMiniAOD']],


                #'ggZH125_ZLL_powheg': [-12503, 1., 'ggZH(b#bar{b})', ['ggZH125_ZLL_powheg']],
            })

    # Direct Application of Wln pT(W) Reweighting
    if channel == 'Wln':
        for sample, options in samples.iteritems():
            if sample.startswith('W') and sample.endswith('_udcsg'):
                options[0] = '({0})*(1.115 - 6.36789e-04*V_pt)'.format(options[0])
            elif sample.startswith('W') and sample.endswith(('_b', '_bb')):
                options[0] = '({0})*(1.337 - 1.56131e-03*V_pt)'.format(options[0])
            elif sample.startswith('ST'):
                options[0] = '(sampleIndex=={0})*(1.337 - 1.56131e-03*V_pt)'.format(options[0])

    return samples

