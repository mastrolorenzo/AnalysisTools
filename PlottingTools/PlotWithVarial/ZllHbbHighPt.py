from samples_2017 import *
import varial


####################
# General Settings #
####################

#varial.settings.max_num_processes = 8
varial.settings.rootfile_postfixes += ['.pdf', '.png']
varial.settings.stacking_order = stacking_order
# try to find output on disk and don't run a step if present
enable_reuse_step = True


#################
# Plot Settings #
#################

name = 'plots_ZllHbbHighPt'

weight = 'weight'

plot_vars = {
    'H_pt_fit_fallback':   ('H_pt_fit_fallback',              ';p_{T}(jj) [GeV];',       35, 0,  350),
    'H_mass':              ('H_mass',                         ';M_{jj} [GeV];',          25, 0,  250),
    'H_mass_fit_fallback': ('H_mass_fit_fallback',            ';fitted M_{jj} [GeV];',   25, 0,  250),
    'MET_Pt':              ('MET_Pt',                         ';MET [GeV];',             25, 0,  250),
    'DeepCSV_max':         ('Jet_btagDeepB[hJetInd1]',        ';DeepCSV_{max};',         20, 0,  1),
    'DeepCSV_min':         ('Jet_btagDeepB[hJetInd2]',        ';DeepCSV_{min};',         20, 0,  1),
    'CMVA_max':            ('Jet_btagCMVA[hJetInd1]',         ';CMVA_{max};',            20, -1, 1),
    'CMVA_min':            ('Jet_btagCMVA[hJetInd2]',         ';CMVA_{min};',            20, -1, 1),
    'HVdPhi_fit_fallback': ('HVdPhi_fit_fallback',            ';#||{#Delta#phi(V, H)};', 16, 0,  3.2),
    'V_pt':                ('V_pt',                           ';p_{T}(V) [GeV];',        40, 0,  400),
    'jjVPtBal_fit':        ('H_pt_fit/V_pt',                  ';fitted jj/V P_{T} Bal;', 30, 0,  2),
    'jjVPtBal':            ('H_pt/V_pt',                      ';jj/V P_{T} Bal;',        30, 0,  2),
    #'BDT':                 ('CMS_vhbb_BDTG_Zll_HighPT_13TeV', ';BDT output;',            20, -1, 1),
    #'DNN':                 ('CMS_vhbb_DNN_Zll_HighPT_13TeV',  ';DNN output;',            20, 0,  1),
}


#######################################
# Samples, Selections, and Categories #
#######################################

# Scale Factors   | HIG-16-044 | 2017 June 12
#--------------------------------------------
# sf_tt           | 1.04       | 0.87774
# sf_zj0b         | 1.02       | 1.2111
# sf_zj1b         | 1.02       | 1.5168
# sf_zj2b         | 1.28       | 1.3206

the_samples_dict = get_samples(
    channel='Zll',
    signal_overlay=False,
    # Note that these scale factors must be
    # appropriate for the high Pt analysis bin.
    #sf_tt  =0.87774,
    #sf_zj0b=1.2111,
    #sf_zj1b=1.5168,
    #sf_zj2b=1.3206,
)

good_ele_1 = '(abs(Electron_eta[lepInd1])<1.44 || abs(Electron_eta[lepInd1])>1.56)'

good_ele_2 = '(abs(Electron_eta[lepInd2])<1.44 || abs(Electron_eta[lepInd2])>1.56)'

CR_ZLF = 'H_mass_fit_fallback>90 && H_mass_fit_fallback<=150 && (controlSample==22 || controlSample==24)'

CR_ZHF = ''.join("""
(H_mass_fit_fallback < 90 || H_mass_fit_fallback >= 150)
&& (controlSample==23
    || controlSample==25
    || (controlSample==0
        && Jet_btagDeepB[hJetInd1]>0.8001 && Jet_btagDeepB[hJetInd2]>0.1522
        && MET_Pt<60
        && HVdPhi_fit_fallback>2.5
        && V_mass>85 && V_mass<=97))
""".split())

SR = 'controlSample==0 && H_pt>100 && H_mass_fit_fallback>90 && H_mass_fit_fallback<150 && Jet_btagDeepB[hJetInd2]>0.1522'

regions = {
    #'SR_ZeeHighPt':     'isZee && sampleIndex!=0 && V_pt>=150 && %s' % SR,
    'CR_ZeeHighPt_TT':  'isZee && V_pt>=150 && controlSample==21 && %s && %s' % (good_ele_1, good_ele_2),
    'CR_ZeeHighPt_ZLF': 'isZee && V_pt>=150 && %s && %s && %s' % (CR_ZLF, good_ele_1, good_ele_2),
    'CR_ZeeHighPt_ZHF': 'isZee && V_pt>=150 && %s && %s && %s' % (CR_ZHF, good_ele_1, good_ele_2),
    #'SR_ZmmHighPt':     'isZmm && sampleIndex!=0 && V_pt>=150 && %s' % SR,
    'CR_ZmmHighPt_TT':  'isZmm && V_pt>=150 && controlSample==21',
    'CR_ZmmHighPt_ZLF': 'isZmm && V_pt>=150 && %s' % CR_ZLF,
    'CR_ZmmHighPt_ZHF': 'isZmm && V_pt>=150 && %s' % CR_ZHF,
}

selections = [
    'Pass_nominal',
    'cutFlow>=2',
    'twoResolvedJets',
    'Jet_bReg[hJetInd2]>30',
    'Jet_bReg[hJetInd1]>30',
]

the_category_dict = {
    'ZllHbbHighPt': [regions, selections, plot_vars],
}

