from samples_2018 import *
#from samples_2017_2 import *
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

#name = 'plot_ZllHbbHighPt_2017data'
name = 'plot_ZllHbbHighPt'

weight = 'weight'

plot_vars = {
    'H_pt_fit_fallback':           ('H_pt_fit_fallback',             ';p_{T}(jj) [GeV];',       35,   0, 350),
    'H_mass':                      ('H_mass',                        ';M(jj) [GeV];',           25,   0, 250),
    'H_pt':                      ('H_pt',                        ';p_{T}(jj) [GeV];',           35,   0, 350),
    'H_mass_fit_fallback':         ('H_mass_fit_fallback',           ';fitted M(jj) [GeV];',    25,   0, 250),
    'MET_Pt':                      ('MET_Pt',                        ';MET [GeV];',             25,   0, 250),
    'JET_Pt':                      ('Jet_pt',                        ';p_T(Jet) [GeV];',        35,   0, 900),
    'leading_JET_Pt':              ('Jet_pt[hJetInd1]',              ';p_T(Jet) [GeV];',        35,   0, 900),
    'subleading_JET_Pt':           ('Jet_pt[hJetInd2]',              ';p_T(Jet) [GeV];',        35,   0, 900),
    'leading_Muon_Pt':              ('Muon_pt[lepInd1]',              ';p_T(l) [GeV];',        35,   0, 600),
    'subleading_Muon_Pt':           ('Muon_pt[lepInd2]',              ';p_T(l) [GeV];',        35,   0, 600),
    'leading_JET_Pt_bReg':              ('Jet_bReg[hJetInd1]',              ';p_T(Jet) [GeV];',        35,   0, 900),
    'subleading_JET_Pt_bReg':           ('Jet_bReg[hJetInd2]',              ';p_T(Jet) [GeV];',        35,   0, 900),
    'leading_JET_eta':              ('Jet_eta[hJetInd1]',             ';#eta(Jet) [GeV];',       60,   -3, 3),
    'subleading_JET_eta':           ('Jet_eta[hJetInd2]',              ';#eta(Jet) [GeV];',       60,   -3, 3),
    'JET_Eta':                     ('Jet_eta',                       ';#eta(Jet) [GeV];',       35,   -3, 3),
    'JET_Phi':                     ('Jet_phi',                       ';#phi(Jet) [GeV];',        35,   0, 5),
    'JET_mass':                     ('Jet_mass',                     ';M(Jet) [GeV];',        35,   0, 300),
    'JET_Pt_bReg':                 ('Jet_bReg',                      ';p_T(Jet) [GeV];',        35,   0, 900),
    'JET_Pt_bReg2':                 ('Jet_bReg',                      ';p_T(Jet) [GeV];',        35,   20, 900),
    'DeepCSV_max':                 ('Jet_btagDeepB[hJetInd1]',       ';DeepCSV_{max};',         20,   0,   1),
    'DeepCSV_max_zoom':            ('Jet_btagDeepB[hJetInd1]',       ';DeepCSV_{max};',         20, 0.8,   1),
    'DeepCSV_min':                 ('Jet_btagDeepB[hJetInd2]',       ';DeepCSV_{min};',         20,   0,   1),
    'CMVA_max':                    ('Jet_btagCMVA[hJetInd1]',        ';CMVA_{max};',            20,  -1,   1),
    'CMVA_max_zoom':               ('Jet_btagCMVA[hJetInd1]',        ';CMVA_{max};',            20, 0.9,   1),
    'CMVA_min':                    ('Jet_btagCMVA[hJetInd2]',        ';CMVA_{min};',            20,  -1,   1),
    'HVdPhi_fit_fallback':         ('HVdPhi_fit_fallback',           ';#Delta#phi(V, H);',      16,   0, 3.2),
    'V_pt':                        ('V_pt',                          ';p_{T}(V) [GeV];',        40,   0, 400),
    'jjVPtBal_fit':                ('H_pt_fit_fallback/V_pt',        ';fitted jj/V P_{T} Bal;', 30,   0,   2),
    'jjVPtBal':                    ('H_pt/V_pt',                     ';jj/V P_{T} Bal;',        30,   0,   2),
    #'BDT':                         ('CMS_vhbb_BDTG_Zll_LowPT_13TeV', ';BDT output;',            20,  -1,   1),
    #'DNN':                         ('CMS_vhbb_DNN_Zll_LowPT_13TeV',  ';DNN output;',            20,   0,   1),
    'mvaInput_H_mass':             ('H_mass_fit_fallback',           ';M(jj) [GeV];',           25,   0, 250),
    'mvaInput_H_pt':               ('H_pt_fit_fallback',             ';p_{T}(jj) [GeV];',       35, 100, 450),
    'mvaInput_V_pt':               ('V_pt',                          ';p_{T}(V) [GeV];',        40,   0, 400),
    'mvaInput_hJets_CMVA_0':       ('hJets_btagged_0',               ';CMVAv2_{max};',          20,  -1,   1),
    'mvaInput_hJets_CMVA_1':       ('hJets_btagged_1',               ';CMVAv2_{min};',          20,  -1,   1),
    'mvaInput_hJets_DeepCSV_0':    ('hJets_btagged_0',               ';DeepCSV_{max};',         20,   0,   1),
    'mvaInput_hJets_DeepCSV_1':    ('hJets_btagged_1',               ';DeepCSV_{min};',         20,   0,   1),
    'mvaInput_nAddJets_2lep':      ('nAddJets_2lep',                 ';N_{aj};',                10,   0,  10),
    'mvaInput_SA5':                ('SA5',                           ';SA5;',                    7,   0,   7),
    'mvaInput_V_mass':             ('V_mass',                        ';M(V) [GeV];',            20,  80, 100),
    'mvaInput_MET_Pt':             ('MET_Pt',                        ';E_{T}^{miss} [GeV];',    25,   0, 250),
    'mvaInput_hJets_leadingPt':    ('hJets_leadingPt',               ';p_{T}(j_{1}) [GeV];',    25,   0, 250),
    'mvaInput_hJets_subleadingPt': ('hJets_subleadingPt',            ';p_{T}(j_{2}) [GeV];',    25,   0, 250),
    'mvaInput_jjVPtRatio':         ('jjVPtRatio_fit_fallback',       ';p_{T}(jj) / p_{T}(V);',  20,   0,   2),
    'mvaInput_HJ1_HJ2_dEta':       ('HJ1_HJ2_dEta',                  ';#Delta#eta(jj);',        40,   0,   4),
    'mvaInput_HVdPhi':             ('HVdPhi_fit_fallback',           ';#Delta#phi(V, H);',      16,   0, 3.2),
    'mvaInput_n_recoil_jets_fit':  ('n_recoil_jets_fit',             ';N_{recoil jets};',       10,   0,  10),
    'mvaInput_H_mass_sigma_fit':   ('H_mass_sigma_fit',              ';#sigma(M(jj)) [GeV];',   25,   0,  50),
    'rho':                         ('fixedGridRhoFastjetAll',        ';#rho;',                  60,   0,  60),
    'CMS_vhbb_DNN_Wen_13TeV':                    ('CMS_vhbb_DNN_Wen_13TeV',        ';DNN;',            20,  -1,   1),
    'CMS_vhbb_DNN_Wmn_13TeV':                    ('CMS_vhbb_DNN_Wmn_13TeV',        ';DNN;',            20,  -1,   1),
    'CMS_vhbb_DNN_Zll_HighPT_13TeV':                    ('CMS_vhbb_DNN_Zll_HighPT_13TeV',        ';DNN;',            20,  -1,   1),
    'CMS_vhbb_DNN_Zll_lowPT_13TeV':                    ('CMS_vhbb_DNN_Zll_lowPT_13TeV',        ';DNN;',            20,  -1,   1),
    'CMS_vhbb_BDTG_Zll_HighPT_13TeV':                    ('CMS_vhbb_BDTG_Zll_HighPT_13TeV',        ';DNN;',            20,  -1,   1),
    'CMS_vhbb_BDTG_Zll_LowPT_13TeV':                    ('CMS_vhbb_BDTG_Zll_LowPT_13TeV',        ';DNN;',            20,  -1,   1),
    'CMS_vhbb_DNN_Znn_13TeV':                    ('CMS_vhbb_DNN_Znn_13TeV',        ';DNN;',            20,  -1,   1),

}


#######################################
# Samples, Selections, and Categories #
#######################################

# Scale Factors | HIG-16-044 | 2016, June 19 | 2017, June 19
#--------------------------------------------|---------------
# sf_tt         | 1.04       | 0.73          | 0.94
# sf_zj0b       | 1.02       | 0.98          | 1.21
# sf_zj1b       | 1.02       | 1.31          | 1.53
# sf_zj2b       | 1.28       | 1.39          | 1.41

the_samples_dict = get_samples(
    channel='Zll',
    signal_overlay=False,
    # Note that these scale factors must be
    # appropriate for the high Pt analysis bin.
    sf_tt  =1,
    sf_zj0b=1,
    sf_zj1b=1,
    sf_zj2b=1,
)

good_ele_1 = '(abs(Electron_eta[lepInd1])<1.44 || abs(Electron_eta[lepInd1])>1.56)'

good_ele_2 = '(abs(Electron_eta[lepInd2])<1.44 || abs(Electron_eta[lepInd2])>1.56)'

CR_ZLF = 'H_mass_fit_fallback>90 && H_mass_fit_fallback<=150 && (controlSample==22 || controlSample==24)'
#CR_ZLF = ' (controlSample==22 || controlSample==24)'

CR_ZHF_2016 = ''.join("""
(H_mass_fit_fallback < 90 || H_mass_fit_fallback >= 150)
&& (controlSample==23
    || controlSample==25
    || (controlSample==0
        && hJets_btagged_0>0.9430 && hJets_btagged_1>-0.5884
        && MET_Pt<60
        && HVdPhi_fit_fallback>2.5
        && V_mass>85 && V_mass<=97))
""".split())

CR_ZHF_2017 = ''.join("""
(H_mass_fit_fallback < 90 || H_mass_fit_fallback >= 150)
&& (controlSample==23
    || controlSample==25
    || (controlSample==0
        && hJets_btagged_0>0.8001 && hJets_btagged_1>0.1522
        && MET_Pt<60
        && HVdPhi_fit_fallback>2.5
        && V_mass>85 && V_mass<=97))
""".split())

#CR_ZHF_2017 = ''.join("""
#(controlSample==23
#    || controlSample==25
#    || (controlSample==0
#        && hJets_btagged_0>0.8001 && hJets_btagged_1>0.1522
#        && MET_Pt<60
#        && HVdPhi_fit_fallback>2.5
#        && V_mass>85 && V_mass<=97))
#""".split())


SR_2016 = 'controlSample==0 && H_mass_fit_fallback>90 && H_mass_fit_fallback<150 && hJets_btagged_0>-0.5884 && hJets_btagged_1>-0.5884'

SR_2017 = 'controlSample==0 && H_mass_fit_fallback>90 && H_mass_fit_fallback<150 && hJets_btagged_0>0.1522 && hJets_btagged_1>0.1522'
#SR_2017 = 'controlSample==0 && hJets_btagged_0>0.1522 && hJets_btagged_1>0.1522'

regions = {
    # 2016
    #'SR_ZeeHighPt':     'isZee && V_pt>=150 && %s' % SR_2016,
    #'SR_ZmmHighPt':     'isZmm && V_pt>=150 && %s' % SR_2016,
    #'CR_ZeeHighPt_ZHF': 'isZee && V_pt>=150 && %s && %s && %s' % (CR_ZHF_2016, good_ele_1, good_ele_2),
    #'CR_ZmmHighPt_ZHF': 'isZmm && V_pt>=150 && %s' % CR_ZHF_2016,
    # 2017
#    'SR_ZeeHighPt':     'isZee && V_pt>=150 && %s' % SR_2017,
#    'SR_ZmmHighPt':     'isZmm && V_pt>=150 && %s' % SR_2017,
    'CR_ZeeHighPt_ZHF': 'isZee && V_pt>=150 && %s && %s && %s' % (CR_ZHF_2017, good_ele_1, good_ele_2),
    'CR_ZmmHighPt_ZHF': 'isZmm && V_pt>=150 && %s' % CR_ZHF_2017,
    'CR_ZeeHighPt_TT':  'isZee && V_pt>=150 && controlSample==21 && %s && %s' % (good_ele_1, good_ele_2),
    'CR_ZeeHighPt_ZLF': 'isZee && V_pt>=150 && %s && %s && %s' % (CR_ZLF, good_ele_1, good_ele_2),
    'CR_ZmmHighPt_TT':  'isZmm && V_pt>=150 && controlSample==21',
    'CR_ZmmHighPt_ZLF': 'isZmm && V_pt>=150 && %s' % CR_ZLF,
}

selections = [
    'Pass_nominal',
    'cutFlow>=2',
    'twoResolvedJets',
    #'Jet_bReg[hJetInd1]>30',
    #'Jet_bReg[hJetInd2]>30',
]

the_category_dict = {
    'ZllHbbHighPt': [regions, selections, plot_vars],
}

# Blinding the Signal Region
def additional_input_hook(wrps):

    @varial.history.track_history
    def blind_H_mass_in_SR(w):
        if w.legend == 'Data' and w.in_file_path.startswith('SR'):
            if 'H_mass' in w.name:
                print 'BLINDING Data in %s' % w.in_file_path
                for i in xrange(w.histo.GetNbinsX() + 1):
                    w.histo.SetBinContent(i, 0.)
                    w.histo.SetBinError(i, 0.)
        return w

    wrps = (blind_H_mass_in_SR(w) for w in wrps)
    return wrps

