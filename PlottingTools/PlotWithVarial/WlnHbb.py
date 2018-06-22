from samples_2017 import *
import varial


####################
# General Settings #
####################

varial.settings.rootfile_postfixes += ['.pdf', '.png']
varial.settings.stacking_order = stacking_order
# try to find output on disk and don't run a step if present
enable_reuse_step = True


#################
# Plot Settings #
#################

name = 'plots_WlnHbb'

weight = 'weight'

plot_vars = {
    'H_pt':                          ('H_pt',                             ';p_{T}(jj) [GeV];',             35,  50, 400),
    'H_mass':                        ('H_mass',                           ';M(jj) [GeV];',                 25,   0, 250),
    'MET_Pt':                        ('MET_Pt',                           ';MET [GeV];',                   20,   0, 200),
    'DeepCSV_max':                   ('Jet_btagDeepB[hJetInd1]',          ';DeepCSV_{max};',               20,   0,   1),
    'DeepCSV_max_zoom':              ('Jet_btagDeepB[hJetInd1]',          ';DeepCSV_{max};',               20, 0.8,   1),
    'DeepCSV_min':                   ('Jet_btagDeepB[hJetInd2]',          ';DeepCSV_{min};',               20,   0,   1),
    'CMVA_max':                      ('Jet_btagCMVA[hJetInd1]',           ';CMVA_{max};',                  20,  -1,   1),
    'CMVA_min':                      ('Jet_btagCMVA[hJetInd2]',           ';CMVA_{min};',                  20,  -1,   1),
    'HVdPhi':                        ('HVdPhi',                           ';#Delta#phi(V, H);',            16,   0, 3.2),
    'rho':                           ('fixedGridRhoFastjetAll',           ';#rho;',                        60,   0,  60),
    'V_pt':                          ('V_pt',                             ';p_{T}(V) [GeV];',              25, 100, 350),
    #'BDT':                           ('CMS_vhbb_BDTG_Wln_13TeV',          ';BDT output;',                  20,  -1,   1),
    #'DNN_Wen':                       ('CMS_vhbb_DNN_Wen_13TeV',           ';DNN (Wen) output;',            20,   0,   1),
    #'DNN_Wmn':                       ('CMS_vhbb_DNN_Wmn_13TeV',           ';DNN (Wmn) output;',            20,   0,   1),
    'mvaInput_H_mass':               ('H_mass',                           ';M(jj) [GeV];',                 25,   0, 250),
    'mvaInput_H_pt':                 ('H_pt',                             ';p_{T}(jj) [GeV];',             35,  50, 400),
    'mvaInput_jjVPtRatio':           ('jjVPtRatio',                       ';p_{T}(jj) / p_{T}(V);',        20,   0,   2),
    'mvaInput_V_pt':                 ('V_pt',                             ';p_{T}(V) [GeV];',              25, 100, 350),
    'mvaInput_hJets_CMVA_0':         ('hJets_btagged_0',                  ';CMVAv2_{max};',                20,  -1,   1),
    'mvaInput_hJets_CMVA_0_zoom':    ('hJets_btagged_0',                  ';CMVAv2_{max};',                20, 0.9,   1),
    'mvaInput_hJets_CMVA_1':         ('hJets_btagged_1',                  ';CMVAv2_{min};',                20,  -1,   1),
    'mvaInput_hJets_DeepCSV_0':      ('hJets_btagged_0',                  ';DeepCSV_{max};',               20,   0,   1),
    'mvaInput_hJets_DeepCSV_0_zoom': ('hJets_btagged_0',                  ';DeepCSV_{max};',               20, 0.8,   1),
    'mvaInput_hJets_DeepCSV_1':      ('hJets_btagged_1',                  ';DeepCSV_{min};',               20,   0,   1),
    'mvaInput_Top_mass':             ('Top1_mass_fromLepton_regPT_w4MET', ';M_{t} [GeV];',                 30,  50, 350),
    'mvaInput_HVdPhi':               ('HVdPhi',                           ';#Delta#phi(V, H);',            16,   0, 3.2),
    'mvaInput_nAddJets302p5_puid':   ('nAddJets302p5_puid',               ';N_{aj};',                      10,   0,  10),
    'mvaInput_SA5':                  ('SA5',                              ';SA5;',                          7,   0,   7),
    'mvaInput_lepMetDPhi':           ('lepMetDPhi',                       ';#Delta#phi(E_{T}^{miss}, l);', 16,   0, 3.2),
    'mvaInput_V_mt':                 ('V_mt',                             ';m_{T} [GeV];',                 25,   0, 250),
    'mvaInput_MET_Pt':               ('MET_Pt',                           ';E_{T}^{miss} [GeV];',          20,   0, 200),
    'mvaInput_hJets_leadingPt':      ('hJets_leadingPt',                  ';p_{T}(j_{1}) [GeV];',          35,   0, 350),
    'mvaInput_hJets_subleadingPt':   ('hJets_subleadingPt',               ';p_{T}(j_{2}) [GeV];',          20,   0, 200),
    'mvaInput_HJ1_HJ2_dEta':         ('HJ1_HJ2_dEta',                     ';#Delta#eta(jj);',              40,   0,   4),
    'mvaInput_HJ1_HJ2_dR':           ('HJ1_HJ2_dR',                       ';#DeltaR(jj);',                 20,   0,   5),
}


#######################################
# Samples, Selections, and Categories #
#######################################

# Scale Factors | HIG-16-044 | 2016, June 19 | 2017, June 19
#--------------------------------------------|---------------
# sf_tt         | 0.91       | 0.84          | 1.02
# sf_wj0b       | 1.14       | 1.06          | 1.13
# sf_wj1b       | 1.66       | 1.13          | 2.00
# sf_wj2b       | 1.49       | 1.60          | 1.75

the_samples_dict = get_samples(
    channel='Wln',
    signal_overlay=False,
    sf_tt  =1.02,
    sf_wj0b=1.13,
    sf_wj1b=2.00,
    sf_wj2b=1.75,
)

regions = {
    # 2016
    #'SR_Wen':     'isWenu && controlSample==0 && (H_mass>90&&H_mass<150) && hJets_btagged_0>0.9430 && hJets_btagged_1>-0.5884',
    #'SR_Wmn':     'isWmunu && controlSample==0 && (H_mass>90&&H_mass<150) && hJets_btagged_0>0.9430 && hJets_btagged_1>-0.5884',
    # 2017
    'SR_Wen':     'isWenu && controlSample==0 && (H_mass>90&&H_mass<150) && hJets_btagged_0>0.8001 && hJets_btagged_1>0.1522',
    'SR_Wmn':     'isWmunu && controlSample==0 && (H_mass>90&&H_mass<150) && hJets_btagged_0>0.8001 && hJets_btagged_1>0.1522',
    'CR_Wen_TT':  'isWenu && controlSample==11',
    'CR_Wen_WLF': 'isWenu && controlSample==12',
    'CR_Wen_WHF': 'isWenu && controlSample==13',
    'CR_Wmn_TT':  'isWmunu && controlSample==11',
    'CR_Wmn_WLF': 'isWmunu && controlSample==12',
    'CR_Wmn_WHF': 'isWmunu && controlSample==13',
}

selections = [
    'Pass_nominal',
    'cutFlow>=2',
    'twoResolvedJets',
    'V_pt>150',
]

the_category_dict = {
    'WlnHbb': [regions, selections, plot_vars],
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

