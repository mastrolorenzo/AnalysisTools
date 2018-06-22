from samples_2016 import *
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

name = 'plots_ZnnHbb'

weight = 'weight'

plot_vars = {
    'H_pt':                          ('H_pt',                           ';p_{T}(jj) [GeV];',             35, 100, 450),
    'H_mass':                        ('H_mass',                         ';M(jj) [GeV];' ,                25,   0, 250),
    'MET_Pt':                        ('MET_Pt',                         ';MET [GeV];',                   30, 150, 450),
    'DeepCSV_max':                   ('Jet_btagDeepB[hJetInd1]',        ';DeepCSV_{max};',               20,   0,   1),
    'DeepCSV_max_zoom':              ('Jet_btagDeepB[hJetInd1]',        ';DeepCSV_{max};',               20, 0.8,   1),
    'DeepCSV_min':                   ('Jet_btagDeepB[hJetInd2]',        ';DeepCSV_{min};',               20,   0,   1),
    'CMVA_max':                      ('Jet_btagCMVA[hJetInd1]',         ';CMVA_{max};',                  20,  -1,   1),
    'CMVA_max_zoom':                 ('Jet_btagCMVA[hJetInd1]',         ';CMVA_{max};',                  20, 0.9,   1),
    'CMVA_min':                      ('Jet_btagCMVA[hJetInd2]',         ';CMVA_{min};',                  20,  -1,   1),
    'HVdPhi':                        ('HVdPhi',                         ';{#Delta#phi(V, H);',           16,   0, 3.2),
    #'BDT':                           ('CMS_vhbb_BDTG_Znn_HighPT_13TeV', ';BDT output;',                  20,  -1,   1),
    #'DNN':                           ('CMS_vhbb_DNN_Znn_13TeV',         ';DNN output;',                  20,   0,   1),
    'mvaInput_H_mass':               ('H_mass',                         ';M(jj) [GeV];',                 25,   0, 250),
    'mvaInput_H_pt':                 ('H_pt',                           ';p_{T}(jj) [GeV];',             35, 100, 450),
    'mvaInput_HVdPhi':               ('HVdPhi',                         ';#Delta#phi(V, H);',            16,   0, 3.2),
    'mvaInput_V_pt':                 ('V_pt',                           ';p_{T}(V) [GeV];',              30, 150, 450),
    'mvaInput_HJ1_HJ2_dEta':         ('HJ1_HJ2_dEta',                   ';#Delta#eta(jj);',              40,   0,   4),
    'mvaInput_hJets_CMVA_0':         ('hJets_btagged_0',                ';CMVAv2_{max};',                20,  -1,   1),
    'mvaInput_hJets_CMVA_0_zoom':    ('hJets_btagged_0',                ';CMVAv2_{max};',                20, 0.9,   1),
    'mvaInput_hJets_CMVA_1':         ('hJets_btagged_1',                ';CMVAv2_{min};',                20,  -1,   1),
    'mvaInput_hJets_DeepCSV_0':      ('hJets_btagged_0',                ';DeepCSV_{max};',               20,   0,   1),
    'mvaInput_hJets_DeepCSV_0_zoom': ('hJets_btagged_0',                ';DeepCSV_{max};',               20, 0.8,   1),
    'mvaInput_hJets_DeepCSV_1':      ('hJets_btagged_1',                ';DeepCSV_{min};',               20,   0,   1),
    'mvaInput_SA5':                  ('SA5',                            ';SA5;',                          7,   0,   7),
    'mvaInput_HJ1_HJ2_dPhi':         ('HJ1_HJ2_dPhi',                   ';#Delta#phi(jj);',              16,   0, 3.2),
    'mvaInput_hJets_leadingPt':      ('hJets_leadingPt',                ';p_{T}(j_{1}) [GeV];',          35,  50, 350),
    'mvaInput_hJets_subleadingPt':   ('hJets_subleadingPt',             ';p_{T}(j_{2}) [GeV];',          25,   0, 250),
    'mvaInput_otherJetsBestCMVA':    ('otherJetsBestBtag',              ';CMVAv2_{add};',                20,  -1,   1),
    'mvaInput_otherJetsBestDeepCSV': ('otherJetsBestBtag',              ';DeepCSV_{add};',               20,  -1,   1),
    'mvaInput_otherJetsHighestPt':   ('otherJetsHighestPt',             ';p_{T}(add) [GeV];',            25,   0, 250),
    'mvaInput_minDPhiFromOtherJets': ('minDPhiFromOtherJets',           ';#Delta#phi(E_{T}^{miss}, j);', 16,   0, 3.2),
}


#######################################
# Samples, Selections, and Categories #
#######################################

# Scale Factors | HIG-16-044 | 2016, June 19 | 2017, June 19
#--------------------------------------------|---------------
# sf_tt         | 0.78       | 0.83          | 1.01
# sf_wj0b       | 1.14       | 1.06          | 1.13
# sf_wj1b       | 1.66       | 1.13          | 2.00
# sf_wj2b       | 1.49       | 1.60          | 1.75
# sf_zj0b       | 1.03       | 1.44          | 0.81
# sf_zj1b       | 1.28       | 2.24          | 0.97
# sf_zj2b       | 1.61       | 1.94          | 1.17

the_samples_dict = get_samples(
    channel='Znn',
    signal_overlay=False,
    sf_tt  =0.83,
    sf_wj0b=1.06,
    sf_wj1b=1.13,
    sf_wj2b=1.60,
    sf_zj0b=1.44,
    sf_zj1b=2.24,
    sf_zj2b=1.94,
)

regions = {
    # 2016
    'SR': 'isZnn && controlSample==0 && (H_mass>60 && H_mass<160) && hJets_btagged_0>0.9430 && hJets_btagged_1>-0.5884',
    # 2017
    #'SR': 'isZnn && controlSample==0 && (H_mass>60 && H_mass<160) && hJets_btagged_0>0.8001 && hJets_btagged_1>0.1522',
    'CR_Znn_TT':  '(isWenu || is Wmunu) && controlSample==1',
    'CR_Znn_ZLF': 'isZnn && controlSample==2',
    'CR_Znn_ZHF': 'isZnn && controlSample==3',
}

selections = [
    'Pass_nominal',
    'cutFlow>=2',
    'twoResolvedJets',
]

the_category_dict = {
    'ZnnHbb': [regions, selections, plot_vars],
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

