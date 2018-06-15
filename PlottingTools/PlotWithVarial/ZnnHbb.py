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

name = 'plots_ZnnHbb'

weight = 'weight'

plot_vars = {
    'H_pt':        ('H_pt',                           ';p_{T}(jj) [GeV];',       35, 100, 450),
    'H_mass':      ('H_mass',                         ';M_{jj} [GeV];',          25, 0,   250),
    'MET_Pt':      ('MET_Pt',                         ';MET [GeV];',             30, 150, 450),
    'DeepCSV_max': ('Jet_btagDeepB[hJetInd1]',        ';DeepCSV_{max};',         20, 0,   1),
    'DeepCSV_min': ('Jet_btagDeepB[hJetInd2]',        ';DeepCSV_{min};',         20, 0,   1),
    'CMVA_max':    ('Jet_btagCMVA[hJetInd1]',         ';CMVA_{max};',            20, -1,  1),
    'CMVA_min':    ('Jet_btagCMVA[hJetInd2]',         ';CMVA_{min};',            20, -1,  1),
    'HVdPhi':      ('HVdPhi',                         ';#||{#Delta#phi(V, H)};', 16, 0,   3.2),
    #'BDT':         ('CMS_vhbb_BDTG_Znn_HighPT_13TeV', ';BDT output;',            20, -1,  1),
    #'DNN':         ('CMS_vhbb_DNN_Znn_13TeV',         ';DNN output;',            20, 0,   1),
}


#######################################
# Samples, Selections, and Categories #
#######################################

# Scale Factors   | HIG-16-044 | 2017 June 12
#--------------------------------------------
# sf_tt           | 0.78       | 1.0047
# sf_wj0b         | 1.14       | 1.1311
# sf_wj1b         | 1.66       | 1.9865
# sf_wj2b         | 1.49       | 1.7138
# sf_zj0b         | 1.03       | 0.82682
# sf_zj1b         | 1.28       | 0.55427
# sf_zj2b         | 1.61       | 0.99100

the_samples_dict = get_samples(
    channel='Znn',
    signal_overlay=False,
    #sf_tt  =1.0047,
    #sf_wj0b=1.1311,
    #sf_wj1b=1.9865,
    #sf_wj2b=1.7138,
    #sf_zj0b=0.82682,
    #sf_zj1b=0.55427,
    #sf_zj2b=0.99100,
)

regions = {
    #'SR': 'isZnn && controlSample==0 && (H_mass>60 && H_mass<160) && Jet_btagDeepB[hJetInd2]>0.1522 && sampleIndex!=0',
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

