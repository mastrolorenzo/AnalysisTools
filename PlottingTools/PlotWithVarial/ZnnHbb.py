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
    'H_pt':        ('H_pt',                    ';p_{T}(jj) [GeV];',       35, 100, 450),
    'H_mass':      ('H_mass',                  ';M_{jj} [GeV];',          25, 0,   250),
    'MET_Pt':      ('MET_Pt',                  ';MET [GeV];',             30, 150, 450),
    'DeepCSV_max': ('Jet_btagDeepB[hJetInd1]', ';DeepCSV_{max};',         20, 0,     1),
    'DeepCSV_min': ('Jet_btagDeepB[hJetInd2]', ';DeepCSV_{min};',         20, 0,     1),
    'HVdPhi':      ('HVdPhi',                  ';#||{#Delta#phi(V, H)};', 16, 0,   3.2),
}


#######################################
# Samples, Selections, and Categories #
#######################################

the_samples_dict = get_samples(
    channel='Znn',
    signal_overlay=False,
    #sf_tt=0.78,
    #sf_wj0b=1.14,
    #sf_wj1b=1.66,
    #sf_wj2b=1.49,
    #sf_zj0b=1.03,
    #sf_zj1b=1.28,
    #sf_zj2b=1.61,
)

regions = {
    #'SR': 'isZnn && controlSample==0 && (H_mass>60 && H_mass<160) && Jet_btagDeepB[hJetInd2]>0.1522',
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

