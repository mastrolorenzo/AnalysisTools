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

name = 'plots_WlnHbb'

weight = 'weight'

plot_vars = {
    'H_pt':        ('H_pt',                    ';p_{T}(jj) [GeV];',       30, 50, 350),
    'H_mass':      ('H_mass',                  ';M_{jj} [GeV];',          25, 0,  250),
    'MET_Pt':      ('MET_Pt',                  ';MET [GeV];',             20, 0,  200),
    'DeepCSV_max': ('Jet_btagDeepB[hJetInd1]', ';DeepCSV_{max};',         20, 0,  1),
    'DeepCSV_min': ('Jet_btagDeepB[hJetInd2]', ';DeepCSV_{min};',         20, 0,  1),
    'HVdPhi':      ('HVdPhi',                  ';#||{#Delta#phi(V, H)};', 16, 0,  3.2),
    'rho':         ('fixedGridRhoFastjetAll',  ';#rho;',                  60, 0,  60),
    'V_pt':        ('V_pt',                    ';p_{T}(V) [GeV];',        30, 50, 350),
}


#######################################
# Samples, Selections, and Categories #
#######################################

the_samples_dict = get_samples(
    channel='Wln',
    signal_overlay=False,
    #sf_tt=0.91,
    #sf_wj0b=1.14,
    #sf_wj1b=1.66,
    #sf_wj2b=1.49,
)

regions = {
    #'SR_Wen':     'isWenu && controlSample==0 && (H_mass>90&&H_mass<150) && Jet_btagDeepB[hJetInd2]>0.1522',
    'CR_Wen_TT':  'isWenu && controlSample==11',
    'CR_Wen_WLF': 'isWenu && controlSample==12',
    'CR_Wen_WHF': 'isWenu && controlSample==13',
    #'SR_Wmn':     'isWmunu && controlSample==0 && (H_mass>90&&H_mass<150) && Jet_btagDeepB[hJetInd2]>0.1522',
    'CR_Wmn_TT':  'isWmunu && controlSample==11',
    'CR_Wmn_WLF': 'isWmunu && controlSample==12',
    'CR_Wmn_WHF': 'isWmunu && controlSample==13',
}

selections = [
    'Pass_nominal',
    'cutFlow>=2',
    'twoResolvedJets',
]

the_category_dict = {
    'WlnHbb': [regions, selections, plot_vars],
}

