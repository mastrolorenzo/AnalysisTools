import nano_plotvariables

veto_MH = ' && (H_mass<90 || H_mass>150)'
good_ele_1 = ' && (abs(Electron_eta[lepInd1])<1.44 || abs(Electron_eta[lepInd1])>1.56)'
good_ele_2 = ' && (abs(Electron_eta[lepInd2])<1.44 || abs(Electron_eta[lepInd2])>1.56)'
low_vpt = ' && V_pt >= 50 && V_pt < 150'
high_vpt = ' && V_pt >= 150'

cats_all_regions = {
    'All_Regions':         '1',  # good for debugging ...
}

cats_SR = {
#    'SR_0Lepton':               '(controlSample==0 && isZnn==1   )',
#    'SR_1Lepton_Mu':            '(controlSample==0 && isWmunu==1 )',
#    'SR_1Lepton_El':            '(controlSample==0 && isWenu==1  )',
#    'SR_2Lepton_lowVpt_Mu':     '(controlSample==0 && isZmm==1 %s)' % low_vpt,
    'SR_2Lepton_highVpt_Mu':    '(controlSample==0 && isZmm==1 %s)' % high_vpt,
#    'SR_2Lepton_lowVpt_El':     '(controlSample==0 && isZee==1 %s)' % low_vpt,
#    'SR_2Lepton_highVpt_El':    '(controlSample==0 && isZee==1 %s)' % high_vpt,
}

cats_CR_0Lepton = {
    'CR_0Lepton_TTbar':    '(controlSample==1 && (isWmunu==1 || isWenu == 1))',
#    'CR_0Lepton_VpLF':     '(controlSample==2 && isZnn==1)',
#    'CR_0Lepton_VpHF':     '(controlSample==3 && isZnn==1)',
}

cats_CR_1Lepton = {
    'CR_1Lepton_TTbar_Mu': '(controlSample==11 && isWmunu==1 && Vtype==2)',
 #   'CR_1Lepton_VpLF_Mu':  '(controlSample==12 && isWmunu==1 && Vtype==2)',
 #   'CR_1Lepton_VpHF_Mu':  '(controlSample==13 && isWmunu==1 && Vtype==2%s)' % veto_MH,
    'CR_1Lepton_TTbar_El': '(controlSample==11 && isWenu==1 && Vtype==3%s)' % good_ele_1,
 #   'CR_1Lepton_VpLF_El':  '(controlSample==12 && isWenu==1 && Vtype==3%s)' % good_ele_1,
 #   'CR_1Lepton_VpHF_El':  '(controlSample==13 && isWenu==1 && Vtype==3%s%s)' % (good_ele_1, veto_MH),
}

cats_CR_2Lepton_lowVpt = {
    'CR_2Lepton_lowVpt_TTbar_Mu': '(controlSample==21 && isZmm==1%s)' % low_vpt,
    'CR_2Lepton_lowVpt_VpLF_Mu':  '(controlSample==22 && isZmm==1%s)' % low_vpt,
    'CR_2Lepton_lowVpt_VpHF_Mu':  '(controlSample==23 && isZmm==1%s)' % low_vpt,
    'CR_2Lepton_lowVpt_TTbar_El': '(controlSample==21 && isZee==1%s%s%s)' % (good_ele_1, good_ele_2, low_vpt),
    'CR_2Lepton_lowVpt_VpLF_El':  '(controlSample==22 && isZee==1%s%s%s)' % (good_ele_1, good_ele_2, low_vpt),
    'CR_2Lepton_lowVpt_VpHF_El':  '(controlSample==23 && isZee==1%s%s%s)' % (good_ele_1, good_ele_2, low_vpt),
}
cats_CR_2Lepton_highVpt = {
    'CR_2Lepton_highVpt_TTbar_Mu': '(controlSample==21 && isZmm==1%s)' % high_vpt,
    'CR_2Lepton_highVpt_VpLF_Mu':  '(controlSample==22 && isZmm==1%s)' % high_vpt,
    'CR_2Lepton_highVpt_VpHF_Mu':  '(controlSample==23 && isZmm==1%s)' % high_vpt,
    'CR_2Lepton_highVpt_TTbar_El': '(controlSample==21 && isZee==1%s%s%s)' % (good_ele_1, good_ele_2, high_vpt),
    'CR_2Lepton_highVpt_VpLF_El':  '(controlSample==22 && isZee==1%s%s%s)' % (good_ele_1, good_ele_2, high_vpt),
    'CR_2Lepton_highVpt_VpHF_El':  '(controlSample==23 && isZee==1%s%s%s)' % (good_ele_1, good_ele_2, high_vpt),
}
cats_CR_TT_boosted = {
        'CR_TT_boosted': '(boostedControlSample==11)'
}
cats_CR_Vlight_boosted = {
        'CR_Vlight_boosted': '(boostedControlSample==12)'
}
cats_CR_Vhf_boosted = {
        'CR_Vhf_boosted': '(boostedControlSample==13)'
}
cats_CR_Vhf_boosted_El = {
        'CR_Vhf_boosted_El': '(boostedControlSample==13)&&(isWenu)'
}
cats_CR_Vhf_boosted_Mu = {
        'CR_Vhf_boosted_Mu': '(boostedControlSample==13)&&(isWmunu)'
}

cats_SR_1Lepton_boosted = {
        'SR_1Lepton_boosted': '(boostedControlSample==0)&&(sampleIndex!=0)'#&&(V_pt>275)&&(FatJetCand_tau32>.4)'#blind the SR
}

cats_SR_1Lepton_boosted_high_btag = {
        'SR_1Lepton_boosted_high_btag': '(boostedControlSample==0)&&(sampleIndex!=0)&&(FatJetCand_doubleB>0.97)'#blind the SR
}

cats_SR_1Lepton_boosted_low_btag = {
        'SR_1Lepton_boosted_low_btag': '(boostedControlSample==0)&&(sampleIndex!=0)&&(FatJetCand_doubleB<=0.97)'#&&(V_pt>300)&&(FatJetCand_tau32>.4)',#blind the SR
}
boosted_sel = ['boostedCategory>0'] #boostedCategory>0 means plot overlap and boosted events together
boosted_SR_sel = boosted_sel+[]

boosted_SR_low_btag_sel = boosted_SR_sel

boosted_SR_high_btag_sel = boosted_SR_sel

no_sel = [
    'Pass_nominal',
    'cutFlow>=2',
    'Jet_bReg[hJetInd2] > 30',
    'Jet_bReg[hJetInd1] > 30',
    #'Jet_Pt[hJetInd1] > 20',
    #'Jet_Pt[hJetInd2] > 20',
]
sr_sel = no_sel + [  # copied from selection.dat
    'H_pt>100',
    'H_mass>90',
    'H_mass<150',
]

# CS 1lep
cs_1lep_sel = no_sel + [
    'H_pt>100',
    'V_pt>100',
    #'selLeptons_relIso_0<0.06',
   # 'MET_Pt<170',
]

#cs_0lep_sel = no_sel + [
#    #'MET_Pt>170',
#    #'H_pt>120',
#]

# CS 2lep
# cs_2lep_sel = no_sel + [
#     'H_pt>100',
#     '(%s)<-0.5884'%nano_plotvariables.CMVAmax,
#     'HVdPhi>2.5',
#     'H_mass>90',
#     'H_mass<150',
# ]

# this dictionary is used for setting up the plotting tools
the_category_dict = {
#   region block     categories               sel.         histograms
  #  'AllRegions':   [cats_all_regions,        no_sel,      nano_plotvariables.vars_used_in_selections],
#    'CR_2LeptLO':   [cats_CR_2Lepton_lowVpt,  no_sel,      nano_plotvariables.test],
#    'CR_2LeptLO':   [cats_CR_2Lepton_lowVpt,  no_sel,      nano_plotvariables.vars_2lepCR],
#  'CR_2LeptHI':   [cats_CR_2Lepton_highVpt, no_sel,      nano_plotvariables.vars_2lepCR],
   # 'CR_2LeptHI':   [cats_CR_2Lepton_highVpt, no_sel,      nano_plotvariables.test],
#  'CR_1Lepton':   [cats_CR_1Lepton,         cs_1lep_sel, nano_plotvariables.vars_1lepCR],
#  'CR_1Lepton':   [cats_CR_1Lepton,         cs_1lep_sel, nano_plotvariables.test],
#   'CR_0Lepton':   [cats_CR_0Lepton,         no_sel,      nano_plotvariables.test],
#   'CR_0Lepton':   [cats_CR_0Lepton,         no_sel,      nano_plotvariables.vars_used_in_selections],
 #   'SR':           [cats_SR,                 sr_sel,      nano_plotvariables.vars_used_in_selections],
   # 'SR':           [cats_SR,                 sr_sel,      nano_plotvariables.test],
  'CR_TT_boosted':   [cats_CR_TT_boosted,         boosted_sel, nano_plotvariables.test],
  'CR_Vlight_boosted':   [cats_CR_Vlight_boosted,         boosted_sel, nano_plotvariables.test],
  'CR_Vhf_boosted':   [cats_CR_Vhf_boosted,         boosted_sel, nano_plotvariables.test],
  'CR_Vhf_boosted_El':   [cats_CR_Vhf_boosted_El,         boosted_sel, nano_plotvariables.test],
  'CR_Vhf_boosted_Mu':   [cats_CR_Vhf_boosted_Mu,         boosted_sel, nano_plotvariables.test],
  'SR_1Lepton_boosted':   [cats_SR_1Lepton_boosted,         boosted_SR_sel, nano_plotvariables.test],
  'SR_1Lepton_boosted_high_btag':   [cats_SR_1Lepton_boosted_high_btag,         boosted_SR_high_btag_sel, nano_plotvariables.test],
  'SR_1Lepton_boosted_low_btag':   [cats_SR_1Lepton_boosted_low_btag,         boosted_SR_low_btag_sel, nano_plotvariables.test],
}
