import nano_2lep_plotvariables

veto_MH = ' && (H_mass<90 || H_mass>150)'
good_ele_1 = ' && (abs(Electron_eta[lepInd1])<1.44 || abs(Electron_eta[lepInd1])>1.56)'
good_ele_2 = ' && (abs(Electron_eta[lepInd2])<1.44 || abs(Electron_eta[lepInd2])>1.56)'
low_vpt = ' && V_pt >= 50 && V_pt < 150'
high_vpt = ' && V_pt >= 150'

cats_debug = {
    'DEBUG': '1',
}

no_sel = [
    'Pass_nominal',
    'twoResolvedJets',
    'cutFlow>=2',
    'Jet_bReg[hJetInd2] > 30',
    'Jet_bReg[hJetInd1] > 30',
]
sr_sel = no_sel + [  # copied from selection.dat
    'H_pt>100',
    'H_mass_fit_fallback>90',
    'H_mass_fit_fallback<150',
    'Jet_btagDeepB[hJetInd2]>0.1522',  # FIXME should be applied in VHbbAanlysis.cc
    # 'H_mass>90',
    # 'H_mass<150',
]
sr_sel_str = '&&'.join(sr_sel)

# CS 2lep
cs_2lep_sel = no_sel + [
# no cuts for all CS
]
cs_sel_str = '&&'.join(cs_2lep_sel)


new_VpLF_CR = 'H_mass_fit_fallback > 90 && H_mass_fit_fallback <= 150 && (controlSample==22 || controlSample==24)'
new_VpHF_CR = '(H_mass_fit_fallback < 90 || H_mass_fit_fallback >= 150) && (controlSample==23 || controlSample==25 || \
    (controlSample==0 && Jet_btagDeepB[hJetInd1] > 0.8001 && Jet_btagDeepB[hJetInd2] > 0.1522 && MET_Pt < 60 && HVdPhi > 2.5 && V_mass > 85 && V_mass <= 97))'

cats = {
# SR
    'SR_2Lepton_lowVpt_Mu':       sr_sel_str+'&& controlSample==0 && isZmm==1 %s' % low_vpt,
    'SR_2Lepton_highVpt_Mu':      sr_sel_str+'&& controlSample==0 && isZmm==1 %s' % high_vpt,
    'SR_2Lepton_lowVpt_El':       sr_sel_str+'&& controlSample==0 && isZee==1 %s' % low_vpt,
    'SR_2Lepton_highVpt_El':      sr_sel_str+'&& controlSample==0 && isZee==1 %s' % high_vpt,
# CR low Vpt
    'CR_2Lepton_lowVpt_TTbar_Mu':  cs_sel_str+'&& controlSample==21 && isZmm==1%s    ' % low_vpt,
    'CR_2Lepton_lowVpt_VpLF_Mu':   cs_sel_str+'&& %s                && isZmm==1%s    ' % (new_VpLF_CR, low_vpt),
    'CR_2Lepton_lowVpt_VpHF_Mu':   cs_sel_str+'&& %s                && isZmm==1%s    ' % (new_VpHF_CR, low_vpt),
    'CR_2Lepton_lowVpt_TTbar_El':  cs_sel_str+'&& controlSample==21 && isZee==1%s%s%s' % (good_ele_1, good_ele_2, low_vpt),
    'CR_2Lepton_lowVpt_VpLF_El':   cs_sel_str+'&& %s                && isZee==1%s%s%s' % (new_VpLF_CR, good_ele_1, good_ele_2, low_vpt),
    'CR_2Lepton_lowVpt_VpHF_El':   cs_sel_str+'&& %s                && isZee==1%s%s%s' % (new_VpHF_CR, good_ele_1, good_ele_2, low_vpt),
# CR high Vpt
    'CR_2Lepton_highVpt_TTbar_Mu': cs_sel_str+'&& controlSample==21 && isZmm==1%s    ' % high_vpt,
    'CR_2Lepton_highVpt_VpLF_Mu':  cs_sel_str+'&& %s                && isZmm==1%s    ' % (new_VpLF_CR, high_vpt),
    'CR_2Lepton_highVpt_VpHF_Mu':  cs_sel_str+'&& %s                && isZmm==1%s    ' % (new_VpHF_CR, high_vpt),
    'CR_2Lepton_highVpt_TTbar_El': cs_sel_str+'&& controlSample==21 && isZee==1%s%s%s' % (good_ele_1, good_ele_2, high_vpt),
    'CR_2Lepton_highVpt_VpLF_El':  cs_sel_str+'&& %s                && isZee==1%s%s%s' % (new_VpLF_CR, good_ele_1, good_ele_2, high_vpt),
    'CR_2Lepton_highVpt_VpHF_El':  cs_sel_str+'&& %s                && isZee==1%s%s%s' % (new_VpHF_CR, good_ele_1, good_ele_2, high_vpt),
}

# this dictionary is used for setting up the plotting tools
the_category_dict = {
#   region block     categories         sel.         histograms
    'AllCats':      [cats,              [],          nano_2lep_plotvariables.vars_2lep],
}
