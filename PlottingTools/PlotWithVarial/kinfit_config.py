import varial
varial.settings.rootfile_postfixes += ['.pdf']
name = 'VHbbPlotsKinfitV3_2'
input_pattern = '/nfs/dust/cms/user/tholenhe/VHbbAnalysisNtuples/VHbb2016_NanoAODPostProcV3_Direct1/%s/*.root'
#weight = 'genWeight*weight_PU*bTagWeight*weight_ptQCD*weight_ptEWK*Lep_SF*weight_mettrigSF'
weight = 'weight'
enable_reuse_step = True  # try to find output on disk and don't run a step if present
treename = 'Events'


from kinfit_samples import *


plotvars = {
    'Vpt':              ('V_pt',            ';p_{T}(V) GeV;',               21,      0,      350     ),
    'Vpt_fit':          ('V_pt_fit',        ';fitted p_{T}(V) GeV;',        21,      0,      350     ),
    'Hpt':              ('H_pt',            ';p_{T}(jj) GeV;',              21,      0,      350     ),
    'Hpt_fit':          ('H_pt_fit',        ';fitted p_{T}(jj) GeV;',       21,      0,      350     ),
    'Hmass':            ('H_mass',          ';M_{jj} GeV;',                 25,      0,      250     ),
    'Hmass_fit':        ('H_mass_fit',      ';fitted M_{jj} GeV;',          25,      0,      250     ),
    'Hmass_fit_fallback':('H_mass_fit_fallback',';fitted M_{jj} GeV;',      25,      0,      250     ),
    'HVdPhi':           ('abs(HVdPhi)',     ';HVdPhi;',                     16,      0,      3.2     ),
    'HVdPhi_fit':       ('abs(HVdPhi_fit)', ';fitted HVdPhi;',              16,      0,      3.2     ),
    'jjVPtBal':         ('H_pt/V_pt',       ';jj/V P_{T} Bal;',             30,      0,      2       ),
    'jjVPtBal_fit':     ('H_pt_fit/V_pt',   ';fitted jj/V P_{T} Bal;',      30,      0,      2       ),

    'kinfit_fit':       ('kinfit_fit',      ';fit result (0:N/A, 1:OK, 2:ERR);',3, -.5,      2.5     ),
    'controlSample':    ('controlSample',   ';controlSample;',              27,   -1.5,      25.5    ),
}


high_vpt = ' && V_pt >= 150'
good_ele_1 = ' && (abs(Electron_eta[lepInd1])<1.44 || abs(Electron_eta[lepInd1])>1.56)'
good_ele_2 = ' && (abs(Electron_eta[lepInd2])<1.44 || abs(Electron_eta[lepInd2])>1.56)'
cats_SR = {
    'SR_2Lepton_highVpt_Mu':    '(controlSample==0 && isZmm==1 %s)' % high_vpt,
    'SR_2Lepton_highVpt_El':    '(controlSample==0 && isZee==1 %s)' % high_vpt,
}

new_VpLF_CR = 'H_mass_fit_fallback > 90 && H_mass_fit_fallback <= 150 && (controlSample==22 || controlSample==24)'
new_VpHF_CR = '(H_mass_fit_fallback < 90 || H_mass_fit_fallback >= 150) && (controlSample==23 || controlSample==25 || \
    (controlSample==0 && Jet_btagCMVA[hJetInd1] > 0.9432 && Jet_btagCMVA[hJetInd2] > -0.5884 && MET_Pt < 60 && HVdPhi > 2.5 && V_mass > 85 && V_mass <= 97))'

# type=2  name=tagWPL             val=-0.5884
# type=2  name=tagWPM             val=0.4432
# type=2  name=tagWPT             val=0.9432


cats_CR = {
    'CR_2Lepton_highVpt_TTbar_Mu': '(controlSample==21 && isZmm==1%s)' % high_vpt,
    'CR_2Lepton_highVpt_VpLF_Mu':  '(%s && isZmm==1%s)' % (new_VpLF_CR, high_vpt),
    'CR_2Lepton_highVpt_VpHF_Mu':  '(%s && isZmm==1%s)' % (new_VpHF_CR, high_vpt),
    'CR_2Lepton_highVpt_TTbar_El': '(controlSample==21 && isZee==1%s%s%s)' % (high_vpt, good_ele_1, good_ele_2),
    'CR_2Lepton_highVpt_VpLF_El':  '(%s && isZee==1%s%s%s)' % (new_VpLF_CR, high_vpt, good_ele_1, good_ele_2),
    'CR_2Lepton_highVpt_VpHF_El':  '(%s && isZee==1%s%s%s)' % (new_VpHF_CR, high_vpt, good_ele_1, good_ele_2),
}

cats_debug = {
    'Debug_isZmm': 'isZmm',
    'Debug_isZee': 'isZee',
}

no_sel = [
    'Pass_nominal',
    'twoResolvedJets',
    'cutFlow>=2',
    'Jet_bReg[hJetInd2] > 30',
    'Jet_bReg[hJetInd1] > 30',
    #'Jet_Pt[hJetInd1] > 20',
    #'Jet_Pt[hJetInd2] > 20',
]
sr_sel = no_sel + [  # copied from selection.dat
    'H_pt>100',
    'H_mass_fit_fallback>90',
    'H_mass_fit_fallback<150',
]

the_category_dict = {
    'SR':           [cats_SR,    sr_sel,      plotvars],
    'CR':           [cats_CR,    no_sel,      plotvars],
    # 'Debug':        [cats_debug, [],          plotvars],
}

from varial_ext.treeprojector import TreeProjectorFileBased
TreeProjector = TreeProjectorFileBased
