

muon_pt_safe = 'lepInd1 > -1 && nMuon > lepInd1 ? Muon_pt[lepInd1] : -1.'
electron_pt_safe = 'lepInd1 > -1 && nElectron > lepInd1 ? Electron_pt[lepInd1] : -1.'
DeepCSVmax = 'Jet_btagDeepB[hJetInd1] > Jet_btagDeepB[hJetInd2] ? Jet_btagDeepB[hJetInd1] : Jet_btagDeepB[hJetInd2]'
DeepCSVmin = 'Jet_btagDeepB[hJetInd1] > Jet_btagDeepB[hJetInd2] ? Jet_btagDeepB[hJetInd2] : Jet_btagDeepB[hJetInd1]'
CMVAmax = 'Jet_btagCMVA[hJetInd1] > Jet_btagCMVA[hJetInd2] ? Jet_btagCMVA[hJetInd1] : Jet_btagCMVA[hJetInd2]'
CMVAmin = 'Jet_btagCMVA[hJetInd1] > Jet_btagCMVA[hJetInd2] ? Jet_btagCMVA[hJetInd2] : Jet_btagCMVA[hJetInd1]'
lepRelIso1_safe = 'isWmunu || isWenu || isZmm || isZee ? selLeptons_relIso_0 : -1'
lepRelIso2_safe = 'isZmm || isZee ? selLeptons_relIso_1 : -1'

test = {
    #'met_pt':           ('MET_Pt',                              ';MET@p_{T}@GeV;',              50,      0,      500     ),
    #'Hpt':              ('H_pt',                                ';p_{T}(jj)@GeV;',              21,      0,      350     ),
    #'Vtype':            ('Vtype',                               ';Vtype;',                      10,      -.5,    9.5     ),
    #'jetsubleadpt':     ('Jet_bReg[hJetInd2]',                  ';Regressed Jet2@P_{T}@GeV;',             25,      0,      250     ),
    'DeepCSVmin':          (DeepCSVmin,                               ';DeepCSV_{min};',                 60,      0.,     1       ),
   # 'DeepCSVmin':          (DeepCSVmin,                               ';DeepCSV_{min};',                 20,      0.,     1       ),
    #'BDTZllHigh':              ('CMS_vhbb_BDT_Zll_HighPT_13TeV',              ';BDT@Output;',                 30,      -1,     1       ),
}

vars_used_in_selections = {
#   histoname            variable                               title;x-axis;y-axis             n-bins   low     up
# signal region
    'Vpt':              ('V_pt',                                ';p_{T}(W)@[GeV];',             30,      0,      500     ),
    'muonPt':           (muon_pt_safe,                          ';Muon@p_{T};',                 30,      0,      300     ),
    'electronPt':       (electron_pt_safe,                      ';Electron@p_{T};',             30,      0,      300     ),
    'jetleadpt':        ('Jet_bReg[hJetInd1]',                  ';Regressed Jet1@P_{T}@GeV;',             50,      0,      500     ),
    'jetsubleadpt':     ('Jet_bReg[hJetInd2]',                  ';Regressed Jet2@P_{T}@GeV;',             25,      0,      250     ),
    'jetleadptNR':        ('Jet_Pt[hJetInd1]',                  ';Jet1@P_{T}@GeV;',             30,      0,      500     ),
    'jetsubleadptNR':     ('Jet_Pt[hJetInd2]',                  ';Jet2@P_{T}@GeV;',             25,      0,      250     ),
    'leadjeteta':       ('Jet_eta[hJetInd1]',                   ';Jet1#eta;',                   30,      -3,     3       ),
    'subjeteta':        ('Jet_eta[hJetInd2]',                   ';Jet2#eta;',                   30,      -3,     3       ),
    'Hpt':              ('H_pt',                                ';p_{T}(jj)@GeV;',              21,      0,      350     ),
    'Hmass':            ('H_mass',                              ';M_{jj}@GeV;',                 25,      0,      250     ),
    'HVdPhi':           ('HVdPhi',                              ';HVdPhi;',                     16,      0,      3.2     ),
    'DeepCSVmax':          (DeepCSVmax,                               ';DeepCSV_{max};',                 60,      0.,     1       ),
    'DeepCSVmin':          (DeepCSVmin,                               ';DeepCSV_{min};',                 60,      0.,     1       ),
    'leadjetDeepCSVZoom':  ('Jet_btagDeepB[hJetInd1]',              ';Jet1@DeepCSV;',                  20,      0.94,   1       ),
    'CMVAmax':          (CMVAmax,                               ';CMVA_{max};',                 60,      -1.,     1       ),
    'CMVAmin':          (CMVAmin,                               ';CMVA_{min};',                 60,      -1.,     1       ),
    'leadjetCMVAZoom':  ('Jet_btagCMVA[hJetInd1]',              ';Jet1@CMVA;',                  20,      0.94,   1       ),
    'nAddJets':         ('nAddJets302p5_puid',                  ';nAddJets302p5_puid;',         11,      -.5,    10.5    ),
    'nAddLeptons':      ('nAddLeptons',                         ';nAddLeptons;',                11,      -.5,    10.5    ),
    'met_pt':           ('MET_Pt',                              ';MET@p_{T}@GeV;',              30,      0,      500     ),
    'dPhi_MET_TkMET':    ('dPhi_MET_TkMET',                        ';MetTkMetDPhi;',               16,      0,      3.2     ),
    'lepMetDPhi':       ('lepMetDPhi',                          ';lepMetDPhi;',                 16,      0,      3.2     ),
    #'lepRelIso1':       (lepRelIso1_safe,                       ';1st Lep@Rel@Iso;',            30,      0,      0.3     ),
    #'lepRelIso2':       (lepRelIso2_safe,                       ';2st Lep@Rel@Iso;',            30,      0,      0.3     ),
    'BDT':              ('CMS_vhbb_BDT_Wln_13TeV',              ';BDT@Output;',                 30,      -1,     1       ),

# additional from 0-lepton CR selection
    'Vtype':            ('Vtype',                               ';Vtype;',                      10,      -.5,    9.5     ),
    'min_dPhi_hJet_MET':      ('min_dPhi_hJet_MET',                         ';MetjDPhi;',                   16,      0,      3.2     ),

# additional from 1-lepton CR selection
    'sigma_met_pt':     ('MET_Pt / sqrt(htJet30)',              ';#sigma(p_{T}^{miss});',       20,      0,      10      ),

# other
    'cutFlow':          ('cutFlow',                             ';cutFlow;',                    15,      -.5,    14.5    ),
    'controlSample':    ('controlSample',                       ';controlSample;',              30,      -1.5,   28.5    ),
    #'nPU':              ('Pileup_nPU',                          ';nPU;',                        50,      0,      50     ),
    'rho':              ('fixedGridRhoFastjetAll',              ';#rho;',                       60,      0,      60     ),
}


vars_1lepCR_only = {  # copied from plotvariables_CS.dat
 # BJet variables
    #'leadjetptNR':      ('Jet_pt[hJetInd1]',                    ';Jet1@P_{T}@GeV;',             30,      0,      500     ),
    #'leadjetpt':        ('Jet_bReg[hJetInd1]',                  ';Jet1@P_{T}@GeV;',             30,      0,      500     ),
    #'leadjeteta':       ('Jet_eta[hJetInd1]',                   ';Jet1#eta;',                   30,      -3,     3       ),
    #'leadjetDeepCSV':      ('Jet_btagDeepB[hJetInd1]',              ';Jet1@DeepCSV;',                  30,      0.,     1       ),
    #'leadjetDeepCSVZoom':  ('Jet_btagDeepB[hJetInd1]',              ';Jet1@DeepCSV;',                  20,      0.94,   1       ),
    #'subjetpt':         ('Jet_bReg[hJetInd2]',                  ';Jet2@P_{T}@GeV;',             25,      0,      250     ),
    #'subjetptNR':       ('Jet_pt[hJetInd2]',                    ';Jet2@P_{T}@GeV;',             25,      0,      250     ),
    #'subjeteta':        ('Jet_eta[hJetInd2]',                   ';Jet2#eta;',                   30,      -3,     3       ),
    #'subjetDeepCSV':       ('Jet_btagDeepB[hJetInd2]',              ';Jet2@DeepCSV;',                  30,      0.,     1       ),

# lepton
    #'lepRelIso':        ('selLeptons_relIso_0',                 ';Lep@Rel@Iso;',                20,      0,      0.2     ),
    'muonPt':           (muon_pt_safe,                          ';Muon@p_{T};',                 30,      0,      300     ),
    'muonPhi':          (muon_pt_safe.replace("_pt","_phi"),    ';Muon@#phi;',                  30,      -3.2,   3.2     ),
    'muonEta':          (muon_pt_safe.replace("_pt","_eta"),     ';Muon@#eta;',                  30,      -3,     3       ),
    'electronPt':       (electron_pt_safe,                      ';Electron@p_{T};',             30,      0,      300     ),
    'electronPhi':      (electron_pt_safe.replace("_pt","_phi"),';Electron@#phi;',              30,      -3.2,   3.2     ),
    'electronEta':      (electron_pt_safe.replace("_pt","_eta"), ';Electron@#eta;',              30,      -3,     3       ),

# DiJet Variables
    'mjj':              ('H_mass',                              ';M_{jj}@GeV;',                 25,      0,      250     ),
    'mjj_noreg':        ('H_mass_noreg',                        ';M_{jj}@NoReg;',               25,      0,      250     ),
    'Hpt':              ('H_pt',                                ';p_{T}(jj)@GeV;',              21,      0,      350     ),
    'dphijj':           ('HJ1_HJ2_dPhi',                        ';dPhi(jj);',                   20,      -3,     3       ),
    'detajj':           ('HJ1_HJ2_dEta',                        ';dEta(jj);',                   20,      0,      3       ),
    'drjj':             ('HJ1_HJ2_dR',                          ';dR(jj);',                     20,      0,      5       ),
    'jjWPtBal':         ('H_pt/V_pt',                           ';jj/W@P_{T}@Bal;',             30,      0,      2       ),
    'HVdPhi':           ('HVdPhi',                              ';HVdPhi;',                     16,      0,      3.2     ),

# MET
    'met_pt':           ('MET_Pt',                              ';MET@p_{T}@GeV;',              30,      0,      500     ),
    'met_phi':          ('MET_Phi',                             ';MET@Phi;',                    30,      -3.15,  3.15    ),
    'lepMetDPhi':       ('lepMetDPhi',                          ';Lep+Met@d#phi;',              30,      0,      3.14    ),
    'sqrtHTJet30':      ('sqrt(htJet30)',                       ';sqrt(htJet30);',              25,      0,      50      ),
    'MetSig':           ('MET_Pt/sqrt(htJet30)',                ';MET@Sig;',                    25,      0,      25      ),
    ##'min_met':          ('min(MET_Pt,mhtJet30)',                ';min(MET,MHT);',               25,      0,      400     ),

# BDT
    'BDT':              ('CMS_vhbb_BDT_Wln_13TeV',              ';BDT@Output;',                 30,      -1,     1       ),

# other
    'TopMass':          ('Top1_mass_fromLepton_regPT_w4MET',    ';m_{top}@[GeV];',              20,      0,      600     ),
    'TopMassNoReg':     ('Top1_mass_fromLepton_w4MET',          ';TopMass@No@Reg;',             20,      0,      600     ),
    'nSAJets5_sel':     ('SA5',                                 ';num.@SA5@Jets;',              10,      0,      10      ),
    'nAddJets':         ('nAddJets252p9_puid',                  ';nAddJets252p9_puid;',         11,      -.5,    10.5    ),
    'Vmt':              ('V_mt',                                ';W@m_{T};',                    20,      0,      200     ),
    'Vpt':              ('V_pt',                                ';p_{T}(W)@[GeV];',             30,      0,      500     ),
}
vars_1lepCR = vars_used_in_selections.copy()
vars_1lepCR.update(vars_1lepCR_only)

vars_2lepCR_only = {
    # these plot definitions have different binning and overwrite the above for the 2-lepton channel
    # 'BDT':              ('CMS_vhbb_BDT_Wln_13TeV',              ';BDT@Output;',                 20,      -1,     1       ),

    'Vpt':              ('V_pt',        ';p_{T} (V) [GeV];',                40,      0,      400    ),
    'Hpt':              ('H_pt',        ';Regressed p_{T} (jj) [GeV];',     40,      0,      400    ),
    'mjj':              ('H_mass',      ';Regressed m(jj) [GeV];',          17,      0,      255    ),
    'DeepCSVmax':          (DeepCSVmax,       ';DeepCSV_{max};',                     20,      0.,     1      ),
    'DeepCSVmin':          (DeepCSVmin,       ';DeepCSV_{min};',                     20,      0.,     1      ),
    'CMVAmax':          (CMVAmax,       ';CMVA_{max};',                     20,      -1.,     1      ),
    'CMVAmin':          (CMVAmin,       ';CMVA_{min};',                     20,      -1.,     1      ),
    'HVdPhi':           ('HVdPhi',      ';HVdPhi [rad];',                   30,      -3.2,   3.2    ),
    'jjWPtBal':         ('H_pt/V_pt',   ';p_{T} balance after regression;', 25,      0,      2.     ),
    'drjj':             ('HJ1_HJ2_dR',  ';reg. Delta R(jj);',               30,      0,      6      ),
    'Vmass':            ('V_mass',                              ';M_{ll}@GeV;',                 50,      0,      250     ),

    'muon1Pt':            ('Muon_pt[lepInd1]',              ';1st Muon@p_{T};',               30,      0,      150     ),
    'muon1Phi':           ('Muon_phi[lepInd1]',             ';1st Muon@#phi;',                30,      -3.2,   3.2     ),
    'muon1Eta':           ('Muon_eta[lepInd1]',             ';1st Muon@#eta;',                30,      -3,     3       ),
    #'muon1ID':          ('Muon_looseIdPOG[lepInd1]',      ';1st Muon mu ID;',               6,      -1.5,    4.5     ),

    'muon2Pt':            ('Muon_pt[lepInd2]',              ';2nd Muon@p_{T};',               30,      0,      150     ),
    'muon2Phi':           ('Muon_phi[lepInd2]',             ';2nd Muon@#phi;',                30,      -3.2,   3.2     ),
    'muon2Eta':           ('Muon_eta[lepInd2]',             ';2nd Muon@#eta;',                30,      -3,     3       ),
    #'muon2ID':          ('Muon_looseIdPOG[lepInd2]',      ';2nd Muon mu ID;',               6,      -1.5,    4.5     ),
    
    'electron1Pt':            ('Electron_pt[lepInd1]',              ';1st Electron@p_{T};',               30,      0,      150     ),
    'electron1Phi':           ('Electron_phi[lepInd1]',             ';1st Electron@#phi;',                30,      -3.2,   3.2     ),
    'electron1Eta':           ('Electron_eta[lepInd1]',             ';1st Electron@#eta;',                30,      -3,     3       ),
    #'electron1ID':          ('Electron_looseIdPOG[lepInd1]',      ';1st Electron mu ID;',               6,      -1.5,    4.5     ),

    'electron2Pt':            ('Electron_pt[lepInd2]',              ';2nd Electron@p_{T};',               30,      0,      150     ),
    'electron2Phi':           ('Electron_phi[lepInd2]',             ';2nd Electron@#phi;',                30,      -3.2,   3.2     ),
    'electron2Eta':           ('Electron_eta[lepInd2]',             ';2nd Electron@#eta;',                30,      -3,     3       ),
    #'electron2ID':          ('Electron_looseIdPOG[lepInd2]',      ';2nd Electron mu ID;',               6,      -1.5,    4.5     ),
    'BDTZllHigh':              ('CMS_vhbb_BDT_Zll_HighPT_13TeV',              ';BDT@Output;',                 30,      -1,     1       ),
    'BDTZllLow':              ('CMS_vhbb_BDT_Zll_LowPT_13TeV',              ';BDT@Output;',                 30,      -1,     1       ),
}
vars_2lepCR = vars_used_in_selections.copy()
vars_2lepCR.update(vars_2lepCR_only)
