

muon_pt_safe = 'lepInd1 > -1 && nMuon > lepInd1 ? Muon_pt[lepInd1] : -1.'
electron_pt_safe = 'lepInd1 > -1 && nElectron > lepInd1 ? Electron_pt[lepInd1] : -1.'
DeepCSVmax = 'Jet_btagDeepB[hJetInd1] > Jet_btagDeepB[hJetInd2] ? Jet_btagDeepB[hJetInd1] : Jet_btagDeepB[hJetInd2]'
DeepCSVmin = 'Jet_btagDeepB[hJetInd1] > Jet_btagDeepB[hJetInd2] ? Jet_btagDeepB[hJetInd2] : Jet_btagDeepB[hJetInd1]'
CMVAmax = 'Jet_btagCMVA[hJetInd1] > Jet_btagCMVA[hJetInd2] ? Jet_btagCMVA[hJetInd1] : Jet_btagCMVA[hJetInd2]'
CMVAmin = 'Jet_btagCMVA[hJetInd1] > Jet_btagCMVA[hJetInd2] ? Jet_btagCMVA[hJetInd2] : Jet_btagCMVA[hJetInd1]'
lepRelIso1_safe = 'isWmunu || isWenu || isZmm || isZee ? selLeptons_relIso_0 : -1'
lepRelIso2_safe = 'isZmm || isZee ? selLeptons_relIso_1 : -1'


test = {
    
        'V_pt':     ("V_pt",                'test;p_{T}(V)@[GeV];', 30, 200, 500),
        'FatJetCand_pt':     ("FatJetCand_pt", 'FatJetCand_pt;p_{T}(fj)@[GeV];', 20, 200, 500),#was 12 bins
        'FatJetCand_eta':     ("FatJetCand_eta", 'FatJetCand_eta;\eta(fj);', 40, -3, 3),
        'FatJetVdPhi':      ("FatJetVdPhi;\Delta\phi(V,FatJetCand);", 30, 2.85, 3.15),
        'V_mt':     ("V_mt",     'V_mt', 30, 0, 150),
        'nBJetsOutsideFatJet':  ("nBJetsOutsideFatJet",             "nBJetsOutsideFatJet", 10, 0, 10),
        'FatJetCand_tau21':     ("FatJetCand_tau21",                'FatJetCand_tau21;#tau_{21};', 20, 000, 1),
        'FatJetCand_tau32':     ("FatJetCand_tau32",                'FatJetCand_tau32;#tau_{32};', 20, 000, 1),
        'FatJetCand_doubleB':     ("FatJetCand_doubleB",                'FatJetCand_doubleB;double btag;', 15, .85, 1),
        'FatJetCand_Msoftdrop_corrected':     ("FatJetCand_Msoftdrop_corrected",      'FatJetCand_Msoftdrop_corrected;m_{SD}@[GeV];', 15, 50, 200),
        #'FatJetCand_mass':  ('FatJetCand_mass', 'FatJetCand_mass', 20, 0, 200),
        'oneMergedJet': ('oneMergedJet', 'oneMergedJet', 2, 0, 2),
        'OutsideJetMaxPt' : ('Max$(Jet_Pt*(sqrt(pow(FatJetCand_eta-Jet_eta,2)+pow((FatJetCand_phi-Jet_phi)*(fabs(FatJetCand_phi-Jet_phi)<3.14)+(2*3.14-fabs(FatJet_phi-Jet_phi))*(fabs(FatJet_phi-Jet_phi)>3.14),2))>0.8))' , 'OutsideJetMaxPt', 50, 100, 400),#Max pt of jet with dR(jet,fatjetcand)>0.8
    #'met_pt':           ('MET_Pt',                              ';MET@p_{T}@GeV;',              50,      0,      500     ),
    #'Hpt':              ('H_pt',                                ';p_{T}(jj)@GeV;',              21,      0,      350     ),
    #'Vtype':            ('Vtype',                               ';Vtype;',                      10,      -.5,    9.5     ),
    #'jetsubleadpt':     ('Jet_bReg[hJetInd2]',                  ';Regressed Jet2@P_{T}@GeV;',             25,      0,      250     ),
    #'DeepCSVmin':          (DeepCSVmin,                               ';DeepCSV_{min};',                 60,      0.,     1       ),
   # 'DeepCSVmin':          (DeepCSVmin,                               ';DeepCSV_{min};',                 20,      0.,     1       ),
    #'BDTZllHigh':              ('CMS_vhbb_BDT_Zll_HighPT_13TeV',              ';BDT@Output;',                 30,      -1,     1       ),
    'GenFatJetCand_partonFlavour': ('GenFatJetCand_partonFlavour', 'GenFatJetCand_partonFlavour', 7,-0.5,6.5),
    
    #'GenJetAK8_partonFlavour': ('GenJetAK8_partonFlavour[hGenFatJetCand]', 'GenJetAK8_partonFlavour', 6, 0, 6),
    #'GenJetAK8_pt': ('GenJetAK8_pt[hGenFatJetCand]', 'GenJetAK8_pt', 20, 0, 500),
    #'GenJetAK8_eta': ('GenJetAK8_eta[hGenFatJetCand]', 'GenJetAK8_eta', 20, -3, 3),
    #'GenJetAK8_phi': ('GenJetAK8_phi[hGenFatJetCand]', 'GenJetAK8_phi', 20, -3.14,3.14),

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
    # MetjDPhi not yet in tree
    # 'MetTkMetDPhi':     ('MetTkMetDPhi',                        ';MetTkMetDPhi;',               16,      0,      3.2     ),
    'lepMetDPhi':       ('lepMetDPhi',                          ';lepMetDPhi;',                 16,      0,      3.2     ),
    #'lepRelIso1':       (lepRelIso1_safe,                       ';1st Lep@Rel@Iso;',            30,      0,      0.3     ),
    #'lepRelIso2':       (lepRelIso2_safe,                       ';2st Lep@Rel@Iso;',            30,      0,      0.3     ),
    'BDT':              ('CMS_vhbb_BDT_Wln_13TeV',              ';BDT@Output;',                 30,      -1,     1       ),

# additional from 0-lepton CR selection
    'Vtype':            ('Vtype',                               ';Vtype;',                      10,      -.5,    9.5     ),
    # 'minMetjDPhi':      ('minMetjDPhi',                         ';MetjDPhi;',                   16,      0,      3.2     ),

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

