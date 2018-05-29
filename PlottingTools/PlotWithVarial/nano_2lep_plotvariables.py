from nano_plotvariables import *

vars_2lep_only = {
    # these plot definitions have different binning and overwrite the above for the 2-lepton channel
    'Vpt':                  ('V_pt',        ';p_{T} (V) [GeV];',                            40,      0,      400    ),
    'Hpt':                  ('H_pt',        ';Regressed p_{T} (jj) [GeV];',                 40,      0,      400    ),
    'DeepCSVmax':           (DeepCSVmax,       ';DeepCSV_{max};',                           20,      0.,     1      ),
    'DeepCSVmin':           (DeepCSVmin,       ';DeepCSV_{min};',                           20,      0.,     1      ),
    'CMVAmax':              (CMVAmax,       ';CMVA_{max};',                                 20,      -1.,    1      ),
    'CMVAmin':              (CMVAmin,       ';CMVA_{min};',                                 20,      -1.,    1      ),
    'HVdPhi':               ('HVdPhi',      ';HVdPhi [rad];',                               30,      -3.2,   3.2    ),
    'jjWPtBal':             ('H_pt/V_pt',   ';p_{T} balance after regression;',             25,      0,      2.     ),
    'drjj':                 ('HJ1_HJ2_dR',  ';reg. Delta R(jj);',                           30,      0,      6      ),
    'Vmass':                ('V_mass',                              ';M_{ll} GeV;',         50,      0,      250     ),

    'muon1Pt':              ('Muon_pt[lepInd1]',              ';1st Muon p_{T};',           30,      0,      150     ),
    'muon1Phi':             ('Muon_phi[lepInd1]',             ';1st Muon #phi;',            30,      -3.2,   3.2     ),
    'muon1Eta':             ('Muon_eta[lepInd1]',             ';1st Muon #eta;',            30,      -3,     3       ),
    #'muon1ID':              ('Muon_looseIdPOG[lepInd1]',      ';1st Muon mu ID;',           6,      -1.5,    4.5     ),

    'muon2Pt':              ('Muon_pt[lepInd2]',              ';2nd Muon p_{T};',           30,      0,      150     ),
    'muon2Phi':             ('Muon_phi[lepInd2]',             ';2nd Muon #phi;',            30,      -3.2,   3.2     ),
    'muon2Eta':             ('Muon_eta[lepInd2]',             ';2nd Muon #eta;',            30,      -3,     3       ),
    #'muon2ID':              ('Muon_looseIdPOG[lepInd2]',      ';2nd Muon mu ID;',           6,      -1.5,    4.5     ),

    'electron1Pt':          ('Electron_pt[lepInd1]',              ';1st Electron p_{T};',   30,      0,      150     ),
    'electron1Phi':         ('Electron_phi[lepInd1]',             ';1st Electron #phi;',    30,      -3.2,   3.2     ),
    'electron1Eta':         ('Electron_eta[lepInd1]',             ';1st Electron #eta;',    30,      -3,     3       ),
    #'electron1ID':          ('Electron_looseIdPOG[lepInd1]',      ';1st Electron mu ID;',   6,      -1.5,    4.5     ),

    'electron2Pt':          ('Electron_pt[lepInd2]',              ';2nd Electron p_{T};',   30,      0,      150     ),
    'electron2Phi':         ('Electron_phi[lepInd2]',             ';2nd Electron #phi;',    30,      -3.2,   3.2     ),
    'electron2Eta':         ('Electron_eta[lepInd2]',             ';2nd Electron #eta;',    30,      -3,     3       ),
    #'electron2ID':          ('Electron_looseIdPOG[lepInd2]',      ';2nd Electron mu ID;',   6,      -1.5,    4.5     ),
    'BDTZllHigh':           ('CMS_vhbb_BDT_Zll_HighPT_13TeV',    ';BDT Output;',            30,      -1,     1       ),
    'BDTZllLow':            ('CMS_vhbb_BDT_Zll_LowPT_13TeV',   ';BDT Output;',              30,      -1,     1       ),

    'Vpt_fit':              ('V_pt_fit',        ';fitted p_{T}(V) GeV;',                    21,      0,      350     ),
    'Hpt_fit':              ('H_pt_fit',        ';fitted p_{T}(jj) GeV;',                   21,      0,      350     ),
    'Hmass_fit':            ('H_mass_fit',      ';fitted M_{jj} GeV;',                      25,      0,      250     ),
    'Hmass_fit_fallback':   ('H_mass_fit_fallback',';fitted M_{jj} GeV;',                   25,      0,      250     ),
    'HVdPhi_fit':           ('abs(HVdPhi_fit)', ';fitted HVdPhi;',                          16,      0,      3.2     ),
    'jjVPtBal_fit':         ('H_pt_fit/V_pt',   ';fitted jj/V P_{T} Bal;',                  30,      0,      2       ),
    'jjVPtBal_fit_fit':     ('H_pt_fit/V_pt_fit', ';fitted jj/ fitted V P_{T} Bal;',        30,      0,      2       ),
    'n_recoil_jets_fit':    ('n_recoil_jets_fit', ';N(recoil jets);',                       12,      -1.5,   10.5    ),
    'llbb_fit_pt':          ('llbb_fit_pt',     ';llbb system p_{T};',                      30,      0,      120     ),
    'llbbr_fit_pt':         ('llbbr_fit_pt',    ';llbbr system p_{T};',                     30,      0,      120     ),

    'kinfit_fit':            ('kinfit_fit',      ';fit result (0:N/A, 1:OK, 2:ERR);',       3,     -.5,      2.5     ),
}
vars_2lep = vars_used_in_selections.copy()
vars_2lep.update(vars_2lep_only)
