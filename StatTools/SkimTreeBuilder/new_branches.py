## New Branch Definitions

def is_signal_Znn_I(event):
    """Flag if the event belongs to a 0-lepton signal sample.
    Used by the rebinner.
    """
    # ZH125_ZNuNu_powheg = -12504
    # ggZH125_ZNuNu_powheg = -12505
    # WminusH125_powheg = -12501
    # WplusH125_powheg = -12500
    if event.sampleIndex in {-12504, -12505, -12501, -12500}:
        return 1
    elif event.sampleIndex == 0:
        # These are data events.
        return -1
    else:
        # These are background events.
        return 0

def bin_index_Znn_I(event):
    """Flag if the event falls within the full range of the 0-lepton BDT
    output [-1, 1] in the 0-lepton signal region. Used by the rebinner.
    """
    if event.isZnn and (60 < event.H_mass < 160) and (-1 < event.CMS_vhbb_BDT_Znn_13TeV < 1):
        return 0
    else:
        return -1

def is_signal_Wln_I(event):
    """Flag if the event belongs to a 1-lepton signal sample.
    Used by the rebinner.
    """
    # WminusH125_powheg = -12501
    # WplusH125_powheg = -12500
    if event.sampleIndex in {-12501, -12500}:
        return 1
    elif event.sampleIndex == 0:
        # These are data events.
        return -1
    else:
        # These are background events.
        return 0

def bin_index_Wln_I(event):
    """Flag if the event falls within the full range of the 1-lepton BDT
    output [-1, 1] in the 1-lepton signal region. Used by the rebinner.
    """
    if (event.isWenu or event.isWmunu) and (90 < event.H_mass < 150) and (-1 < event.CMS_vhbb_BDT_Wln_13TeV < 1):
        return 0
    else:
        return -1

def is_signal_Zll_I(event):
    """Flag if the event belongs to a 2-lepton signal sample.
    Used by the rebinner.
    """
    # ZH125_ZLL_powheg = -12502
    # ggZH125_ZLL_powheg = -12503
    if event.sampleIndex in {-12502, -12503}:
        return 1
    elif event.sampleIndex == 0:
        # These are data events.
        return -1
    else:
        # These are background events.
        return 0

def bin_index_Zll_lowPt_I(event):
    """Flag if the event falls within the full range of the 2-lepton BDT
    output [-1, 1] in the 2-lepton low pT(V) signal region. Used by the rebinner.
    """
    if (event.isZee or event.isZmm) and (50 < event.V_pt < 150) and (90 < event.H_mass < 150) and event.hJets_btagged_1 > 0.1522 and (-1 < event.CMS_vhbb_BDT_Zll_LowPT_13TeV < 1):
        return 0
    else:
        return -1

def bin_index_Zll_highPt_I(event):
    """Flag if the event falls within the full range of the 2-lepton BDT
    output [-1, 1] in the 2-lepton high pT(V) signal region. Used by the rebinner.
    """
    if (event.isZee or event.isZmm) and event.V_pt > 150 and (90 < event.H_mass < 150) and event.hJets_btagged_1 > 0.1522 and (-1 < event.CMS_vhbb_BDT_Zll_HighPT_13TeV < 1):
        return 0
    else:
        return -1

