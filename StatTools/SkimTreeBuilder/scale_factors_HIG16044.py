## Scale Factors from HIG-16-044

def sf_tt(event):
    """The scale factors for TTbar."""
    if event.isZnn:
        return event.weight * 0.78
    elif event.isWenu or event.isWmunu:
        return event.weight * 0.91
    elif event.isZee or event.isZmm:
        # The scale factor for the 2-lepton channel low V_pt bin is unity.
        if event.V_pt > 150:
            return event.weight * 1.04
        else:
            return event.weight
    else:
        return event.weight

def sf_zjets(event):
    """The scale factors for Z+Jets processes."""
    sampleIndex = str(event.sampleIndex)
    if event.isZnn:
        if sampleIndex.endswith('00'):
            return event.weight * 1.03
        elif sampleIndex.endswith('01'):
            return event.weight * 1.28
        elif sampleIndex.endswith('02'):
            return event.weight * 1.61
        else:
            return event.weight
    elif event.isZee or event.isZmm:
        # Low V_pt Bin
        if 50 < event.V_pt < 150:
            if sampleIndex.endswith('00'):
                return event.weight * 1.01
            elif sampleIndex.endswith('01'):
                return event.weight * 0.98
            elif sampleIndex.endswith('02'):
                return event.weight * 1.09
            else:
                return event.weight
        # High V_pt Bin
        elif event.V_pt > 150:
            if sampleIndex.endswith('00'):
                return event.weight * 1.02
            elif sampleIndex.endswith('01'):
                return event.weight * 1.02
            elif sampleIndex.endswith('02'):
                return event.weight * 1.28
            else:
                return event.weight
        else:
            return event.weight
    else:
        return event.weight

def sf_wjets(event):
    """The scale factors for W+Jets processes."""
    sampleIndex = str(event.sampleIndex)
    if event.isZnn or event.isWenu or event.isWmunu:
        if sampleIndex.endswith('00'):
            return event.weight * 1.14
        elif sampleIndex.endswith('01'):
            return event.weight * 1.66
        elif sampleIndex.endswith('02'):
            return event.weight * 1.49
        else:
            return event.weight
    else:
        return event.weight

