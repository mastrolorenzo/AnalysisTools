import varial
varial.settings.rootfile_postfixes += ['.pdf']
name = 'VHbbPlots2Lep2016_2'
input_pattern = '/pnfs/desy.de/cms/tier2/store/user/htholen/VHbbNTuple/2016V4_RunALL_try2/%s/*.root'
weight = 'weight'
enable_reuse_step = True  # try to find output on disk and don't run a step if present
treename = 'Events'

from nano_2lep_samples_2016 import *
from nano_2lep_selections import the_category_dict




# BLINDING
##########

import varial
def additional_input_hook(wrps):

    @varial.history.track_history
    def blind_H_mass_and_BDT_in_SR(w):
        if w.legend == 'Data' and w.in_file_path.startswith('SR_'):
            if w.name.startswith('Hmass'):
                print 'BLINDING %s in Data' % w.in_file_path
                for i in xrange(w.histo.GetNbinsX() + 1):
                    w.histo.SetBinContent(i, 0.)
                    w.histo.SetBinError(i, 0.)

            if w.name.startswith('BDT'):
                print 'BLINDING %s in Data' % w.in_file_path
                for i in xrange(19, w.histo.GetNbinsX() + 1):
                    w.histo.SetBinContent(i, 0.)
                    w.histo.SetBinError(i, 0.)

        return w

    wrps = (blind_H_mass_and_BDT_in_SR(w) for w in wrps)
    return wrps
