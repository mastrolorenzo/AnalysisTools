
name = 'VHbbPlots'
input_pattern = '/nfs/dust/cms/user/tholenhe/VHbbAnalysisNtuples/V25_VHbb_runonskim_20180212/haddjobs/sum_%s.root'
weight = 'weight*puWeight*sign(genWeight)'
enable_reuse_step = True  # try to find output on disk and don't run a step if present


from main_samples import the_samples_dict, sample_colors
the_samples_dict = dict(
    (sample, the_samples_dict[sample])
    for sample in ['TT_powheg', 'ZH125']
)

# from main_selections import the_category_dict
from main_selections import *
the_category_dict = {
    'SR': [cats_SR, sr_sel, main_plotvariables.vars_used_in_selections],
}





# Howto manipulate single histograms while plotting??
#####################################################


import varial

# a single function should be patched with '@varial.history.track_history', so that
# it appears in the history in the webcreator (rrrreeeally helps with debugging).
@varial.history.track_history
def set_y_min(w):

    # do your thing here .. (setting y_min for single histograms)
    if w.name == 'BDT':
        w.val_y_min = 50 # 1e-2
    elif w.name == 'HVdPhi':
        w.val_y_min = 1e-1

    # the histogram wrapper _must_ be returned
    return w

# then any single-histogram function should be plugged into this is a generator function
# it must be called "additional_input_hook", so that
def additional_input_hook(wrps):

    # this is how to make a generator from a single-argument-function
    # don't do []-brackets, as this will make a list, which slower and wastes memory
    wrps = (set_y_min(w) for w in wrps)

    # could add more functions here, e.g. you could print the names of the histograms flying by..
    @varial.history.track_history
    def print_histo_name(w):
        print w.name
        return w

    wrps = (print_histo_name(w) for w in wrps)

    # if no list or generator is returned, it will complain. Please return.
    return wrps
