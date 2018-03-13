
name = 'VHbbPlots'
#Luca input_pattern = '/nfs/dust/cms/user/tholenhe/VHbbAnalysisNtuples/V25_VHbb_runonskim_20180212/haddjobs/sum_%s.root'
#input_pattern = '/nfs/dust/cms/user/lmastrol/VHbbAnalysisNtuples/VHbb_March18_v0/haddjobs/sum_%s.root'
input_pattern = '/nfs/dust/cms/user/lmastrol/VHbbAnalysisNtuples/VHbb_March18_v2/haddjobs/sum_%s.root'
weight = 'weight*puWeight*sign(genWeight)'
enable_reuse_step = True  # try to find output on disk and don't run a step if present

from main_samples import the_samples_dict, sample_colors
the_samples_dict = dict(
    (sample, the_samples_dict[sample])
#Luca    for sample in ['TT_powheg', 'ZH125']
    for sample in ['T_tW','Tbar_tW','W_udcsg_HT100To200','W_udcsg_HT200To400','W_udcsg_HT400To600','W_b_HT100To200','W_b_HT200To400','W_b_HT400To600','W_bb_HT100To200','W_bb_HT200To400','W_bb_HT400To600','ZH125','WH125p','WH125m']

)

# from main_selections import the_category_dict
from main_selections import *
the_category_dict = {
    'SR': [cats_SR, sr_sel, main_plotvariables.vars_used_luca],
}
