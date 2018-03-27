
name = 'VHbbPlotsNanoAOD'
#input_pattern = '/nfs/dust/cms/user/tholenhe/VHbbAnalysisNtuples/V25_VHbb_runonskim_20180121/haddjobs/sum_%s.root'
#input_pattern = '/eos/uscms/store/user/sbc01/VHbbAnalysisNtuples/Nano2016V1_March14_v5/haddjobs/sum_%s.root'
input_pattern = '/eos/uscms/store/group/lpchbb/VHbbAnalysisNtuples/2017V2_March22/%s/*.root'
weight = 'weight'
enable_reuse_step = True  # try to find output on disk and don't run a step if present
treename = 'Events'


from nano_samples import the_samples_dict, sample_colors
from nano_selections import the_category_dict
