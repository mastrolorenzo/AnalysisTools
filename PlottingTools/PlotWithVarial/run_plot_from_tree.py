#!/usr/bin/env python

try:
    import varial.tools
except ImportError:
    print 'You need Varial for this script.'
    print 'Checkout https://github.com/HeinerTholen/Varial'
    exit(-1)

from varial_ext.treeprojector import TreeProjector
from varial.util import iterableize
import collections
import itertools
import glob
import sys
import os

if len(sys.argv) < 2:
    print 'Please provide a config file. Exit.'
    exit(-1)

sys.dont_write_bytecode = True

cwd = os.getcwd()
full_cfg_path = sys.argv[1]
print 'cwd,fukll_cfg: ',cwd, full_cfg_path
cfg_path, cfg = os.path.split(full_cfg_path.replace('.py', ''))
print 'cfg_path after removing of .py', cfg_path
try:
    os.chdir(cfg_path)
    config = __import__(cfg)
    print 'trying...'
except ImportError, e:
    print 'Could not load config. Exit.'
    print e
    exit(-1)
finally:
    os.chdir(cwd)

# use ~80% of cores ...
varial.settings.max_num_processes = int(varial.settings.max_num_processes * 0.8) + 1
varial.raise_root_error_level()

# gather input tokens from sample definitions
all_input_tokens_and_file_patterns = dict(
    (tok, config.input_pattern % tok)
    for items in config.the_samples_dict.itervalues()
    for tok in items[-1]
)
print all_input_tokens_and_file_patterns
filenames = dict(
    (
        s_name,
        list(  # collect all files from all input tokens for a given sample here
            filename
            for input_token in s_items[-1]
            for filename in glob.glob(all_input_tokens_and_file_patterns[input_token])        
        )           
    )
    for s_name, s_items in config.the_samples_dict.iteritems()
)
print 'filename: ', filenames
   
# structure of filenames: {'samplename': ['/path/file1.root', ...]}

# order filenames by size of the sample, so that the largest samples are started first
ordered_s_names = sorted(
    filenames.keys(),
    key=lambda s: sum(os.path.getsize(f) for f in filenames[s]),
    reverse = True,
)
filenames = collections.OrderedDict(
    (s_name, filenames[s_name])
    for s_name in ordered_s_names
)

def make_weight_string(index):
    if index == 0:  # Data
        return '(sampleIndex==0)'
    if index != None:
        return config.weight + '*(sampleIndex==%i)' % index
    return config.weight

# make sample weights_dict
weights_dict = dict(
    (s_name, make_weight_string(s_items[0]))
    for s_name, s_items in config.the_samples_dict.iteritems()
)
# structure of weights_dict: {'samplename': 'weight*(sampleIndex==12345)*scale', ...}


def mk_tree_projector(region_block):
    cats, base_sel, histos = config.the_category_dict[region_block]

    # plug a block of parameter together (including histogram def's)
    params = {
        'histos': histos,
#Luca        'treename': getattr(config, 'treename', 'tree'),
        'treename': getattr(config, 'treename', 'Events'),
        'nm1': region_block.endswith('N-1'),
    }
    params.update(getattr(config, 'params', {}))

    if callable(base_sel):
        params['tree_prep'] = base_sel
        base_sel = []

    # define selections to be run with (section name, selection list, weights)
    nm1_token = '_N-1' if params['nm1'] else ''
    sec_sel_weight = list(
        (cat_name+nm1_token, list(iterableize(cat_sel))+base_sel, weights_dict)
        for cat_name, cat_sel in cats.iteritems()
    )

    # create and feed the tool for the plot-from-tree step
    TP_class = getattr(config, 'TreeProjector', TreeProjector)
    return TP_class(
        filenames, params, sec_sel_weight,
        name=region_block,
        add_aliases_to_analysis=False,
    )

treeplotters = varial.tools.ToolChainParallel(
    'HistosFromTree',
    list(mk_tree_projector(block) for block in config.the_category_dict),
    n_workers=getattr(config, 'n_parallel_treeplotters', 1),
)

# protect against category override:
pltr_combos = itertools.combinations(treeplotters.tool_chain, 2)
for pltr1, pltr2 in pltr_combos:
    sections1 = set(s for s, _, _ in pltr1.sec_sel_weight)
    sections2 = set(s for s, _, _ in pltr2.sec_sel_weight)
    intersect = sections1 & sections2
    if intersect:
        raise RuntimeError(
            'In config "%s": the blocks "%s" and "%s" have common sections: %s' % (
                full_cfg_path, pltr1.name, pltr2.name, intersect))

# plotting...
varial.settings.colors.update(config.sample_colors)

def input_hook(wrps):
    '''Histograms are wrapped. This function sets some meta data.'''
    wrps = varial.gen.gen_add_wrp_info(
        wrps,
        sample = lambda w: os.path.basename(w.file_path)[:-5]
    )
    wrps = varial.gen.gen_add_wrp_info(
        wrps,
        is_data = lambda w: config.the_samples_dict[w.sample][2] == 'Data', # legend == 'Data'?
        is_signal = lambda w: 'H125' in w.sample,
        legend = lambda w: config.the_samples_dict[w.sample][2],
        scale = lambda w: config.the_samples_dict[w.sample][1],
    )

    # scale histograms if the scale is different from 1.
    @varial.history.track_history
    def scale(wrp):
        wrp.obj.Scale(wrp.scale)
        return wrp

    wrps = varial.gen.imap_conditional(wrps, lambda w: w.scale != 1., scale)

    # the cutFlow histogram should show the integral of events up to a certain step
    wrps = varial.gen.imap_conditional(
        wrps,
        lambda w: w.name == 'cutFlow',
        varial.op.int_r,
    )

    return wrps


def mk_plotter_and_webcreator():
    '''
    Instanciate plotter after the plot-from-tree step.

    The plotters cannot be created in the beginning, because the directory structure of the plot
    directory is derived from the files containing the histograms.
    '''
    plotter = varial.tools.mk_rootfile_plotter(
        name='Plots',
        pattern=config.name+'/HistosFromTree/*/*.root',
        hook_loaded_histos=input_hook,
        combine_files=True,
        stack=True,
    )
    return [plotter, varial.tools.WebCreator()]


tc = varial.tools.ToolChain(
    config.name,
    [treeplotters],  # tools to be run immediately
    lazy_eval_tools_func=mk_plotter_and_webcreator,  # to be run later
)


varial.tools.Runner(tc, default_reuse=getattr(config, 'enable_reuse_step', False))
