import collections
import glob
import os

import ROOT
import jinja2


MODULE_DIR = os.path.abspath(os.path.dirname(__file__))

TMPDIR = os.path.join(MODULE_DIR, 'tmp')

TEMPLATE_DIR = os.path.join(MODULE_DIR, 'templates')

JINJA_ENV = jinja2.Environment(loader=jinja2.FileSystemLoader(TEMPLATE_DIR), trim_blocks=True)


# Clean up dynamically generated build products upon import.
for pattern in ('*.cxx', '*.d', '*.pcm'):
    for match in glob.glob(os.path.join(TMPDIR, pattern)):
        os.remove(match)


def _get_branchtypes(path, treename, branches):
    """Return a dictionary of branch names mapped to their type."""
    f = ROOT.TFile.Open(path)
    t = f.Get(treename)
    branch_types = collections.OrderedDict()
    for branch in branches:
        branch_types[branch] = t.GetBranch(branch).GetListOfLeaves()[0].GetTypeName()
    f.Close()
    return branch_types


def loadEvents(events, features, src, treename, aliases={}, debug=False):
    """Load events from a ROOT file.

    Parameters
    ----------
    events : autocategorizer.Events
        The vector of rebinning events.
    features : autocategorizer.Features
        The vector of feature branch names.
    src : path
        The path to the ROOT file.
    treename : string
        The name of the TTree from which to load events.
    aliases : dict, optional
        A dictionary mapping the default names of branches required
        by the rebinner to new aliases.
    debug : bool, optional
        If True, the source code of the dynamically generated
        C++ macro is printed for debugging. The default is False.
    """
    # We use the memory address of the features argument
    # because it is guaranteed unique during its lifetime.
    funcname = 'loadEvents_{0!s}'.format(id(features))
    func = globals().get(funcname, None)
    if func is None:
        branchtypes = _get_branchtypes(src, treename, features)
        macro_source = JINJA_ENV.get_template('loadEvents').render(funcname=funcname, branchtypes=branchtypes, aliases=aliases)
        if debug:
            print macro_source
        macro_path = os.path.join(TMPDIR, '{0}.cxx'.format(funcname))
        with open(macro_path, 'wb') as f:
            f.write(macro_source)
        ROOT.gSystem.CompileMacro(macro_path, 'fOs')
        func = getattr(ROOT, funcname)
        globals()[funcname] = func
    return func(events, features, src, treename)

