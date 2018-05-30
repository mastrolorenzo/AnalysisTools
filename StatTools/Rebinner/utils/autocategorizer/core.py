import os

import ROOT


MODULE_DIR = os.path.abspath(os.path.dirname(__file__))

# Update interpreter include paths.
ROOT.gInterpreter.AddIncludePath(os.path.join(MODULE_DIR, 'bdt', 'src', 'bdtlib'))
ROOT.gInterpreter.AddIncludePath(os.path.join(MODULE_DIR, 'bdt', 'src', 'tinyxml2'))

# Include autocategorizer headers.
AUTOCAT_HEADERS = [
    'CategoryReader.h',
    'Event.h',
    'Function.h',
    'Node.h',
    'SignificanceMetrics.hxx',
    'Tree.h',
    'Utilities.h',
    'tinyxml2.h',
]

for _ in AUTOCAT_HEADERS:
    ROOT.gInterpreter.ProcessLine('#include "{0}"'.format(_))

# Load the shared object library containing necessary symbols.
ROOT.gSystem.Load(os.path.join(MODULE_DIR, 'libAutocategorizer.so'))

# Load h2mumu specific macros.
ROOT.gSystem.CompileMacro(os.path.join(MODULE_DIR, 'bdt', 'studies', 'h2mumu', 'LoadEvents.hxx'), 'kOs')

