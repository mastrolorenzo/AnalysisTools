from .core import ROOT as _ROOT
from .load_events import loadEvents


# Classes from SignificanceMetrics.hxx
SignificanceMetric = _ROOT.SignificanceMetric
AsimovSignificance = _ROOT.AsimovSignificance
PoissonSignificance = _ROOT.PoissonSignificance

# Classes from Event.h
Event = _ROOT.Event

# Classes from Tree.h
Tree = _ROOT.Tree

# Classes from CategoryReader.h
Category = _ROOT.Category
Categorizer = _ROOT.Categorizer
CategoryNode = _ROOT.CategoryNode
XMLCategorizer = _ROOT.XMLCategorizer

# Classes from dynamically compiled macros
loadEventsCSV = _ROOT.loadEventsCSV
loadEventsROOT = _ROOT.loadEventsROOT

# Factory function that returns an empty std::vector<Event*>
def Events():
    return _ROOT.vector('Event*')()

# Factory function to collect feature names into a std::vector<string>
def Features(*features):
    v = _ROOT.vector('string')()
    for feature in features:
        v.push_back(feature)
    return v

# Factory function to create an empty std::vector<string> which
# is filled with feature names ranked by net error reduction
def FeatureRanking():
    return _ROOT.vector('string')()

