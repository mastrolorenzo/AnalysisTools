# Author: Silvio Donato
#

from ROOT import *

//prefix = "/eos/uscms/store/group/lpchbb/HeppyNtuples/V24/"
prefix = ""

presel = "(Vtype==2 || Vtype==3)"

def getWeight(fileInc, fileB, region):
    region += "&&" + presel
    print region
    tree = TChain("tree")
    for fileIncName in fileInc.split(','):
        print tree.Add(prefix+fileIncName+"*.root")
        print prefix+fileIncName+"*.root"
    countInc    = 1.* tree.Draw("",region)
    print countInc
    tree.Reset()

    tree2 = TChain("tree")
    for fileBName in fileB.split(','):
        print tree2.Add(prefix+fileBName+"*.root")
        print prefix+fileBName+"*.root"
    countB      = 1.* tree2.Draw("",region)
    print countB
    tree2.Reset()

    if countInc == 0:
        weight = 0.
    else:
        weight = countInc/(countB+countInc)
    print weight
    return weight

WjetsHT100       = "/eos/uscms/store/group/lpchbb/HeppyNtuples/V25/WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/VHBB_HEPPY_V25_WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-Py8__RunIISummer16MAv2-PUMoriond17_80r2as_2016_TrancheIV_v6-v1/170128_122742/0000/,/eos/uscms/store/group/lpchbb/HeppyNtuples/V25/WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/VHBB_HEPPY_V25_WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-Py8__RunIISummer16MAv2-PUMoriond17_80r2as_2016_TrancheIV_v6_ext1-v1/170128_122829/0000/,/eos/uscms/store/group/lpchbb/HeppyNtuples/V25/WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/VHBB_HEPPY_V25_WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-Py8__RunIISummer16MAv2-PUMoriond17_80r2as_2016_TrancheIV_v6_ext2-v1/170128_125240/0000/"
WjetsHT200       = "/eos/uscms/store/group/lpchbb/HeppyNtuples/V25/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/VHBB_HEPPY_V25_WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-Py8__RunIISummer16MAv2-PUMoriond17_80r2as_2016_TrancheIV_v6-v1/170128_125327/0000/,/eos/uscms/store/group/lpchbb/HeppyNtuples/V25/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/VHBB_HEPPY_V25_WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-Py8__RunIISummer16MAv2-PUMoriond17_80r2as_2016_TrancheIV_v6_ext1-v1/170128_123055/0000/,/eos/uscms/store/group/lpchbb/HeppyNtuples/V25/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/VHBB_HEPPY_V25_WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-Py8__RunIISummer16MAv2-PUMoriond17_80r2as_2016_TrancheIV_v6_ext2-v1/170128_123142/0000"
WjetsHT400       = "/eos/uscms/store/group/lpchbb/HeppyNtuples/V25/WJetsToLNu_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/VHBB_HEPPY_V25_WJetsToLNu_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-Py8__RunIISummer16MAv2-PUMoriond17_80r2as_2016_TrancheIV_v6-v1/170128_123406/0000/,/eos/uscms/store/group/lpchbb/HeppyNtuples/V25/WJetsToLNu_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/VHBB_HEPPY_V25_WJetsToLNu_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-Py8__RunIISummer16MAv2-PUMoriond17_80r2as_2016_TrancheIV_v6_ext1-v1/170128_123458/0000/"
WjetsHT600       = "/eos/uscms/store/group/lpchbb/HeppyNtuples/V25/WJetsToLNu_HT-600To800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/VHBB_HEPPY_V25_WJetsToLNu_HT-600To800_TuneCUETP8M1_13TeV-madgraphMLM-Py8__RunIISummer16MAv2-PUMoriond17_80r2as_2016_TrancheIV_v6-v1/170128_123558/0000/,/eos/uscms/store/group/lpchbb/HeppyNtuples/V25/WJetsToLNu_HT-600To800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/VHBB_HEPPY_V25_WJetsToLNu_HT-600To800_TuneCUETP8M1_13TeV-madgraphMLM-Py8__RunIISummer16MAv2-PUMoriond17_80r2as_2016_TrancheIV_v6_ext1-v1/170206_171038/0000/"
WjetsHT800       = "/eos/uscms/store/group/lpchbb/HeppyNtuples/V25/WJetsToLNu_HT-800To1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/VHBB_HEPPY_V25_WJetsToLNu_HT-800To1200_TuneCUETP8M1_13TeV-madgraphMLM-Py8__RunIISummer16MAv2-PUMoriond17_80r2as_2016_TrancheIV_v6_ext1-v1/170128_123824/0000/"
WjetsHT1200      = "/eos/uscms/store/group/lpchbb/HeppyNtuples/V25/WJetsToLNu_HT-1200To2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/VHBB_HEPPY_V25_WJetsToLNu_HT-1200To2500_TuneCUETP8M1_13TeV-madgraphMLM-Py8__RunIISummer16MAv2-PUMoriond17_80r2as_2016_TrancheIV_v6-v1/170128_122917/0000/,/eos/uscms/store/group/lpchbb/HeppyNtuples/V25/WJetsToLNu_HT-1200To2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/VHBB_HEPPY_V25_WJetsToLNu_HT-1200To2500_TuneCUETP8M1_13TeV-madgraphMLM-Py8__RunIISummer16MAv2-PUMoriond17_80r2as_2016_TrancheIV_v6_ext1-v1/170128_123006/0000/"
WjetsHT2500      = "/eos/uscms/store/group/lpchbb/HeppyNtuples/V25/WJetsToLNu_HT-2500ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/VHBB_HEPPY_V25_WJetsToLNu_HT-2500ToInf_TuneCUETP8M1_13TeV-madgraphMLM-Py8__RunIISummer16MAv2-PUMoriond17_80r2as_2016_TrancheIV_v6-v1/170128_123231/0000/,/eos/uscms/store/group/lpchbb/HeppyNtuples/V25/WJetsToLNu_HT-2500ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/VHBB_HEPPY_V25_WJetsToLNu_HT-2500ToInf_TuneCUETP8M1_13TeV-madgraphMLM-Py8__RunIISummer16MAv2-PUMoriond17_80r2as_2016_TrancheIV_v6_ext1-v1/170128_123317/0000/"

ZjetsHT100       = "ZJetsToNuNu_HT-100To200_13TeV-madgraph"
ZjetsHT200       = "ZJetsToNuNu_HT-200To400_13TeV-madgraph"
ZjetsHT400       = "ZJetsToNuNu_HT-400To600_13TeV-madgraph"
ZjetsHT600       = "ZJetsToNuNu_HT-600ToInf_13TeV-madgraph"

ZLLjetsHT0     = "DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8"
ZLLjetsHT100     = "DYJetsToLL_M-50_HT-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8"
ZLLjetsHT200     = "DYJetsToLL_M-50_HT-200to400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8"
ZLLjetsHT400     = "DYJetsToLL_M-50_HT-400to600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8"
ZLLjetsHT600     = "DYJetsToLL_M-50_HT-600toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8"

ZLLBjets        = "DYBJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8"

ZBjets          = "DYBJetsToNuNu_Zpt-40toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8"

WBjets          = "/eos/uscms/store/group/lpchbb/HeppyNtuples/V25/WBJetsToLNu_Wpt-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/VHBB_HEPPY_V25_WBJetsToLNu_Wpt-100to200_TuneCUETP8M1_13TeV-madgraphMLM-Py8__RunIISummer16MAv2-PUMoriond17_80r2as_2016_TrancheIV_v6-v1/170130_121737/0000/"
WjetsBgen       = "/eos/uscms/store/group/lpchbb/HeppyNtuples/V25/WJetsToLNu_BGenFilter_Wpt-200toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/VHBB_HEPPY_V25_WJetsToLNu_BGenFilter_Wpt-200toInf_TuneCUETP8M1_13TeV-madgraphMLM-Py8__RunIISummer16MAv2-PUMoriond17_80r2as_2016_TrancheIV_v6-v1/170128_122607/0000/"

DYBJets             = "(lheNb>0 && lheV_pt>40)"
DYLightJets         = "!(lheNb>0 && lheV_pt>40)"

#WBJets             = " (lheNb>0  && lheV_pt>40)"
WBJets             = " (lheNb>0  && lheV_pt>100 && lheV_pt<200)"
#WJetsBGen           = " (lheNb==0 && nGenStatus2bHad>0  && lheV_pt>40)"
WJetsBGen           = " (lheNb==0 && nGenStatus2bHad>0  && lheV_pt>200)"
WLightJets          = "((lheNb==0 && nGenStatus2bHad==0 && lheV_pt>40) || (lheV_pt<40))"

HT0            = "(lheHT<100)"
HT100          = "(lheHT>100&&lheHT<200)"
HT200          = "(lheHT>200&&lheHT<400)"
HT400          = "(lheHT>400&&lheHT<600)"
HT600          = "(lheHT>600&&lheHT<800)"
HT800          = "(lheHT>800&&lheHT<1200)"
HT1200          = "(lheHT>1200&&lheHT<2500)"
HT2500          = "(lheHT>2500)"

print "weightWBjetsHT100=\t%.2f" %getWeight(WjetsHT100,   WBjets, HT100+"&&"+WBJets)
print "weightWBjetsHT200=\t%.2f" %getWeight(WjetsHT200,   WBjets, HT200+"&&"+WBJets)
print "weightWBjetsHT400=\t%.2f" %getWeight(WjetsHT400,   WBjets, HT400+"&&"+WBJets)
print "weightWBjetsHT600=\t%.2f" %getWeight(WjetsHT600,   WBjets, HT600+"&&"+WBJets)
print "weightWBjetsHT800=\t%.2f" %getWeight(WjetsHT800,   WBjets, HT800+"&&"+WBJets)
print "weightWBjetsHT1200=\t%.2f" %getWeight(WjetsHT1200,   WBjets, HT1200+"&&"+WBJets)
print "weightWBjetsHT25000=\t%.2f" %getWeight(WjetsHT2500,   WBjets, HT2500+"&&"+WBJets)
print ""
print "weightWjetsBgenHT100=\t%.2f" %getWeight(WjetsHT100,   WjetsBgen, HT100+"&&"+WJetsBGen)
print "weightWjetsBgenHT200=\t%.2f" %getWeight(WjetsHT200,   WjetsBgen, HT200+"&&"+WJetsBGen)
print "weightWjetsBgenHT400=\t%.2f" %getWeight(WjetsHT400,   WjetsBgen, HT400+"&&"+WJetsBGen)
print "weightWjetsBgenHT600=\t%.2f" %getWeight(WjetsHT600,   WjetsBgen, HT600+"&&"+WJetsBGen)
print "weightWjetsBgenHT800=\t%.2f" %getWeight(WjetsHT800,   WjetsBgen, HT800+"&&"+WJetsBGen)
print "weightWjetsBgenHT1200=\t%.2f" %getWeight(WjetsHT1200,   WjetsBgen, HT1200+"&&"+WJetsBGen)
print "weightWjetsBgenHT2500=\t%.2f" %getWeight(WjetsHT2500,   WjetsBgen, HT2500+"&&"+WJetsBGen)
print ""
#print "weightZBjetsHT100=\t%.2f" %getWeight(ZjetsHT100,   ZBjets, HT100+"&&"+DYBJets)
#print "weightZBjetsHT200=\t%.2f" %getWeight(ZjetsHT200,   ZBjets, HT200+"&&"+DYBJets)
#print "weightZBjetsHT400=\t%.2f" %getWeight(ZjetsHT400,   ZBjets, HT400+"&&"+DYBJets)
#print "weightZBjetsHT600=\t%.2f" %getWeight(ZjetsHT600,   ZBjets, HT600+"&&"+DYBJets)
#print ""
#print "weightZLLBjetsHT0=\t%.2f" %getWeight(ZLLjetsHT0,   ZLLBjets, HT0+"&&"+DYBJets)
#print "weightZLLBjetsHT100=\t%.2f" %getWeight(ZLLjetsHT100,   ZLLBjets, HT100+"&&"+DYBJets)
#print "weightZLLBjetsHT200=\t%.2f" %getWeight(ZLLjetsHT200,   ZLLBjets, HT200+"&&"+DYBJets)
#print "weightZLLBjetsHT400=\t%.2f" %getWeight(ZLLjetsHT400,   ZLLBjets, HT400+"&&"+DYBJets)
#print "weightZLLBjetsHT600=\t%.2f" %getWeight(ZLLjetsHT600,   ZLLBjets, HT600+"&&"+DYBJets)


