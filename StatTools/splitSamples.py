import ROOT
import sys
from math import sqrt
from math import pow
import numpy
import argparse
import re

## Create histogram with histograms of a given BDT shape, one histogram for each sample. Used in creation of datacards.
##
## Author: Stephane Cooperstein
##

parser = argparse.ArgumentParser("Create histograms for datacards")
parser.add_argument('-i', '--inputfile', type=str, default="", help="The input root file (ntuple)")
parser.add_argument('-c', '--catName', type=str, default="Cat", help="The category label for datacard")
parser.add_argument('-p', '--preselection', type=str, default="", help="The category selection cuts on top of ntuple")
parser.add_argument('-s', '--systematics', type=str, default="", help="The systematics config file")
parser.add_argument('-w', '--weights', type=str, default="", help="Comma-separated list of weights to apply on top of nominal 'weight'")
parser.add_argument('-d', '--doData', type=int, default=0, help="If non-zero run on real data, if 0 compute data as sum of background pdf's (default 0)")
parser.add_argument('-t', '--dataTree', type=str, default="", help="If doData is true, you specify here the ntuple with the real data events")
parser.add_argument('-tt', '--ttbarTree', type=str, default="", help="Specify separate file with ttbar powheg so it only takes forever when running tt")
parser.add_argument('-wj', '--wjetsTree', type=str, default="", help="Specify separate file with W+jets MC")
parser.add_argument('-vz', '--vzTree', type=str, default="", help="Specify separate file with VZ MC")
parser.add_argument('-v', '--varname', type=str, default="CMS_vhbb_BDT_Wln_13TeV", help="The name of the variable shape which goes into the histograms and will be fitted")
parser.add_argument('-xl','--xlow', type=float, default=-1.0, help="Lowest bin edge of fitted distribution")
parser.add_argument('-xh','--xhigh', type=float, default=1.0, help="Highest bin edge of fitted distribution")
parser.add_argument('-bh','--binhigh',type=float,default=1.0, help="manually set low bin edge of most sensitive bin")
parser.add_argument('-bl','--binlow',type=float,default=-1.0, help="manually set low bin edge of least sensitive bin")
parser.add_argument('-o', '--ofilename',type=str,default="", help="output histogram file name (default hists_[cat].root)")
parser.add_argument('-b', '--binstats', type=str, default="", help="Text file listing all the individual bin. stat. uncertainties to include")
parser.add_argument('-n', '--nbins', type=int, default=20, help="number of bins in datacard shapes")
parser.add_argument('-r', '--doRebin',type=int, default=1, help="If not 0 then rebin histogram so that bins are sufficiently populated in background (default=1)")
parser.add_argument('-tol','--tolerance',type=float,default=0.5,help="Tolerance threshold on dB/sqrt(B) for inclusion of bin-by-bin stat. shape uncertainties (default=0.50)")
parser.add_argument('-bb','--binBoundaries',type=str,default="",help="Manually specify the bin edges for the histogram that goes into the histogram")
parser.add_argument('-vv','--doVV', type=bool, default=False, help="If true do VV analysis (default False)")
parser.add_argument('-donlo','--doNLOWJets', type=bool, default=False, help="If true use NLO W+jets samples, if false use LO (default False)")
parser.add_argument('-rwtt','--reweightTT', type=bool, default=False, help="If true assume full TT statistics, if false use only half (default False)")
parser.add_argument('-si','--oldSI', type=bool, default=False, help="If true use old sampleIndex def., if false use new (correct) definition (default False)")
parser.add_argument('-sa', '--sample', type=str, default="", help="Specify to run on only one particular sample")
parser.add_argument('--drawFromNom', type=bool, default=False, help="If true, do not try to draw args.varname_systname but simply use nominal args.varname")
args = parser.parse_args()
print args

#if (len(sys.argv) != 4 and len(sys.argv) != 5):
#    print "Takes three or four arguments:"
#    print "python splitSamples.py [input_root_file] [catName] [presel]"
#    print "python splitSamples.py [input_root_file] [catName] [presel] [sys.txt] (for shape systematics besides bin-by-bin)"
#    sys.exit(1)

ROOT.gROOT.SetBatch(True)

#print "going to open this file now: ",args.inputfile
ipath = args.inputfile
#ifile = ROOT.TFile.Open(args.inputfile, "r")
catName = args.catName
presel = args.preselection

if args.ofilename == "":
    ofile = ROOT.TFile.Open("hists_%s.root" % catName, "RECREATE")
else:
    ofile = ROOT.TFile.Open(args.ofilename, "RECREATE")

systematics = {}
if (args.systematics != ""):
    # prepare shape systematics
    sysfile = open(args.systematics)
    for line in sysfile:
        if (line[0] == '#'): continue
        line = line.strip()
        params = line.split(' ')
        paramsToKeep = []
        for param in params:
            if (param != ''):
                paramsToKeep.append(param)
        print paramsToKeep
        if (paramsToKeep[1] != "shape"): continue
        sysSamples = paramsToKeep[3].split(',')
        # map systematic name to the corresponding weight in the ntuple and possibly an extension if it has a different BDT/Mjj output than the nominal
        # In that case, the BDT name used is [bdtname]_[paramsToKeep[5]]UpDown
        if (len(paramsToKeep) > 5):
            systematics[paramsToKeep[0]] = (paramsToKeep[4], paramsToKeep[5], sysSamples)
        else:
            systematics[paramsToKeep[0]] = (paramsToKeep[4], "", sysSamples)
    
print systematics

set_of_weights = []
if (args.weights != ""):
    set_of_weights = args.weights.split(',')
print set_of_weights
#weight_string = "weight"
#weight_string = "weight*(bTagWeightEF/bTagWeight)*(1+isWenu*(-1+(SF_HLT_Ele23_WPLoose[lepInd]/EffHLT_Ele27_WPLoose_Eta2p1[lepInd])))"
#weight_string = "(1./CS_SF)*weight*(1 + (sampleIndex==120)*(-1 + 0.5))"
if args.reweightTT:
    weight_string = "weight*(1 + (sampleIndex==120)*(-1 + 0.5))"
else:
    # HACK because I screwed up the norm on the NLO Z+jets samples...
    #weight_string = "weight*(1+(sampleIndex==200)*(-1+(0.0002014*abs(genWeight))))*(1+(sampleIndex==201)*(-1+(0.01679*abs(genWeight))))*(1+(sampleIndex==202)*(-1+(0.1387*abs(genWeight))))*(1+(sampleIndex==203)*(-1+(1.148*abs(genWeight))))*(1+(sampleIndex==204)*(-1+(0.0003781*abs(genWeight))))*(1+(sampleIndex==205)*(-1+(0.01439*abs(genWeight))))*(1+(sampleIndex==206)*(-1+(0.119579*abs(genWeight))))*(1+(sampleIndex==207)*(-1+(1.007*abs(genWeight))))"
    weight_string = "weight"
#weight_string = "(1./CS_SF)*weight*(1 + (sampleIndex==17)*(-1 + 3.07))"
#weight_string = "(1./CS_SF)*(bTagWeightICHEP/bTagWeight)*weight"
#weight_string = "(1./CS_SF)*weight*(1+isWenu*(-1+(SF_HLT_Ele23_WPLoose[lepInd]/EffHLT_Ele27_WPLoose_Eta2p1[lepInd])))"
#weight_string = "weight*(1+isWenu*(-1+(1./SF_egammaEffi_tracker[lepInd])))*(1+isWenu*(-1+(SF_HLT_Ele23_WPLoose[lepInd]/EffHLT_Ele27_WPLoose_Eta2p1[lepInd])))"
#weight_string = "weight*(1./SF_egammaEffi_tracker[lepInd])*(1+isWenu*(-1+(SF_HLT_Ele23_WPLoose[lepInd]/EffHLT_Ele27_WPLoose_Eta2p1[lepInd])))"
#weight_string = "weight*(puWeight/weight_PU)*(1./CS_SF)"
#weight_string = "weight*(1./CS_SF)*(puWeight/weight_PU)*(1 + (sampleIndex>=24&&sampleIndex<=31)*(-0.63)) * ( 1 + (sampleIndex>=24&&sampleIndex<=31&&isWmunu)*(-1 + 1./(0.0673 *SF_SMuTrig_Block1[lepInd] + 0.9327 * SF_SMuTrig_Block2[lepInd]))  ) * ( 1 + (sampleIndex>=24&&sampleIndex<=31&&isWenu)*(-1 + 1./(SF_HLT_Ele23_WPLoose[lepInd]))  )"
#weight_string = "weight*(1./CS_SF)*(puWeight/weight_PU)*(1 + (sampleIndex>=24&&sampleIndex<=31)*(-0.63)) * ( 1 + (sampleIndex>=24&&sampleIndex<=31&&isWmunu)*(-1 + 1./(0.0673 *SF_SMuTrig_Block1[lepInd] + 0.9327 * SF_SMuTrig_Block2[lepInd]))  ) * ( 1 + (sampleIndex>=24&&sampleIndex<=31&&isWenu)*(-1 + 1./(SF_HLT_Ele23_WPLoose[lepInd]))  )"
#weight_string = "weight*(CS_SF_new/CS_SF)"
#weight_string = "weight*((sampleIndex!=50&&sampleIndex!=51&&sampleIndex!=52) + (sampleIndex==50||sampleIndex==51||sampleIndex==52)*(245.79/831.76))" ## FIXME: change it back!!!!!
#weight_string = "weight*((sampleIndex!=50&&sampleIndex!=51&&sampleIndex!=52) + (sampleIndex==50||sampleIndex==51||sampleIndex==52)*(245.79/831.76)*1.956)" ## FIXME: change it back!!!!!
for weight in set_of_weights:
    weight_string += "*" + weight
print weight_string

tree_mc = ROOT.TChain("Events")
#tree_mc = ifile.Get("Events")
#bdtname = "BDT_wMass_Dec14_3000_5"
#bdtname = "BDT_wMass_Dec4"
#bdtname = "CMS_vhbb_BDT_Wln_13TeV"
bdtname = args.varname # not necessarily the bdt shape, can fit whatever shape is specified
#nBins = 1000
nBins = args.nbins  # number of bins in final histogram after rebinning
nBinsFine = 1000 # number of candidate bin edges to begin with
#tolerance = 1.5 # dB/B tolerance for bin-by-bin stat. uncertainties
#tolerance = 0.5 # dB/B tolerance for bin-by-bin stat. uncertainties
#tolerance = 0.75 # dB/B tolerance for bin-by-bin stat. uncertainties
tolerance = args.tolerance

sampleMap = {} # map sampleNames to list of sampleIndex's
sampleMapAltModel = {} # alternate MC samples for model shape systematics
sampleNameMap = {}

from nano_samples import the_samples_dict
#from nano_samples_znn2016 import the_samples_dict
for sample in the_samples_dict:
    sampleMap[the_samples_dict[sample][2]] = []
    sampleNameMap[the_samples_dict[sample][2]] = []
for sample in the_samples_dict:
    sampleMap[the_samples_dict[sample][2]].append(the_samples_dict[sample][0])
    sampleNameMap[the_samples_dict[sample][2]].extend(the_samples_dict[sample][3])
print "sampleMap = "
print sampleMap
print "sampleNameMap = "
print sampleNameMap


allSampInd = [] # list of all indices for all backgrounds
sigSamps = ["WH_hbb","ZH_hbb"]
if args.doVV:
    sigSamps = ["VVHF","VVLF"]
print "sigSamps = "
print sigSamps
for sample in sampleMap:
    if (sample in sigSamps or sample == "data_obs"): continue
    allSampInd.extend(sampleMap[sample])
sampleMap["Bkg"] = allSampInd
print "all Bkg sampleInd = "
print allSampInd
if not args.doData:
    sampleMap["data_obs"] = allSampInd # fake data as sum of all background MC

def makeCutString(sample, sMap):
    cutString = presel + "&&("
    for index in sMap[sample]:
            if isinstance(index,str):
                cutString += "(%s)||" % index
            else:
                cutString += "(sampleIndex==%i)||" % index
    cutString = cutString[0:len(cutString)-2]
    cutString += ")"
    return cutString

print "Bkg string = ",makeCutString("Bkg", sampleMap)
#print "Wj0b string = ",makeCutString("Wj0b", sampleMap)
print "Thank you for your attention!"

binBoundaries = numpy.zeros(nBins+1,dtype=float)
binBoundaries[0] = args.xlow
binBoundaries[nBins] = args.xhigh
hBkg = ROOT.TH1F("hBkg","hBkg",nBinsFine,args.xlow,args.xhigh)
hSig = ROOT.TH1F("hSig","hSig",nBinsFine,args.xlow,args.xhigh)
if args.doRebin:
    # first rebin the histogram so that the first and last bins are not empty and have less than 35% stat. uncertainty
    #tree.Draw("%s>>hBkg" % bdtname,"((%s)&&sampleIndex>0)*weight*(2.2/1.28)" % presel)
    bkgCutString = makeCutString("Bkg", sampleMap)
    print bkgCutString
    if args.doVV:
        sigCutString = makeCutString("VVHF", sampleMap)
    else:
        sigCutString = makeCutString("WH_hbb", sampleMap)

    fnames = [] # avoid adding the same tree twice
    tree_sig = ROOT.TChain("Events")
    tree_bkg = ROOT.TChain("Events")
    for sample in sampleNameMap:
        if sample in sigSamps:
            for sname in sampleNameMap[sample]:
                #fname = ipath + "/sum_" + sname + "_3.root"
                #fname = ipath + "/sum_" + sname + "_weighted2.root"
                fname = ipath + "/sum_" + sname + ".root"
                tree_sig.Add(fname)
        elif sample != "data_obs":
            for sname in sampleNameMap[sample]:
                #fname = ipath + "/sum_" + sname + "_3.root"
                #fname = ipath + "/sum_" + sname + "_weighted2.root"
                fname = ipath + "/sum_" + sname + ".root"
                if fname not in fnames:
                    print "adding  %s for background" % fname
                    tree_bkg.Add(fname)
                    fnames.append(fname)
         
    tree_bkg.Draw("%s>>hBkg" % bdtname,"((%s)&&(%s))*%s" % (presel,bkgCutString,weight_string))
    tree_sig.Draw("%s>>hSig" % bdtname,"((%s)&&(%s))*%s" % (presel,sigCutString,weight_string))

    print "total signal, total background:",hSig.Integral(),", ",hBkg.Integral()

    foundLowBinEdge = False
    foundHighBinEdge = False
    if (args.binlow != -1):
        foundLowBinEdge = True
        binBoundaries[1] = args.binlow
    if (args.binhigh != 1):
        foundHighBinEdge = True
        binBoundaries[nBins-1] = args.binhigh
    B_err2_low = 0.
    B_err2_high = 0.
    for ibin in range(1,nBinsFine):
        if not foundLowBinEdge:
            B_low = hBkg.Integral(1,ibin)
            B_err2_low += pow(hBkg.GetBinError(ibin),2) 
            #S_low = hSig.Integral(1,ibin)
            #print "edge = ",hBkg.GetBinLowEdge(ibin+1)
            #print "B_low = ",B_low
            #if (B_low > 0):
            #    print "1./sqrt(B_low) = ",1./sqrt(B_low)
            #print B_low,ibin,hBkg.Integral(),args.xlow,args.xhigh
            #print "B_low,sqrt(B_err2_low)",B_low,sqrt(B_err2_low)
            #if (B_low > 0 and 1./sqrt(B_low) < 0.35):
            if (B_low > 0 and sqrt(B_err2_low)/B_low < 0.35):
                binBoundaries[1] = hBkg.GetBinLowEdge(ibin+1)
                foundLowBinEdge = True
                #print "found low bin!",ibin
        if not foundHighBinEdge:
            B_high = hBkg.Integral(nBinsFine-ibin+1, nBinsFine) 
            B_err2_high += pow(hBkg.GetBinError(nBinsFine-ibin+1),2)
            S_high = hSig.Integral(nBinsFine-ibin+1, nBinsFine) 
            #S_high = hSig.Integral(nBinsFine-ibin+1, nBinsFine) 
            #print "B_high,sqrt(B_err2_high)",B_high,sqrt(B_err2_high)
            #if (B_high > 0 and 1./sqrt(B_high) < 0.35):
            #if (B_high > 0 and sqrt(B_err2_high)/B_high < 0.35 and S_high >= 1.0 and abs(hBkg.GetBinLowEdge(nBinsFine-ibin+1)-0.608)<0.001):
            #if (B_high > 0 and sqrt(B_err2_high)/B_high < 0.35 and S_high >= 1.0):
            #if (B_high > 0 and sqrt(B_err2_high)/B_high < 0.35 and S_high >= 6.0):
            if (B_high > 0 and sqrt(B_err2_high)/B_high < 0.35):
                binBoundaries[nBins-1] = hBkg.GetBinLowEdge(nBinsFine-ibin+1)
                foundHighBinEdge = True
                #print "found high bin!",(nBinsFine-ibin+1)
    # split the middle bins equidistantly
    for i in range(2,nBins-1):
        binBoundaries[i] = binBoundaries[1] + (i-1)*((binBoundaries[nBins-1] - binBoundaries[1])/(nBins-2)) 

if not args.doRebin:
    for i in range(1,nBins):
        binBoundaries[i] = binBoundaries[0] + (i)*((binBoundaries[nBins] - binBoundaries[0])/(nBins))
if (args.binBoundaries != ""):
    boundaryString = args.binBoundaries
    print "parsing boundaries from ",boundaryString
    boundaries = boundaryString.split(',')
    nBins = len(boundaries)-1
    binBoundaries = numpy.zeros(nBins+1,dtype=float)
    for i in range(len(boundaries)):
        binBoundaries[i] = float(boundaries[i].replace('m','-'))
    print binBoundaries
    
hBkg = hBkg.Rebin(nBins, "", binBoundaries)
print binBoundaries
if args.binstats == "":
    otextfile = open("binStats_%s.txt" % catName, "w")
else:
    otextfile = open(args.binstats,"w")
tree = ROOT.TChain("Events")
#tree = ROOT.TTree("Events","Events")
#hBkg.Write()
for sample in sampleMap:
    if (args.sample != "" and sample != args.sample): continue
    #if (sample!="s_Top" and sample!="VVHF" and sample!="VVLF"): continue
    #if (sample != "Wj0b" and sample!="Wj1b" and sample!="Wj2b"): continue
    #if (sample != "Bkg" and sample != "data_obs" and sample != "WH_hbb" and sample != "ZH_hbb"): continue
    #if (sample != "VVHF"): continue
    #if (sample != "Bkg"): continue
    #if (args.ttbarTree != "" and (sample == "TT" or sample == "Bkg")): 
    #    ifile_tt = ROOT.TFile.Open(args.ttbarTree,"r")
    #    tree_tt = ifile_tt.Get("Events")
    #    if (sample == "TT"):
    #        tree = tree_tt
    #    else:
    #        chain = ROOT.TChain("Events")
    #        chain.Add(args.inputfile)
    #        chain.Add(args.ttbarTree)
    #        tree = chain
    #    #ofile.cd()
    #    #ifile_tt.Close()
    #elif (sample == "WH_hbb" or sample == "ZH_hbb"):
    #    ifile_sig = ROOT.TFile.Open(args.inputfile.replace("output_mc","output_signal"))
    #    tree_sig = ifile_sig.Get("Events")
    #    tree = tree_sig
    #elif (sample == "VVHF" or sample == "VVLF"):
    #    ifile_vv = ROOT.TFile.Open(args.inputfile.replace("output_mc","output_vv"))
    #    tree_vv = ifile_vv.Get("Events")
    #    tree = tree_vv
    #elif (args.wjetsTree != "" and (sample == "Wj0b" or sample == "Wj1b" or sample=="Wj2b")):
    #    ifile_wjets = ROOT.TFile.Open(args.wjetsTree)
    #    tree_wjets = ifile_wjets.Get("Events")
    #    tree = tree_wjets
    #elif (args.vzTree != "" and (sample == "VVHF" or sample == "VVLF")):
    #    ifile_vz = ROOT.TFile.Open(args.vzTree)
    #    tree_vz = ifile_vz.Get("Events")
    #    tree = tree_vz
    #else:
    #    tree = tree_mc 
    ##if sample != "Bkg" and sample != "data_obs": 
    if sample != "Bkg": 
        for sname in sampleNameMap[sample]:
            fname = ipath + "/sum_" + sname + ".root"
            #fname = ipath + "/sum_" + sname + "_weighted2.root"
            #fname = ipath + "/sum_" + sname + "_3.root"
            #print tree.GetEntries()
            tree.Add(fname)
            print "Added for sample %s: %s" % (sample,fname)
            #print sname,tree.GetEntries()
    elif sample == "Bkg":
        for sample in sampleNameMap:
            if sample in sigSamps: continue
            for sname in sampleNameMap[sample]: 
                fname = ipath + "/sum_" + sname + ".root"
                #fname = ipath + "/sum_" + sname + "_weighted2.root"
                #fname = ipath + "/sum_" + sname + "_3.root"
                tree.Add(fname)
     
    #cutString = presel + "&&("
    #for index in sampleMap[sample]:
    #        cutString += "(sampleIndex==%i)||" % index
    #cutString = cutString[0:len(cutString)-2]
    #cutString += ")"
    #print cutString
    cutString = makeCutString(sample, sampleMap)
    print cutString
    #hBDT = ROOT.TH1F(sample,sample,nBins,-1,1)
    hBDT = ROOT.TH1F(sample,sample,nBins,binBoundaries)
    #tree.Draw("%s>>%s" % (bdtname, sample),"(%s)*weight" % cutString)
    #tree.Draw("%s>>%s" % (bdtname, sample),"(%s)*weight*(2.2/1.28)" % cutString)  
    if (sample == "data_obs" and args.doData):
        #ifile_data = ROOT.TFile.Open(args.dataTree,"r")
        #tree_data = ifile_data.Get("Events")
        ofile.cd()
        # make sure we don't weight actual data by puWeight, SF's, etc.
        #tree_data.Draw("%s>>%s" % (bdtname, sample),"((%s)&&Pass_nominal)*sb_weight3" % (cutString)) ## FIXME
        tree.Draw("%s>>%s" % (bdtname, sample),"((%s)&&Pass_nominal)" % (cutString)) 
        #tree_data.Draw("BDT_V25_March27_400_5>>%s" % (sample),"((%s)&&Pass_nominal)" % (cutString)) # HACK FIXME PLEASE!!
        #ifile_data.Close()
    elif (sample == "data_obs"):
        # fake data which is sum of all MC
        #hBkg.Write(sample)
        ofile.cd()
        hBkg.Write("BDT_%s_%s" % (catName,sample))
         
    #    continue
    #elif (sample == "TT" and args.ttbarTree != ""):
    #    print "got here"
    #    ifile_tt = ROOT.TFile.Open(args.ttbarTree,"r")
    #    tree_tt = ifile_tt.Get("Events")
    #    ofile.cd()
    #    # make sure we don't weight actual data by puWeight, SF's, etc.
    #    tree_tt.Draw("%s>>%s" % (bdtname, sample),"((%s)&&Pass_nominal)*%s" % (cutString,weight_string)) 
    #    ifile_tt.Close() 
    else:
        tree.Draw("%s>>%s" % (bdtname, sample),"((%s)&&Pass_nominal)*%s" % (cutString,weight_string)) 
        #print "just drew tree for sample %s, now chain.Print()" % sample
        #tree.Print()
    print "tree.Draw(\"%s>>%s\",\"((%s)&&Pass_nominal)*%s\"" % (bdtname,sample,cutString,weight_string) 
    #hBDT = hBDT.Rebin(nBins, "", binBoundaries)
    # Add bin-by-bin stat. uncertainties
    if (sample not in ["data_obs"] and (tolerance <= 0.5 or sample not in ["QCD","VVLF","VVHF","WH_hbb","ZH_hbb"]) ): # assuming for SR the tolerance does not go above 0.50
    #if (sample not in ["QCD","VVLF","VVHF","WH_hbb","ZH_hbb","data_obs"]): # can exclude these when running on background to reduce the number of nuisances
        for ibin in range(1, hBDT.GetNbinsX()+1):
            B = hBDT.GetBinContent(ibin)
            B_err = hBDT.GetBinError(ibin)
            #print "ibin = %i" % ibin
            #print "B = %f" % B
            #print "B_err = %f" % B_err
            #print "GetEntries = %i" % hBDT.GetEntries()
            #print "GetSumOfWeights = %f" % hBDT.GetSumOfWeights()
            #if (B > 0):
                #print "B_err/sqrt(B) = %f" % (B_err/sqrt(B))
                #print "B_err/B = %f" % (B_err/B)
            #B_eff = hBDT.GetEffectiveEntries() # effective statistical number of entries in bin considering weights

            if ( B > 0 and ( ( B >=1 and (B_err/sqrt(B)) > tolerance) or (B < 1 and B_err/B > tolerance) ) ):
                
                ## FIXME: just for dB/sqrt(B) sensitivity study
                #scale = 0.1
                #if (B>=1):
                #    B_err = scale * sqrt(B)
                #else:
                #    B_err = scale * B

                print "qualified as a bin stat. shape uncertainty! Bin %i, sample %s" % (ibin,sample)
                if (B >=1): print B_err/sqrt(B)
                else: print B_err/B
                print "Bin content: %f +/- %f" % (B,B_err)
                #hBinStat = ROOT.TH1F("CMS_vhbb_stat%s_%s_bin%i_13TeV" % (sample,catName,ibin),"CMS_vhbb_stat%s_%s_bin%i_13TeV" % (sample,catName,ibin),nBins,-1,1)
                hBinStatUp = hBDT.Clone()
                hBinStatDown = hBDT.Clone()
                hBinStatUp.SetName("BDT_%s_%s_CMS_vhbb_stat%s_%s_bin%i_13TeVUp" % (catName,sample,sample,catName,ibin))
                hBinStatDown.SetName("BDT_%s_%s_CMS_vhbb_stat%s_%s_bin%i_13TeVDown" % (catName,sample,sample,catName,ibin))
                #hBinStatUp.SetBinContent(ibin, B + sqrt(B))
                #hBinStatDown.SetBinContent(ibin, max(B - sqrt(B),0.000001))
                hBinStatUp.SetBinContent(ibin, B + B_err)
                hBinStatDown.SetBinContent(ibin, max(B - B_err,0.000001))
                otextfile.write("CMS_vhbb_stat%s_%s_bin%i_13TeV\n" % (sample,catName,ibin))
                ofile.cd()
                hBinStatUp.Write()
                hBinStatDown.Write()
    # other shape systematics
    for syst in systematics:
        sysWeight, sysName, sysSamples = systematics[syst]
        if sample not in sysSamples: continue
        #if (syst.find("LHE")==-1): continue
        ofile.cd()
        hBDTSystUp = ROOT.TH1F("BDT_%s_%s_%sUp" % (catName,sample,syst), "%s_%sUp" % (sample,syst),nBins,binBoundaries)
        hBDTSystDown = ROOT.TH1F("BDT_%s_%s_%sDown" % (catName,sample,syst), "%s_%sDown" % (sample,syst),nBins,binBoundaries)
        sysBDTNameUp = bdtname
        sysBDTNameDown = bdtname
        passSysup = "Pass_nominal"
        passSysdown = "Pass_nominal"
        cutStringdown = cutString
        cutStringup = cutString
        if (sysName != ""):
            if not args.drawFromNom:
                if (bdtname.find('[') != -1):
                    # drawing an array variable
                    sysBDTNameUp = bdtname[:bdtname.find('[')] + "_" + sysName + "Up" + bdtname[bdtname.find('['):] 
                else:
                    sysBDTNameUp += "_%sUp" % sysName
                    sysBDTNameDown += "_%sDown" % sysName

            passSysup = "Pass_%sUp" % sysName
            passSysdown = "Pass_%sDown" % sysName
            cutStringup = cutString.replace("controlSample","controlSample_%sUp"%sysName)
            cutStringdown = cutString.replace("controlSample","controlSample_%sDown"%sysName)
        if (sysWeight != "1.0" and sysWeight.find(".root") == -1):
            if (sysWeight.find(',') == -1):
                if (sysWeight.find("bTagWeight") != -1):
                    bTagWeightNom = sysWeight[:sysWeight.find('_')]
                    print "bTagWeightNom = "+bTagWeightNom
                    tree.Draw("%s>>BDT_%s_%s_%sUp" % (sysBDTNameUp, catName,sample, syst),"((%s)&&%s)*%s*(%sUp)" % (cutStringup,passSysup,weight_string,sysWeight))   
                    tree.Draw("%s>>BDT_%s_%s_%sDown" % (sysBDTNameDown, catName, sample, syst),"((%s)&&%s)*%s*(%sDown)" % (cutStringdown,passSysdown,weight_string,sysWeight))  
                elif (sysWeight.find("VPtCorrFactorSplit") != -1):
                     tree.Draw("%s>>BDT_%s_%s_%sUp" % (sysBDTNameUp, catName,sample, syst),"((%s)&&%s)*(1./VPtCorrFactorSplit3)*%s*(%sUp)" % (cutStringup,passSysup,weight_string,sysWeight)) 
                     tree.Draw("%s>>BDT_%s_%s_%sDown" % (sysBDTNameDown, catName, sample, syst),"((%s)&&%s)*(1./VPtCorrFactorSplit3)*%s*(%sDown)" % (cutStringdown,passSysdown,weight_string,sysWeight))
                     #if hBDTSystUp.Integral() > 0:
                     #    hBDTSystUp.Scale(hBDT.Integral()/hBDTSystUp.Integral())                
                     #if hBDTSystDown.Integral() > 0: 
                     #    hBDTSystDown.Scale(hBDT.Integral()/hBDTSystDown.Integral())                 
                #elif (sysWeight.find("weight_PU") != -1):
                #    tree.Draw("%s>>BDT_%s_%s_%sUp" % (sysBDTNameUp, catName,sample, syst),"((%s)&&%s)*(1./weight_PU)*%s*(%sUp)" % (cutString,passSys,weight_string,sysWeight))   
                #    tree.Draw("%s>>BDT_%s_%s_%sDown" % (sysBDTNameDown, catName, sample, syst),"((%s)&&%s)*(1./weight_PU)*%s*(%sDown)" % (cutString,passSys,weight_string,sysWeight))   
                #elif (sysWeight.find("weight_PU") != -1): ## FIXMEEE
                #     if (sample != "Wj0b" and sample != "Wj1b" and sample != "Wj2b"):
                #         tree.Draw("%s>>BDT_%s_%s_%sUp" % (sysBDTNameUp, catName,sample, syst),"((%s)&&%s)*(1./weight_PU)*%s*(%sUp)" % (cutString,passSys,weight_string,sysWeight))   
                #         tree.Draw("%s>>BDT_%s_%s_%sDown" % (sysBDTNameDown, catName, sample, syst),"((%s)&&%s)*(1./weight_PU)*%s*(%sDown)" % (cutString,passSys,weight_string,sysWeight))
                #     else:
                #         tree.Draw("%s>>BDT_%s_%s_%sUp" % (sysBDTNameUp, catName,sample, syst),"((%s)&&%s)*%s*(%sUp)" % (cutString,passSys,weight_string,sysWeight))   
                #         tree.Draw("%s>>BDT_%s_%s_%sDown" % (sysBDTNameDown, catName, sample, syst),"((%s)&&%s)*%s*(%sDown)" % (cutString,passSys,weight_string,sysWeight))
                #elif (sysWeight.find("puWeight") != -1):
                #    tree.Draw("%s>>BDT_%s_%s_%sUp" % (sysBDTNameUp, catName,sample, syst),"((%s)&&%s)*(1./puWeight)*%s*(%sUp)" % (cutString,passSys,weight_string,sysWeight))   
                #    tree.Draw("%s>>BDT_%s_%s_%sDown" % (sysBDTNameDown, catName, sample, syst),"((%s)&&%s)*(1./puWeight)*%s*(%sDown)" % (cutString,passSys,weight_string,sysWeight))   
                else:
                    tree.Draw("%s>>BDT_%s_%s_%sUp" % (sysBDTNameUp, catName,sample, syst),"((%s)&&%s)*%s*(%sUp)" % (cutStringup,passSysup,weight_string,sysWeight))   
                    tree.Draw("%s>>BDT_%s_%s_%sDown" % (sysBDTNameDown, catName, sample, syst),"((%s)&&%s)*%s*(%sDown)" % (cutStringdown,passSysdown,weight_string,sysWeight))  
            else:
                # separate branches for up/down variation
                sysWeightUp = sysWeight.split(',')[0]
                sysWeightDown = sysWeight.split(',')[1]

                #if (syst.find("LHE_weights_pdf")!=-1):
                #    if (sample != "TT"):
                #        sysWeightUp += "*LHE_weights_pdf_normwgt[%s]" % re.findall('\d+', sysWeightUp)[0] 
                #        sysWeightDown += "*LHE_weights_pdf_normwgt[%s]" % re.findall('\d+', sysWeightDown)[0]
                #    else:
                #        # have to do something special here since we didn't hadd everything on the ttbar jobs so that it doesn't take 100 million years to run the jobs
                #        ifile_ttpowhegcounts = ROOT.TFile.Open(args.inputfile.replace("output_mc","ttpowheg_counts"))
                #        CountWeightedLHEWeightPdf_TT_powheg = ifile_ttpowhegcounts.Get("CountWeightedLHEWeightPdf_TT_powheg")
                #        normweightUp = CountWeightedLHEWeightPdf_TT_powheg.GetBinContent(CountWeightedLHEWeightPdf_TT_powheg.FindBin(int(re.findall('\d+', sysWeightUp)[0])))
                #        normweightDown = CountWeightedLHEWeightPdf_TT_powheg.GetBinContent(CountWeightedLHEWeightPdf_TT_powheg.FindBin(int(re.findall('\d+', sysWeightDown)[0])))
                #        sysWeightUp += "*(nProcEvents/%f)" % normweightUp
                #        sysWeightDown += "*(nProcEvents/%f)" % normweightDown
                #        ifile_ttpowhegcounts.Close()
                #        ofile.cd()
                #print sysWeightUp,sysWeightDown 
                      
                #print "tree.Draw(\"%s>>BDT_%s_%s_%sUp\",\"((%s)&&%s)*%s*(%s)\")" % (sysBDTNameUp, catName,sample, syst,cutString,passSys,weight_string,sysWeightUp)
                tree.Draw("%s>>BDT_%s_%s_%sUp" % (sysBDTNameUp, catName,sample, syst),"((%s)&&%s)*%s*(%s)" % (cutStringup,passSysup,weight_string,sysWeightUp)) 
                #print hBDTSystUp.Integral() 
                #print "tree.Draw(\"%s>>BDT_%s_%s_%sUp\",\"((%s)&&%s)*%s*(%s)\")" % (sysBDTNameDown, catName,sample, syst,cutString,passSys,weight_string,sysWeightDown)
                tree.Draw("%s>>BDT_%s_%s_%sDown" % (sysBDTNameDown, catName, sample, syst),"((%s)&&%s)*%s*(%s)" % (cutStringdown,passSysdown,weight_string,sysWeightDown))  
                #print hBDTSystDown.Integral()
                
                if (syst.find("LHE_weights_pdf")!=-1): 
                    # new prescription for LHE PDF nuisances, treat each variation like an alternative generator shape
                    #if (hBDTSystDown.Integral() != 0):
                    #    hBDTSystDown.Scale(hBDT.Integral()/hBDTSystDown.Integral())
                    # hnom + (hnom - halt) = 2*hnom - halt
                    hBDTSystUp = hBDT.Clone("BDT_%s_%s_%sUp" %(catName, sample, syst))
                    hBDTSystUp.Scale(2.0)
                    hBDTSystUp.Add(hBDTSystDown, -1.0)
                    # hnom - (hnom - halt) = halt
                    for ibin in range(hBDTSystUp.GetNbinsX()):
                        B_Up = hBDTSystUp.GetBinContent(ibin)
                        B_Down = hBDTSystDown.GetBinContent(ibin)
                        if (B_Up < 0):
                            B_Up_Err = hBDTSystUp.GetBinError(ibin)
                            hBDTSystUp.SetBinContent(ibin, 0.0)
                            hBDTSystUp.SetBinError(ibin, B_Up + B_Up_Err)
                        if (B_Down < 0):
                            B_Down_Err = hBDTSystDown.GetBinError(ibin)
                            hBDTSystDown.SetBinContent(ibin, 0.0)
                            hBDTSystDown.SetBinError(ibin, B_Down + B_Down_Err)
        elif (sysWeight != "1.0"):
            print "made it here"
            # read in up/down weights from a root file
            weightfilename,weighthistname = sysWeight.split(':')
            weightfile = ROOT.TFile.Open(weightfilename)
            weighthist = weightfile.Get(weighthistname)
            ofile.cd()
            hBDTSystUp = hBDT + (hBDT * weighthist)
            hBDTSystDown = hBDT - (hBDT * weighthist) 
            # name gets set to the name of hBDT when you do as above, need to change it back
            hBDTSystUp.SetName("BDT_%s_%s_%sUp" % (catName,sample,syst))
            hBDTSystDown.SetName("BDT_%s_%s_%sDown" % (catName,sample,syst))
            print hBDTSystUp.GetName(),hBDTSystDown.GetName()
            weightfile.Close()
        else:
            #tree.Draw("%s>>%s_%sUp" % (sysBDTNameUp, sample, syst),"(%s)*weight*(2.2/1.28)" % (cutString))   
            #tree.Draw("%s>>%s_%sUp" % (sysBDTNameUp, sample, syst),"((%s)&&%s)*%s" % (cutString,passSys,weight_string))   
            print passSysup
            tree.Draw("%s>>BDT_%s_%s_%sUp" % (sysBDTNameUp, catName,sample, syst),"((%s)&&%s)*%s" % (cutStringup,passSysup,weight_string))   
            #tree.Draw("%s>>%s_%sDown" % (sysBDTNameDown, sample, syst),"(%s)*weight*(2.2/1.28)" % (cutString))   
            #tree.Draw("%s>>%s_%sDown" % (sysBDTNameDown, sample, syst),"((%s)&&%s)*%s" % (cutString,passSys,weight_string))   
            print passSysdown
            tree.Draw("%s>>BDT_%s_%s_%sDown" % (sysBDTNameDown, catName,sample, syst),"((%s)&&%s)*%s" % (cutStringdown,passSysdown,weight_string))   
            if (syst.find("Model") != -1):
                # shape from different sample (amc@NLO vs. POWHEG, for example)
                cutStringAltModel = makeCutString(sample,sampleMapAltModel)
                hBDTAltModel = ROOT.TH1F("%s_%sAltModel" % (sample,syst), "%s_%sAltModel" % (sample,syst),nBins,binBoundaries)
                tree.Draw("%s>>%s_%sAltModel" % (sysBDTNameUp, sample, syst),"((%s)&&%s)*%s" % (cutStringAltModel,passSysup,weight_string))
                if (hBDTAltModel.Integral() != 0):
                    hBDTAltModel.Scale(hBDT.Integral()/hBDTAltModel.Integral())
                # hnom + (hnom - halt) = 2*hnom - halt
                hBDTSystUp.Scale(2.0)
                hBDTSystUp.Add(hBDTAltModel, -1.0)
                # hnom - (hnom - halt) = halt
                hBDTSystDown = hBDTAltModel
                hBDTSystDown.SetName("BDT_%s_%s_%sDown" % (catName,sample, syst))
                for ibin in range(hBDTSystUp.GetNbinsX()):
                    B_Up = hBDTSystUp.GetBinContent(ibin)
                    B_Down = hBDTSystDown.GetBinContent(ibin)
                    if (B_Up < 0):
                        B_Up_Err = hBDTSystUp.GetBinError(ibin)
                        hBDTSystUp.SetBinContent(ibin, 0.0)
                        hBDTSystUp.SetBinError(ibin, B_Up + B_Up_Err)
                    if (B_Down < 0):
                        B_Down_Err = hBDTSystDown.GetBinError(ibin)
                        hBDTSystDown.SetBinContent(ibin, 0.0)
                        hBDTSystDown.SetBinError(ibin, B_Down + B_Down_Err)
        # Combine cannot handle shapes with negative norms
        if (hBDT.Integral() < 0):
            hBDT.Scale(-0.00001) ## if it was just 0.0 we would lose the shape
        if (hBDTSystUp.Integral() < 0):
            hBDTSystUp.Scale(-0.00001) ## if it was just 0.0 we would lose the shape
        if (hBDTSystDown.Integral() < 0):
            hBDTSystDown.Scale(-0.00001) ## if it was just 0.0 we would lose the shape

        ofile.cd()
        hBDTSystUp.Write()
        hBDTSystDown.Write()
        print syst
        print sysWeight
        print hBDT.Integral()
        print hBDTSystUp.Integral(),hBDTSystDown.Integral()
        hBDTSystUp.Delete()
        hBDTSystDown.Delete()
    ofile.cd()
    if (sample != "data_obs"):
        print True
    else:
        print False
    print args.doData
    print sample
    if (sample != "data_obs" or args.doData):
        print "doing this"
        hBDT.Write("BDT_%s_%s" % (catName,sample))
        print hBDT.Integral()
    else: print hBkg.Integral()
    print bdtname
    print "((%s)&&Pass_nominal)*%s" % (cutString,weight_string)
    if (sample == "TT" and args.ttbarTree != ""): 
        ifile_tt.Close()
    tree.Reset()
#ifile.Close()
ofile.Close()
otextfile.close()
   
