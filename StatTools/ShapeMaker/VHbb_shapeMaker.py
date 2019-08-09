# ======================================================
# Shape Maker for Datacards, for VHcc 2016 Analysis
# Author: Spandan Mondal
# RWTH Aachen University, Germany
#
# This code aims to optimize file access by adopting
# a per-file processing approach. It implements
# ROOT RDataFrame to create histograms (instead of
# TTree::Draw) for all pocesses, channels and
# systematics, so that the program loops only once
# over the input chain, thereby reducing disk access.
# ======================================================

from datetime import datetime
launchtime = datetime.now()
lasttime = launchtime

import rootpy
from ROOT import *
gROOT.SetBatch(True)
import argparse, os, sys
from array import array
#from Dictionaries.VHbb_fileDict import fileDict
from Dictionaries.VHbb_fileDict import Zll_fileDict, Wln_fileDict, Znn_fileDict
from Dictionaries.VHbb_channelDict import channelDict
from Dictionaries.VHbb_sampleDict import Zll_sampleDict, Wln_sampleDict, Znn_sampleDict
from Dictionaries.VHbb_binsDict import binsDict


fileDict = dict(Zll_fileDict.items() + Wln_fileDict.items() + Znn_fileDict.items())

usingRDF = False
rootver = gROOT.GetVersion()
if '/' in rootver: rootver = rootver.split('/')[0]
if float(rootver) >= 6.14:
    TH1DModel = ROOT.RDF.TH1DModel
    ROOTDataFrame = ROOT.RDataFrame
    print "Detected ROOT %s. Successfully loaded RDataFrame."%rootver
    usingRDF = True
else:
    TH1DModel = ROOT.Experimental.TDF.TH1DModel
    ROOTDataFrame = ROOT.Experimental.TDataFrame
    print "***** WARNING *****: Detected ROOT %s. Loaded TDataFrame instead of RDataFrame. Multithreading will also be forced off. This is usually slower than RDataFrame. Please use ROOT 6.14 and above for better performance.\n"%rootver 

# Message logger
def elapsed(string,level=0):
    nowtime = datetime.now()
    global lasttime
    print ' '*(level*4)+"T = %.1f mins; dT = %.3f mins: "%(float((nowtime-launchtime).seconds)/60, float((nowtime-lasttime).seconds)/60) + string
    lasttime = nowtime

parser = argparse.ArgumentParser("Produce shape histograms for datacards.")
parser.add_argument('-i', '--inputfile', type=str, default="", help="The input root file (ntuple).")
parser.add_argument('-d', '--inputdir', type=str, default="/nfs/dust/cms/user/lmastrol/VHbbAnalysisNtuples/VHccNtuple_March16_forApproval_fixZnnHF/haddjobs/", help="The directory where all the ntuples are stored.")
parser.add_argument('-samp', '--samplename', type=str, default="", help="The name of the sample, eg, 'DYToLL_HT100to200'.")
parser.add_argument('-n', '--nOuts', type=str, default="0", help="The max. number of output files to produce. Use 0 to process full input file.")
parser.add_argument('-io', '--iOut', type=str, default="", help="The index of the first output file that will be produced if iOut!=0. So this instance of the program will produce ith to (i+n-1)th outputs.")
parser.add_argument('-ld', '--lockDir', type=str, default=".", help="The directory where lock files will be stored to avoid parallel access while reading the same input file")
parser.add_argument('--vptcut', type=str, default="", help="The cut on V_pT, required for combination with boosted analysis.")
parser.add_argument('--skipVZ', action="store_true", default=False, help="Skip VZ selections.")
parser.add_argument('--skipsysts', action="store_true", default=False, help="Skip systematics.")
parser.add_argument('--dochannel',  type=str, default="", help="Do only specified channel.")
parser.add_argument('--multithread', action="store_true", default=False, help="Enable Multithreading for RDataFrame.")
args = parser.parse_args()
print args

if args.multithread and usingRDF: ROOT.EnableImplicitMT(2)

ifile = args.inputfile
ipath = args.inputdir
isamp = args.samplename

ilist = []
if ',' in ifile:
    ilist = ifile.split(',')
    ifile = ilist[0]

nOuts = int(args.nOuts)
if nOuts != 0:
    if args.iOut=="" or int(args.iOut) < 0:
        print "\n***** ERROR *****: Invalid input for iOut.\n"
        sys.exit(1)
    iOut = int(args.iOut)
    lockDir = args.lockDir

if isamp == "" and ifile != "":
    sampMatch = [sampname for sampname in fileDict.keys() if sampname in ifile]
    if len(sampMatch) == 0:
        print "\n***** ERROR *****: No sample matches with the input file, %s.\n"%(ifile)
        sys.exit(1)
    else:
        isamp = sorted(sampMatch,key=len)[-1]
        elapsed("Detected sample: %s"%isamp)
        fileglob = ipath+"/*"+ifile+"*"
elif ifile == "" and isamp != "":
    ifile = isamp
    fileglob = ipath+"/*"+ifile+"*.root"
elif ifile != "" and isamp != "":
    fileglob = ipath+"/*"+ifile+"*"
else:
    print "\n***** ERROR *****: Please provide inputfile and/or samplename as input."
    sys.exit(1)

if not (isamp in ifile):
    print "\n***** WARNING *****: The sample name and file name may not be consistent. I will proceed, but please double check.\n"

outSuffix = ''.join(ifile.rstrip('.root').lstrip('sum_').split('_'))

if ilist==[]:
    elapsed("Using *%s* as input file(s) and %s as input sample name."%(ifile,isamp))
    mainDF = ROOTDataFrame("Events", fileglob)
else:
    elapsed("Using %s as input files and %s as input sample name."%(str(ilist),isamp))
    chain = TChain("Events")
    for rootfile in ilist:
        chain.Add(ipath+"/"+rootfile)
    mainDF = ROOTDataFrame(chain)
elapsed("Created pointer to input files containing %s in the filename."%ifile)

processArrays = fileDict[isamp]
processNames = [i[0] for i in processArrays]
elapsed("I see this file contributes to the following processes: "+str(processNames))

#allChannels = sorted(channelDict.keys())
#only apply channels to the relevant samples
channelList = []
if isamp in Zll_fileDict:
    channelList.extend([channel for channel in channelDict if channelDict[channel]["channelName"]=="Zll"])
if isamp in Wln_fileDict:
    channelList.extend([channel for channel in channelDict if channelDict[channel]["channelName"]=="Wln"])
if isamp in Znn_fileDict:
    channelList.extend([channel for channel in channelDict if channelDict[channel]["channelName"]=="Znn"])

if args.skipVZ: channelList = [i for i in channelList if not i.startswith('VZ')]
elapsed("I am going to run over %i channels: "%(len(channelList)) + str(channelList))
print "Therefore, I will produce %ix%i=%i outputs."%(len(processNames),len(channelList),len(processNames)*len(channelList))
if nOuts != 0: print "However, I see this instance of me has been asked to produce only %i of these outputs starting with index %i."%(nOuts,iOut)

collist = list(mainDF.GetColumnNames())
if len(collist) == 0:
    print "\n***** ERROR *****: Wrong input file. Failed to load any data."
    sys.exit(1)
#if "CMS_vhcc_BDT_Znn_13TeV"         in collist: mainDF = mainDF.Alias("CMS_vhcc_BDTG_Znn_HighPT_13TeV", "CMS_vhcc_BDT_Znn_13TeV")
#if "CMS_vhcc_BDT_Wln_13TeV"         in collist: mainDF = mainDF.Alias("CMS_vhcc_BDTG_Wln_13TeV",        "CMS_vhcc_BDT_Wln_13TeV")
#if "CMS_vhcc_BDT_Zll_HighPT_13TeV"  in collist: mainDF = mainDF.Alias("CMS_vhcc_BDTG_Zll_HighPT_13TeV", "CMS_vhcc_BDT_Zll_HighPT_13TeV")
#if "CMS_vhcc_BDT_Zll_LowPT_13TeV"   in collist: mainDF = mainDF.Alias("CMS_vhcc_BDTG_Zll_LowPT_13TeV",  "CMS_vhcc_BDT_Zll_LowPT_13TeV")

outputHistoDict = {}

nMyCol = 0
columnDict = {}
def addcolumn(formula):
    global mainDF, nMyCol, columnDict
    collist = list(mainDF.GetColumnNames())
    if formula not in columnDict.keys() and formula not in collist:
        nMyCol+=1
        mainDF = mainDF.Define("MyCol%i"%nMyCol,formula)
        columnDict[formula] = "MyCol%i"%nMyCol
    if formula not in columnDict.keys() and formula in collist:
        columnDict[formula] = formula

nFile = -1

for iprocess, processArray in enumerate(processArrays):
    processName = processArray[0]
    processIdx = processArray[1]
    elapsed("(%i of %i) Beginning with %s process."%(iprocess+1,len(processArrays),processName))
    if type(processIdx)==int:
        processCut = "sampleIndex==%i"%processIdx
    elif type(processIdx)==str:
        processCut = processIdx
    else:
        raise Exception("Bad sampleIndex value")

    for ichannel, channel in enumerate(channelList):
        nFile += 1
        if nOuts != 0:
            if nFile < iOut or nFile >= iOut+nOuts: continue
        if args.dochannel != "":
            if not ',' in args.dochannel:
                if channel != args.dochannel: continue
            else:
                channelsToDo = args.dochannel.split(',')
                if channel not in channelsToDo: continue
        
        elapsed("(%i of %i) Processing %s channel."%(ichannel+1,len(channelList),channel),1)
        subdict = channelDict[channel]
        presel = subdict["cutstring"]
        varName = subdict["varName"]
        doEven = bool(subdict["doEven"])
        xlow = float(subdict["xlow"])
        xhigh = float(subdict["xhigh"])
        nBins = int(subdict["nBins"])
        channelName = subdict["channelName"]
        drawFromNom = bool(subdict["drawFromNom"])
        doVV = bool(subdict["doVV"])
        doRebin = bool(subdict["doRebin"])
        
        if args.vptcut == "":
            vptcutStr = ""
        else:
            vptcutStr = " && V_pt < %s"%(args.vptcut)
        presel = presel.replace("VPTCUT",vptcutStr)
        if type(processIdx)==int:
            outfilename = "hists_%s_%s_idx%i_file-%s.root"%(channel,processName,processIdx,outSuffix)
        elif type(processIdx)==str:
            print processIdx
            outfilename = "hists_%s_%s_idx%s_file-%s.root"%(channel,processName,processIdx,outSuffix)
        outputHistoDict[outfilename] = []

        def addHisto(filterStr, hName, hTitle, varToPlot, weight="", info=""):
            global outputHistoDict
            if doRebin:
                binning = binsDict[channel]
                if len(binning)-1 != nBins:
                    print "\n***** WARNING *****: The binning in the file and 'nBins' are not consistent for channel %s. I will proceed, but please double check.\n"%channel
                histoModel = TH1DModel(hName,hTitle,len(binning)-1,array('d',binning))
            else:
                histoModel = TH1DModel(hName,hTitle,nBins,xlow,xhigh)
            if weight=="":
                outputHistoDict[outfilename].append([mainDF.Filter(filterStr)        \
                                                    .Histo1D(   histoModel,
                                                                columnDict[varToPlot]
                                                            ),
                                                    info ]
                                                   )
            else:
                outputHistoDict[outfilename].append([mainDF.Filter(filterStr)        \
                                                    .Histo1D(   histoModel,
                                                                columnDict[varToPlot],
                                                                columnDict[weight]
                                                            ),
                                                    info ]
                                                   )

        if doEven and processName!="data_obs":
            weight_string = "2.0*weight"
            cutstring = presel + " && (event%2==0)"
        else:
            weight_string = "weight"
            cutstring = presel

        cutstring += "&& ("+processCut+")"

        #================================= Nominal =================================
        addcolumn(varName)

        #print cutstring+" && Pass_nominal"
        #print weight_string

        if processName == "data_obs":
            addHisto(cutstring+" && Pass_nominal", "BDT_%s_%s"%(channel,processName), "nominal", varName)
        else:
            addcolumn(weight_string)
            addHisto(cutstring+" && Pass_nominal", "BDT_%s_%s"%(channel,processName), "nominal", varName, weight_string)
        #---------------------------------------------------------------------------
        if args.skipsysts: continue

        if channelName == "Zll":
            sampleDict = Zll_sampleDict
        elif channelName == "Wln":
            sampleDict = Wln_sampleDict
        elif channelName == "Znn":
            sampleDict = Znn_sampleDict

        

        if processName == "data_obs":
            elapsed("Created pointer to nominal histogram. Since this is data, there are no systematics to process.",1)
            continue
        systematicsProcesses = [sampleDict[k][3] for k in sampleDict]
        systematicsProcesses = set([x for sub in systematicsProcesses for x in sub])
        if processName not in systematicsProcesses:
        #if processName not in set([x for x in sampleDict[k][3] for k in sampleDict]):
        #if processName not in sampleDict:
            elapsed("I could not find the systematics relevant to %s channel for %s process. Perhaps this process is not required for this channel...? I have created the nominal histogram, but skipping systematics."%(channel,processName))
            continue
        systListofList = [[x,sampleDict[x][0],sampleDict[x][1],sampleDict[x][2]] for x in sampleDict if processName in sampleDict[x][3]] 
        #systListofList = sampleDict[processName]
        systNameList = [i[0] for i in systListofList]
        elapsed("Created pointer to nominal histogram, next I will plot %i systematics. "%(len(systNameList)),1)

        for isyst, systList in enumerate(systListofList):
            print "\n"
            if "Zee" in channel or "Wen" in channel:
                systName = systList[0].replace("eff_m","eff_e")
            else:
                systName = systList[0]
            elapsed("(%i of %i) Processing %s systematic."%(isyst+1,len(systListofList),systName),2)

            systWeight = systList[3]
            
            #print "systList"
            #print systList
            if '.root' in systWeight:
                elapsed("***** ERROR *****: I have not been programmed to load weights from .root files for systematics. Exiting." )
                # Weights in .root files are not used in VHcc 2016, so I didn't bother to code it... yet. *shrug*
                sys.exit(1)
            if "Model" in systName:
                elapsed("***** ERROR *****: I have not been programmed to handle alternative 'Models' in the systematics. Exiting." )
                # None of the systematics in VHcc 2016 have "Model" in their names
                sys.exit(1)


            bS = ""
            if len(systList) == 5:
                bS = systList[4]
            varNameUp = varName
            varNameDown = varName
            pSUp = "Pass_nominal"
            pSDown = "Pass_nominal"
            cutstringUp = cutstring
            cutstringDown = cutstring
            weight_stringUp = weight_string
            weight_stringDown = weight_string            
            #FIXME when bS was called branchSuffix, the variable would get unset and I have no clue why
            if bS != "":
                if not drawFromNom:
                    if '[' in varName:
                        #bTagWeight_HFUp_pt1_eta0 et al break this convention
                        #btagweight_XUp_... or btagweight_XDown_...
                        varNameUp = varName.replace('[',"_"+bS+"Up[")
                        varNameDown = varName.replace('[',"_"+bS+"Down[")


                    if '_13TeV' in varName:
                        varNameUp = varName.replace('_13TeV',"_13TeV_"+bS+"Up")
                        varNameDown = varName.replace('_13TeV',"_13TeV_"+bS+"Down")
                    if '[' not in varName and '_13TeV' not in varName:
                        varNameUp = varName + "_%sUp"%bS
                        varNameDown = varName + "_%sDown"%bS

                pSUp = "Pass_%sUp"%bS
                pSDown = "Pass_%sDown"%bS

                for var in ["controlSample", "Jet_Pt", "hJetInd1", "hJetInd2",
                            "HVdPhi", "HVdR", "H_pt", "V_pt", "MET_Pt",
                            "H_mass", "H_mass_fit_fallback",
                            "JetsCloseToMET", "dPhi_MET_TkMET",
                            "htJet30", "HJ1_HJ2_dR",
                            "nAddJetsFSRsub302p5_puid"]:
                    cutstringUp = cutstringUp.replace(var, var+"_%sUp"%bS)
                    cutstringDown = cutstringDown.replace(var, var+"_%sDown"%bS)

                weight_stringUp = weight_string.replace("weight","weight_%sUp"%bS)
                weight_stringDown = weight_string.replace("weight","weight_%sDown"%bS)

            else:
            # systWeight is not 1.0
                if ',' not in systWeight:
                    # up and down branch are not separate
                    if systWeight == "weight_PU":
                        weight_stringUp = "%s*(1./weight_PU)*(%sUp)"%(weight_stringUp,systWeight)
                        weight_stringDown = "%s*(1./weight_PU)*(%sDown)"%(weight_stringDown,systWeight)
                    elif systWeight == "WJetNLOWeight":
                        weight_stringUp = "%s*(WJetNLOWeight/1.153)"%weight_stringUp
                        weight_stringDown = "%s*(1.153/WJetNLOWeight)"%weight_stringDown
                        #Renormalizing to nominal.Integral() is done while writing the histograms (at the end of this code).
                    elif systWeight == "recoZReWeight":
                        weight_stringUp = "%s*(1./recoZReWeight)*(%sUp)"%(weight_stringUp,systWeight)
                        weight_stringDown = "%s*(1./recoZReWeight)*(%sDown)"%(weight_stringDown,systWeight)
                    elif systWeight == "recoWReWeight":
                        weight_stringUp = "%s*(1./recoWReWeight)*(%sUp)"%(weight_stringUp,systWeight)
                        weight_stringDown = "%s*(1./recoWReWeight)*(%sDown)"%(weight_stringDown,systWeight)
                    #check if weight is of form btagweight_ptX_etaY and correct name to btagweightUp_ptX_etaY
                    elif systWeight[-8:-6] == "pt" and systWeight[-4:-1] == "eta":
                        #e.g. systWeight: bTagWeight_JES_pt0_eta0+Up
                        #to bTagWeight_JESUp_pt0_eta0
                        #so if _ptX_etaY in systWeight, remove suffix, add Up, reappend
                        eta_pt_suffix = systWeight[-9:]
                        systWeight_root = systWeight[:-9]

                        weight_stringUp = "%s*(%sUp%s)"%(weight_stringUp, systWeight_root,eta_pt_suffix)
                        weight_stringDown = "%s*(%sDown%s)"%(weight_stringDown, systWeight_root,eta_pt_suffix)

                    else:
                        weight_stringUp = "%s*(%sUp)"%(weight_stringUp,systWeight)
                        weight_stringDown = "%s*(%sDown)"%(weight_stringDown,systWeight)

                else:
                    # up and down branch may be separate
                    weight_stringUp = "%s*(%s)"%(weight_stringUp,systWeight.split(',')[1])
                    weight_stringDown = "%s*(%s)"%(weight_stringDown,systWeight.split(',')[0])
            #print "branchSuffix " + (bS)
            #print "varname " + (varName)
            #print "varnameUp " + (varNameUp)
            #print "varnameDown " + (varNameDown)
            #print "weight_stringUp " + (weight_stringUp)
            #print "weight_stringDown " + (weight_stringDown)
            #print "systWeight " + (systWeight)
            addcolumn(varNameUp)
            addcolumn(varNameDown)
            addcolumn(weight_stringUp)
            addcolumn(weight_stringDown)
            addHisto(cutstringUp+" && "+pSUp, "BDT_%s_%s_%sUp"%(channel,processName,systName), systName+"Up", varNameUp, weight_stringUp, systWeight)
            addHisto(cutstringDown+" && "+pSDown, "BDT_%s_%s_%sDown"%(channel,processName,systName), systName+"Down", varNameDown, weight_stringDown, systWeight)
            elapsed("Created pointers to both Up and Down histograms for this systematic.",2)

print
elapsed("==== Now triggering the actual loop over all events. This will take a while. ====")

nEvents = int(mainDF.Count())   # The loop gets triggered while assigning nEvents.
elapsed("Done! I looped over %i events and produced all the histograms."%nEvents)

for outfilename in outputHistoDict:
    ofile = TFile.Open(outfilename,"RECREATE")
    ofile.cd()
    histcount = 0
    for histtuple in outputHistoDict[outfilename]:
        hist = histtuple[0]
        if hist.GetTitle() == "nominal": histNominal = hist
        if hist.Integral() < 0:
            hist.Scale(-0.00001)
        if histtuple[1]=="WJetNLOWeight" and hist.Integral() > 0:
            hist.Scale(histNominal.Integral()/hist.Integral())
            elapsed("Reweighted %s histogram."%(hist.GetName()),1)
        hist.Write()
        histcount += 1
    ofile.Close()
    elapsed("Written %i histograms to %s."%(histcount,outfilename))
elapsed("I produced a total of %i output .root files. And now my watch is ended."%(len(outputHistoDict.keys())))
