import ROOT
import sys
import glob
import array

directory=sys.argv[1]

def PrintAllBins(hist):
    NbinX=hist.GetNbinsX()
    NbinY=hist.GetNbinsY()
    print "\t",
    for iY in range(NbinY+1):
        print iY,
    print
    for iX in range(NbinX+1):
        print iX,
        for iY in range(NbinY+1):
            print hist.GetBinContent(iX,iY),
        print

year=2016
#if len(sys.argv)>3:
#    year=int(sys.argv[2])
if year==2018:
    tuneKey="TuneCP5"
elif year==2017:
    tuneKey="TuneCP5"
elif year==2016:
    tuneKey="TuneCUETP8M1"
else:
    print "what year is ",year,"?"
    sys.exit()


listOfSamples={}
listOfSamples[0]={}
listOfSamples[0]["HT"]=["ZJetsToNuNu_HT-100To200_13TeV-madgraph","ZJetsToNuNu_HT-200To400_13TeV-madgraph","ZJetsToNuNu_HT-400To600_13TeV-madgraph","ZJetsToNuNu_HT-600To800_13TeV-madgraph","ZJetsToNuNu_HT-800To1200_13TeV-madgraph","ZJetsToNuNu_HT-1200To2500_13TeV-madgraph","ZJetsToNuNu_HT-2500ToInf_13TeV-madgraph"]
listOfSamples[0]["BENR"]=["ZJetsToNuNu_BGenFilter_Zpt-100to200_"+tuneKey+"_13TeV-madgraphMLM-pythia8","ZJetsToNuNu_BGenFilter_Zpt-200toInf_"+tuneKey+"_13TeV-madgraphMLM-pythia8","ZBJetsToNuNu_Zpt-100to200_"+tuneKey+"_13TeV-madgraphMLM-pythia8","ZBJetsToNuNu_Zpt-200toInf_"+tuneKey+"_13TeV-madgraphMLM-pythia8"]

listOfSamples[1]={}
listOfSamples[1]["Inclu"]=["WJetsToLNu_"+tuneKey+"_13TeV-madgraphMLM-pythia8"]
listOfSamples[1]["HT"]=["WJetsToLNu_HT-100To200_"+tuneKey+"_13TeV-madgraphMLM-pythia8","WJetsToLNu_HT-200To400_"+tuneKey+"_13TeV-madgraphMLM-pythia8","WJetsToLNu_HT-400To600_"+tuneKey+"_13TeV-madgraphMLM-pythia8","WJetsToLNu_HT-600To800_"+tuneKey+"_13TeV-madgraphMLM-pythia8","WJetsToLNu_HT-800To1200_"+tuneKey+"_13TeV-madgraphMLM-pythia8","WJetsToLNu_HT-1200To2500_"+tuneKey+"_13TeV-madgraphMLM-pythia8","WJetsToLNu_HT-2500ToInf_"+tuneKey+"_13TeV-madgraphMLM-pythia8"]
listOfSamples[1]["BENR"]=["WBJetsToLNu_Wpt-100to200_"+tuneKey+"_13TeV-madgraphMLM-pythia8","WBJetsToLNu_Wpt-200toInf_"+tuneKey+"_13TeV-madgraphMLM-pythia8","WJetsToLNu_BGenFilter_Wpt-100to200_"+tuneKey+"_13TeV-madgraphMLM-pythia8","WJetsToLNu_BGenFilter_Wpt-200toInf_"+tuneKey+"_13TeV-madgraphMLM-pythia8"]



listOfSamples[2]={}
if year == 2018:
    listOfSamples[2]["BENR"]=["DYBJetsToLL_M-50_Zpt-100to200_"+tuneKey+"_13TeV-madgraphMLM-pythia8","DYBJetsToLL_M-50_Zpt-200toInf_"+tuneKey+"_13TeV-madgraphMLM-pythia8","DYJetsToLL_BGenFilter_Zpt-100to200_M-50_"+tuneKey+"_13TeV-madgraphMLM-pythia8","DYJetsToLL_BGenFilter_Zpt-200toInf_M-50_"+tuneKey+"_13TeV-madgraphMLM-pythia8"]
    listOfSamples[2]["HT"]=["DYJetsToLL_M-50_HT-70to100_"+tuneKey+"_PSweights_13TeV-madgraphMLM-pythia8","DYJetsToLL_M-50_HT-100to200_"+tuneKey+"_PSweights_13TeV-madgraphMLM-pythia8","DYJetsToLL_M-50_HT-200to400_"+tuneKey+"_PSweights_13TeV-madgraphMLM-pythia8","DYJetsToLL_M-50_HT-400to600_"+tuneKey+"_PSweights_13TeV-madgraphMLM-pythia8","DYJetsToLL_M-50_HT-600to800_"+tuneKey+"_PSweights_13TeV-madgraphMLM-pythia8","DYJetsToLL_M-50_HT-800to1200_"+tuneKey+"_PSweights_13TeV-madgraphMLM-pythia8","DYJetsToLL_M-50_HT-1200to2500_"+tuneKey+"_PSweights_13TeV-madgraphMLM-pythia8","DYJetsToLL_M-50_HT-2500toInf_"+tuneKey+"_PSweights_13TeV-madgraphMLM-pythia8"]
else:
    listOfSamples[2]["BENR"]=["DYBJetsToLL_M-50_Zpt-100to200_"+tuneKey+"_13TeV-madgraphMLM-pythia8","DYBJetsToLL_M-50_Zpt-200toInf_"+tuneKey+"_13TeV-madgraphMLM-pythia8","DYJetsToLL_BGenFilter_Zpt-100to200_M-50_"+tuneKey+"_13TeV-madgraphMLM-pythia8","DYJetsToLL_BGenFilter_Zpt-200toInf_M-50_"+tuneKey+"_13TeV-madgraphMLM-pythia8"]
    listOfSamples[2]["HT"]=["DYJetsToLL_M-50_HT-70to100_"+tuneKey+"_13TeV-madgraphMLM-pythia8","DYJetsToLL_M-50_HT-100to200_"+tuneKey+"_13TeV-madgraphMLM-pythia8","DYJetsToLL_M-50_HT-200to400_"+tuneKey+"_13TeV-madgraphMLM-pythia8","DYJetsToLL_M-50_HT-400to600_"+tuneKey+"_13TeV-madgraphMLM-pythia8","DYJetsToLL_M-50_HT-600to800_"+tuneKey+"_13TeV-madgraphMLM-pythia8","DYJetsToLL_M-50_HT-800to1200_"+tuneKey+"_13TeV-madgraphMLM-pythia8","DYJetsToLL_M-50_HT-1200to2500_"+tuneKey+"_13TeV-madgraphMLM-pythia8","DYJetsToLL_M-50_HT-2500toInf_"+tuneKey+"_13TeV-madgraphMLM-pythia8"]

listOfSamples[2]["Inclu"]=["DYJetsToLL_M-50_"+tuneKey+"_13TeV-madgraphMLM-pythia8"]

if year==2017:
    for iChan in range(3):
        for iName in range(len(listOfSamples[iChan]["BENR"])):
            listOfSamples[iChan]["BENR"][iName]=listOfSamples[iChan]["BENR"][iName]+"_newgridpack"

enrichmentCats="((LHE_Nb>0&&nGenStatus2bHad>0&&LHE_Vpt<200&&LHE_Vpt>100)*1+(LHE_Nb>0&&nGenStatus2bHad>0&&LHE_Vpt>=200)*2 +(LHE_Nb==0&&nGenStatus2bHad>0&&LHE_Vpt<200&&LHE_Vpt>100)*3+(LHE_Nb==0&&nGenStatus2bHad>0&&LHE_Vpt>=200)*4 - 1)"
HTbins=[0,70,100,200,400,600,800,1200,2500,100000]
HTbins_array=array.array("d",HTbins)
print HTbins_array
#HTBinString="-1 + "
#for iBin in len(HTbins):
#    thisBinString="(LHE_HT>="+str(HTbin)
#    if HTbin is HTbins[-1]:
#        thisBinString=thisBinString+")"
#    else:
#        thisBinString=thisBinString+"&&LHE_HT<"+str(HTbin)")"
#    HTBinString=HTBinString+thisBinString+"*"+str(iBin+1)

hists={}
ratios={}
can=ROOT.TCanvas("can","",700,700)
can.SetLogy(1)
ROOT.gStyle.SetOptStat(0)
channels=listOfSamples.keys()

for chan in channels:
    print "channel",chan
    for sampleType in listOfSamples[chan]:
        print "type",sampleType
        histName=sampleType+"_BenrVsHTBin_chan"+str(chan)
        hists[histName]=ROOT.TH2F(histName,";b-enrichment bin;HT bin",5,-1,4,len(HTbins)-1,HTbins_array)
        for sample in listOfSamples[chan][sampleType]:
            print sample
            files=glob.glob(directory+"/"+sample+"/*/*/*/*.root")
            print len(files)
            histName=sample+"_BenrVsHTBin"
            hists[histName]=ROOT.TH2F(histName,";b-enrichment bin;HT bin",5,-1,4,len(HTbins)-1,HTbins_array)
            localHist=ROOT.TH2F("localHist",";b-enrichment bin;HT bin",5,-1,4,len(HTbins)-1,HTbins_array)
            
            iFile=0
            for fileName in files:
                tfile=ROOT.TFile.Open(fileName)
                tree=tfile.Get("Events")
                print iFile,tree.GetEntries(),
                #hists[histName].SetDirectory(ROOT.gDirectory)
                localHist.SetDirectory(ROOT.gDirectory)
                tree.Draw("LHE_HT:"+enrichmentCats+">>localHist")
                localHist.SetDirectory(0)
                print localHist.GetEntries(),localHist.Integral(),
                tfile.Close()
                hists[histName].Add(localHist)
                print hists[histName].GetEntries(),hists[histName].Integral()
                iFile=iFile+1

            hists[histName].GetYaxis().SetRangeUser(5,HTbins[-1])
            hists[histName].Draw("colz")
            can.Update()
            can.SaveAs(histName+".png")
            #PrintAllBins(hists[histName])
            hists[sampleType+"_BenrVsHTBin_chan"+str(chan)].Add(hists[histName])


    totalName="totalChan"+str(chan)
    hists[totalName]=ROOT.TH2F(totalName,";b-enrichment bin;HT bin",5,-1,4,len(HTbins)-1,HTbins_array)
    for sampleType in listOfSamples[chan]:
        hists[totalName].Add(hists[sampleType+"_BenrVsHTBin_chan"+str(chan)])
        ratioName="ratio"+sampleType+"Chan"+str(chan)
        ratios[ratioName]=ROOT.TH2F(ratioName,"HT/(HT+BENR);b-enrichment bin;HT bin",5,-1,4,len(HTbins)-1,HTbins_array)
        ratios[ratioName].Add(hists[sampleType+"_BenrVsHTBin_chan"+str(chan)])
        
    for sampleType in listOfSamples[chan]:
        ratioName="ratio"+sampleType+"Chan"+str(chan)
        ratios[ratioName].Divide(hists[totalName])
        ratios[ratioName].Draw("text colz")
        can.Update()
        can.SaveAs(ratioName+".png")
    
    hists[totalName].Draw("colz")
    can.Update()
    can.SaveAs(totalName+".png")


newTFile=ROOT.TFile.Open("output.root","RECREATE")
for histName in ratios.keys():
    ratios[histName].Write()
newTFile.Close()
