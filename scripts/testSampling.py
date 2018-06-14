import ROOT
f1 = ROOT.TFile.Open("/eos/uscms/store/group/lpchbb/VHbbAnalysisNtuples/2017V5ResampTest2/haddjobs/sum_WJets-HT200To400.root")
f2 = ROOT.TFile.Open("/eos/uscms/store/group/lpchbb/VHbbAnalysisNtuples/2017V5ResampTest2Nom/haddjobs/sum_WJets-HT200To400.root")
tree1 = f1.Get("Events")
tree2 = f2.Get("Events")
ROOT.gStyle.SetOptStat(0)

#vars_to_plot = [("CMS_vhbb_DNN_Wmn_13TeV",30,0,1),("CMS_vhbb_BDTG_Wln_13TeV",20,-1,1),("V_pt",20,100.,300.),("hJets_btagged_0",15,0.,1.),("hJets_btagged_1",15,0.,1.)]
vars_to_plot = [("V_pt",20,100.,300.),("hJets_btagged_0",15,0.,1.),("hJets_btagged_1",15,0.,1.)]
canv = ROOT.TCanvas("canv")

for var,nbins,xmin,xmax in vars_to_plot:
    print var
    h1 = ROOT.TH1F("h1","h1",nbins,xmin,xmax)
    h2 = ROOT.TH1F("h2","h2",nbins,xmin,xmax)
    tree1.Draw("%s>>h1" %var,"(hJets_btagged_0>0.8001&&hJets_btagged_1>0.1522&&twoResolvedJets&&V_pt>=150&&(isWmunu||isWenu)&&controlSample==0&&H_mass>90&&H_mass<150)*weight")
    tree2.Draw("%s>>h2" %var,"(hJets_btagged_0>0.8001&&hJets_btagged_1>0.1522&&twoResolvedJets&&V_pt>=150&&(isWmunu||isWenu)&&controlSample==0&&H_mass>90&&H_mass<150)*weight")
    print h1.Integral(),h2.Integral()
    h1.Scale(1./h1.Integral())
    h2.Scale(1./h2.Integral())
    h2.SetTitle("")
    h2.Draw("ep")
    h1.Draw("hist same")
    leg = ROOT.TLegend(0.2,0.2,0.5,0.5)
    leg.AddEntry(h1,"W+jets HT200-400 Resampled")
    leg.AddEntry(h2,"W+jets HT200-400 Nominal")
    #leg.Draw("same")
    canv.SetLogy(False)
    canv.SaveAs("rsvalplots/%s.png" % var)
    canv.SetLogy(True)
    canv.SaveAs("rsvalplots/%s_log.png" % var)
