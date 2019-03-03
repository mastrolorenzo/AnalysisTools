//
//  Written by Chris Palmer
//
//  VHbb analysis
//  Inheriting from Analysis Manager
//

#ifndef ANALYSISTOOLS_PLUGINS_VHBBANALYSIS_H_
#define ANALYSISTOOLS_PLUGINS_VHBBANALYSIS_H_

#include <iostream>
#include "TLorentzVector.h"
#include "../AnalysisManager.h"
#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "RooFunctor.h"


class VHbbAnalysis : public AnalysisManager {

    public:
        VHbbAnalysis();
        virtual ~VHbbAnalysis();

        virtual void InitAnalysis();
        virtual bool Preselection();
        virtual bool Analyze();
        virtual void FinishEvent();
        virtual void TermAnalysis(bool skim=false);

        std::string taggerName;
        void SetTaggerName(float taggerType);

        TH1D* ewkCorrHist_wp;
        TH1D* ewkCorrHist_wm;
        TH1D* ewkCorrHist_zll;
        TH1D* ewkCorrHist_znn;

        RooFunctor* met_trigger_sf120_func;
        RooFunctor* met_trigger_eff_2016data_func;

        void UpdateJetPts();

        std::pair<int,int> HighestPtGoodElectronsOppCharge(float min_pt, float max_rel_iso, float idcut, bool isOneLepton);
        std::pair<int,int> HighestPtGoodMuonsOppCharge(float min_pt, float max_rel_iso,bool isOneLepton=false);
        bool ElectronSelection(int);
        bool MuonSelection(int);
        int UpdatedVType();
        bool PassVTypeAndTrigger(int vtype);
        bool SelectLeptonChannel();
        bool SelectJets();
        bool ReconstructHiggsCand();
        bool ReconstructVCand();
        void ComputeVHKinematics();
        void ComputeOtherEventKinematics();
        void ControlSampleSelection();
        float ReWeightMC(int nPU=0);
        float puWeight_2016(int i=0);
        float puWeight_2016Up(int i=0);
        float puWeight_2016Down(int i=0);
        float puWeight_2016to2017(int i=0);
        float puWeight_ichep(int i=0);
        float puWeight_ichep_up(int i=0);
        float puWeight_ichep_down(int i=0);
        std::pair<int,int> HighestPtBJets();
        std::pair<int,int> HighestTaggerValueBJets(float j1ptCut, float j2ptCut, std::string taggerName);
        std::pair<int,int> HighestDeepCSVBJets(float j1ptCut, float j2ptCut);
        std::pair<int,int> HighestCMVABJets(float j1ptCut, float j2ptCut);
        std::pair<int,int> HighestCSVBJets(float j1ptCut, float j2ptCut);
        std::pair<int,int> HighestPtJJBJets();
        double GetRecoTopMass(TLorentzVector Jet, bool isJet=true, int useMET=0, bool regPT=true, bool smearedPT=true);
        float ptWeightQCD(int nGenVbosons=0, float lheHT=0., int GenVbosons_pdgId=0);
        float ptWeightEWK(int nGenVbosons=0, float GenVbosons_pt=0., int VtypeSim=0, int GenVbosons_pdgId=0);
        TLorentzVector getNu4Momentum(const TLorentzVector& TLepton, const TLorentzVector& TMET);
        double LOtoNLOWeightBjetSplitEtabb(double etabb=10., int njets=0);
        float GetVHEWKCorrFactor( float V_pt, TH1D* hist );
        float getVPtCorrFactor(float V_pt=0., int si=0, int sysVar=0);
        void  smearJets(float JERScale=1.0);
        float evaluateRegression(int i=0);
        void SetupFactorizedJECs(std::string variation="nominal");

        void FatJetSelection();
        void BoostedSelection();
        void ComputeBoostedVariables();
        bool atLeastOnePreselFatJet;
	bool enableFSRRecovery;

        TLorentzVector HJ1, HJ2, Hbb;
        TLorentzVector HJ1_noFSR, HJ2_noFSR, Hbb_noFSR;
        TLorentzVector HJ1_noreg, HJ2_noreg, Hbb_noreg;
        TLorentzVector fatJetCand;
        TLorentzVector V;
        TLorentzVector W_withNuFromMWCon; //I think we could make this the V in the 1-lepton case, but I didn't want to change logic in this PR.
};

#endif // ANALYSISTOOLS_PLUGINS_VHBBANALYSIS_H_
