//
//  Written by Chris Palmer, Stephane Cooperstein
//
//  VHbb analysis
//  Inheriting from Analysis Manager
//

#include "VHbbAnalysis.h"
#include "HelperClasses/EquationSolver.h"

// initialize parameters
VHbbAnalysis::VHbbAnalysis() {
    if (debug > 10) std::cout << "Constructing VHbbAnalysis" << std::endl;
}

// remove any 'new'
VHbbAnalysis::~VHbbAnalysis() {
    if (debug > 10) std::cout << "Deconstructing VHbbAnalysis" << std::endl;
}

//probably want to do bdtSetup here
void VHbbAnalysis::InitAnalysis() {
    //SetupBDT();
    /*Get the corrections workspace */
    TFile* wsp_file = new TFile("./aux/vhbb_metsf.root","READ");
    RooWorkspace *wsp = (RooWorkspace*)wsp_file->Get("w");
    wsp_file->Close();
    met_trigger_sf120_2018_func = wsp->function("met_trig_sf120_2018")->functor(wsp->argSet("met_mht"));
    met_trigger_sf120_2018_func_up = wsp->function("met_trig_sf120_up_2018")->functor(wsp->argSet("met_mht"));
    met_trigger_sf120_2018_func_down = wsp->function("met_trig_sf120_down_2018")->functor(wsp->argSet("met_mht"));
    met_trigger_sf120_2017_func = wsp->function("met_trig_sf120_2017")->functor(wsp->argSet("met_mht"));
    met_trigger_sf120_2017_func_up = wsp->function("met_trig_sf120_up_2017")->functor(wsp->argSet("met_mht"));
    met_trigger_sf120_2017_func_down = wsp->function("met_trig_sf120_down_2017")->functor(wsp->argSet("met_mht"));
    met_trigger_sf110OR170_2016_func = wsp->function("met_trig_sf110OR170_2016")->functor(wsp->argSet("met_mht"));
    met_trigger_sf110OR170_2016_func_up = wsp->function("met_trig_sf110OR170_up_2016")->functor(wsp->argSet("met_mht"));
    met_trigger_sf110OR170_2016_func_down = wsp->function("met_trig_sf110OR170_down_2016")->functor(wsp->argSet("met_mht"));

    /* Open the files with the EWK correction factor */
    TFile* wpfile = new TFile("./aux/Wp_nloEWK_weight_unnormalized.root","READ");
    ewkCorrHist_wp = (TH1D*)wpfile->Get("SignalWeight_nloEWK");
    TFile* wmfile = new TFile("./aux/Wm_nloEWK_weight_unnormalized.root","READ");
    ewkCorrHist_wm = (TH1D*)wmfile->Get("SignalWeight_nloEWK");
    TFile* zllfile = new TFile("./aux/Zll_nloEWK_weight_unnormalized.root","READ");
    ewkCorrHist_zll = (TH1D*)zllfile->Get("SignalWeight_nloEWK");
    TFile* znnfile = new TFile("./aux/Znn_nloEWK_weight_unnormalized.root","READ");
    ewkCorrHist_znn = (TH1D*)znnfile->Get("SignalWeight_nloEWK");

    ewkCorrHist_wp->Rebin(4);
    ewkCorrHist_wp->Scale(1./4.);
    ewkCorrHist_wm->Rebin(4);
    ewkCorrHist_wm->Scale(1./4.);
    ewkCorrHist_zll->Rebin(4);
    ewkCorrHist_zll->Scale(1./4.);
    ewkCorrHist_znn->Rebin(4);
    ewkCorrHist_znn->Scale(1./4.);

    return;
}

// Check if events pass preselection.
// Returns true for events passing the preselection and false otherwise.
bool VHbbAnalysis::Preselection() {
    bool doCutFlowInPresel = int(m("doCutFlow")) < 0;

    //Set the b-tagger
    SetTaggerName(m("taggerType"));

    if (m("onlyEvenEvents") && m("onlyOddEvents")) {
        std::cout << "Cannot set both onlyEvenEvents and onlyOddEvents to true!!" << std::endl;
        return false;
    }

    if (m("onlyEvenEvents")) {
        if (mInt("event") % 2 == 1) return false;
    } else if (m("onlyOddEvents")) {
        if (mInt("event") % 2 == 0) return false;
    }

    // stitch WJets inclusive sample to HT-binned samples
    //if (cursample->sampleNum == 40 && m("LHE_HT") > 100 ) return false;
    
    // stitch ZJets inclusive sample to HT-binned samples
    // Asking the inclusive sample to cover the whole phase space for 2018 
    
    if(m("dataYear") !=2018)
      {
	if (cursample->sampleNum == 110 && m("LHE_HT") > 100) return false;
      }
    
    if(m("dataYear") ==2018)
      {
	if (cursample->sampleNum == 110 && ( (m("LHE_HT") > 100 && m("LHE_HT") <= 200) || (m("LHE_HT") > 400 && m("LHE_HT") <= 1200) ) ) return false;
      }


    // use W+jets b-enriched samples but make sure all samples are orthogonal
    int nGenStatus2bHad = 0;
    if ( (cursample->sampleNum >= 40 && cursample->sampleNum <=54) || (cursample->sampleNum>=150 && cursample->sampleNum<=163) || (cursample->sampleNum>=110 && cursample->sampleNum<=142)){
        for(int indGP=0; indGP<mInt("nGenPart"); indGP++){
            if(mInt("GenPart_status",indGP)!=2) continue;
            if(((std::abs(mInt("GenPart_pdgId",indGP))/100)%10 ==5) || ((std::abs(mInt("GenPart_pdgId",indGP))/1000)%10==5)){
              nGenStatus2bHad+=1;
            }
        }
    }
    if(m("dataYear")==2016){
        if (cursample->sampleNum >= 40 && cursample->sampleNum <=47) {
            if (m("LHE_Vpt") > 100) {
                if (mInt("LHE_Nb") != 0 || nGenStatus2bHad != 0) return false;
            }
        } else if (cursample->sampleNum == 50) {
            if (m("LHE_Vpt") < 100 || m("LHE_Vpt") > 200 || mInt("LHE_Nb") == 0) return false;
        } else if (cursample->sampleNum == 51) {
            if (m("LHE_Vpt") < 200 || mInt("LHE_Nb") == 0) return false;
        } else if (cursample->sampleNum == 53) {
            //if (m("LHE_Vpt") < 100 || m("LHE_Vpt") > 200 ) return false;
            if (m("LHE_Vpt") < 100 || m("LHE_Vpt") > 200 || nGenStatus2bHad == 0) return false;
        } else if (cursample->sampleNum == 54) {
            //if (m("LHE_Vpt") < 200) return false;
            if (m("LHE_Vpt") < 200 || nGenStatus2bHad == 0) return false;
        }
        if (cursample->sampleNum >= 110 && cursample->sampleNum<=117){
            if (m("LHE_Vpt") > 100) {
                if (mInt("LHE_Nb") != 0 || nGenStatus2bHad != 0) return false;
            }
        } else if (cursample->sampleNum == 121){
            if (m("LHE_Vpt") < 100 || m("LHE_Vpt") > 200 || mInt("LHE_Nb") == 0) return false;
        } else if (cursample->sampleNum == 122){
            if (m("LHE_Vpt") < 200 || mInt("LHE_Nb") == 0) return false;
        } else if (cursample->sampleNum == 141){
            if (m("LHE_Vpt") < 100 || m("LHE_Vpt") > 200 || nGenStatus2bHad == 0) return false;
        } else if (cursample->sampleNum == 142){
            if (m("LHE_Vpt") < 200 || nGenStatus2bHad == 0) return false;
        }
        if (cursample->sampleNum >= 150 && cursample->sampleNum<=156){
            if (m("LHE_Vpt") > 100) {
                if (mInt("LHE_Nb") != 0 || nGenStatus2bHad != 0) return false;
            }
       } else if (cursample->sampleNum == 160){
            if (m("LHE_Vpt") < 100 || m("LHE_Vpt") > 200 || mInt("LHE_Nb") == 0) return false;
        } else if (cursample->sampleNum == 161){
            if (m("LHE_Vpt") < 200 || mInt("LHE_Nb") == 0) return false;
        } else if (cursample->sampleNum == 162){
            if (m("LHE_Vpt") < 100 || m("LHE_Vpt") > 200 || nGenStatus2bHad == 0) return false;
        } else if (cursample->sampleNum == 163){
            if (m("LHE_Vpt") < 200 || nGenStatus2bHad == 0) return false;
        }
    }

    //if (cursample->sampleNum == 0) {
    //    //if (*f["json"] != 1 && !doCutFlowInPresel) return false;
    //    if (!doCutFlowInPresel) return false;
    //}

    // Trigger Selections
    if(m("reVType")){
        *in["Vtype"]= UpdatedVType();
    }

    //if (!PassVTypeAndTrigger(*in["Vtype"]) && !doCutFlowInPresel) {
    //    *in["controlSample"] = -1;
    //}


    if (mInt("Vtype")<2) {
        if (m("V_pt") < m("vptPreselCut") && !doCutFlowInPresel) { // preserve 50 < Vpt < 100 for the 2-lepton channels
            return false;
        }
    } else if (m("V_pt") < m("vptcut") && !doCutFlowInPresel) {    // else cut away Vpt < 100
        return false;
    }

    // Heppy jet corrections for JER/JEC are full correction, it's easier to just use the
    // correction on top of the nominal

    // FIXME configure b-tagged via config?
    //    if (int(*f["doCMVA"]) != 0) {
    //        // use CMVAv2 discriminator instead of CSV
    //        f["Jet_btagCSVV2"][i] = f["Jet_btagCMVA"][i];
    //    }

    if (m("smearJets")) {
        for (int i = 0; i < mInt("nJet"); i++) {
            float JERScale = m("JERScale"); // apply JER smearing x times the nominal smearing amount
            if (JERScale != 1.0) {
                smearJets(JERScale);
            }
        }
    }

    // FIXME why removed?
    //SetupFactorizedJECs(cursyst->name);

    // FIXME add Factorized JEC systs for MET
    //*f["MET_pt_JECUp_ratio"] = *f["met_shifted_JetEnUp_pt"] / *f["MET_pt"];
    //*f["MET_pt_JECDown_ratio"] = *f["met_shifted_JetEnDown_pt"] / *f["MET_pt"];

    // Preselect for two jets and one lepton which pass some minimum pt threshold
    int nPreselJets = 0;
    for (int i = 0; i < mInt("nJet"); i++) {
        if (m("Jet_PtReg",i) > m("JetPtPresel") && (mInt("Jet_puId",i) > 6 || m("Jet_Pt",i)>50) && mInt("Jet_lepFilter",i)) nPreselJets++;
    }



    atLeastOnePreselFatJet = false;
    if(int(m("doBoost"))==1){
        //take min pt
        float fatJetPtCut;
        fatJetPtCut       = std::min(m("fatJetPtCut_0lepchan")      ,m("fatJetPtCut_1lepchan"));
        fatJetPtCut       = std::min((double)fatJetPtCut            ,m("fatJetPtCut_2lepchan"));

        for(int iFatJet=0; iFatJet<mInt("nFatJet"); iFatJet++){
            if(m("FatJet_pt",iFatJet)>fatJetPtCut
                && fabs(m("FatJet_eta",iFatJet))<2.5){
                atLeastOnePreselFatJet=true;
                break;
            }
        }
    }

    if (nPreselJets < 2 && !atLeastOnePreselFatJet) return false;

    return true;
}

bool VHbbAnalysis::Analyze() {
    *in["sampleIndex"] = cursample->sampleNum;
    bool doCutFlow = int(m("doCutFlow")) != 0;
    *in["controlSample"] = 0;
    *in["cutFlow"] = 0;
    *b["twoResolvedJets"]=false;
    *b["oneMergedJet"]=false;
    *in["FatJetCand_index"]=-1;
    *f["MET_Pt_Nano"] = m("MET_Pt");
    
    // move fast requirements to the front
    if (debug>1000) std::cout << "Imposing json and trigger requirements" << std::endl;

    if (!PassVTypeAndTrigger(mInt("Vtype"))) {
        *in["controlSample"] = -1;
    }

    if (mInt("controlSample") > -1) {
        *in["cutFlow"] += 1; // json and trigger selection
    }

    bool passLeptonChannel=SelectLeptonChannel();
    
    if (passLeptonChannel && mInt("controlSample") > -1) {
        *in["cutFlow"] += 1; // lepton selection
    }

    // asking for doCutFlow instead of doOnlySignalRegion, as events with bad leptons will not go into CR
    if (!doCutFlow && mInt("controlSample") < 0) {
        return false;
    }

    bool selectJets=SelectJets();
    if(debug>1000) std::cout<<"selectJets "<<selectJets<<std::endl;

    if(!selectJets) return false;

    if (m("doOnlySignalRegion")>0 && mInt("controlSample") < 0) {
        return false;
    }

    if (mInt("controlSample") > -1) {
        *in["cutFlow"] += 1; // selected jets
    }
  
    // asking for doCutFlow instead of doOnlySignalRegion, as events with bad leptons will not go into CR
    if (!doCutFlow && mInt("controlSample") < 0) {
        return false;
    }


    bool recoedH=ReconstructHiggsCand();
    if(!recoedH) return false;

    bool recoedV=ReconstructVCand();
    if(!recoedV) return false;

    ComputeVHKinematics();

    ComputeOtherEventKinematics();

    ComputeBoostedVariables();
    
    ControlSampleSelection();
   
    //if HTXS_stage1_1_cat_pTjet30GeV isn't defined, set it to 0 for ShapeMaker
    if( m("sampleIndex")!=0){
        if( m("HTXS_stage1_1_cat_pTjet30GeV") != 0){
            *in["HTXS_stage1_1_cat_pTjet30GeV_defined"] = m("HTXS_stage1_1_cat_pTjet30GeV"); 
        }
        else{
            *in["HTXS_stage1_1_cat_pTjet30GeV_defined"] =0;
        }
    }

    
    if (doCutFlow && mInt("cutFlow") >= m("doCutFlow")) {
        // keep all preselected events for cutflow
        return true;
    } else {
        
        if( mInt("boostedControlSample")==-1 && mInt("controlSample") > -1) *in["boostedCategory"] = 0; //resolved
        if( mInt("boostedControlSample")!=-1 && mInt("controlSample") > -1) *in["boostedCategory"] = 1; //overlap
        if( mInt("boostedControlSample")!=-1 && mInt("controlSample") == -1) *in["boostedCategory"] = 2; //boosted
        
        if( (int(m("doBoost")) == 1) && mInt("boostedControlSample")!=-1) return true;
        return mInt("controlSample") > -1;
    }
}


void VHbbAnalysis::FinishEvent() {


// move the variable computation and BDT look up to top


//                                                                    _|
//        _|_|_|    _|_|    _|_|_|      _|_|    _|  _|_|    _|_|_|  _|_|_|_|    _|_|    _|  _|_|
//      _|    _|  _|_|_|_|  _|    _|  _|_|_|_|  _|_|      _|    _|    _|      _|    _|  _|_|
//      _|    _|  _|        _|    _|  _|        _|        _|    _|    _|      _|    _|  _|
//        _|_|_|    _|_|_|  _|    _|    _|_|_|  _|          _|_|_|      _|_|    _|_|    _|
//            _|
//        _|_|

    // Compare gen kinematics for b jets for signal vs. ttbar
    if (mInt("sampleIndex") != 0 && int(m("reRunGenInfo"))==1) {
        // Store number of gen b hadrons with status 2:
        *in["nGenStatus2bHad_re"] = 0;
        for(int indGP=0; indGP<mInt("nGenPart"); indGP++){
            if(mInt("GenPart_status",indGP)!=2) continue;
            if(((std::abs(mInt("GenPart_pdgId",indGP))/100)%10 ==5) || ((std::abs(mInt("GenPart_pdgId",indGP))/1000)%10==5)){
              *in["nGenStatus2bHad_re"]=mInt("nGenStatus2bHad_re")+1;
            }
        }

        //Only compute gen b hadrons in final state
        *in["nGenStatus2bHadFinal"]=0;
        std::vector<int> bMotherIds;
        for(int indGP=mInt("nGenPart")-1; indGP>-1; indGP--){//Loop starting at the end
            if(mInt("GenPart_status",indGP)==2
                && (((std::abs(mInt("GenPart_pdgId",indGP))/100)%10 ==5) || ((std::abs(mInt("GenPart_pdgId",indGP))/1000)%10==5)) ){
                if(std::find(bMotherIds.begin(),bMotherIds.end(),indGP) == bMotherIds.end()){//if particle is not a mother of another composite b particle
                    *in["nGenStatus2bHadFinal"]=mInt("nGenStatus2bHadFinal")+1;
                }
                bMotherIds.push_back(mInt("GenPart_genPartIdxMother",indGP));
            }
        }
        //Count gen b's in fatjet to split ttbar
        int nGenBJetsInFatJetCand = 0;
        for(int indGJ=0; indGJ<mInt("nGenJet"); indGJ++){
            if (mInt("GenJet_hadronFlavour",indGJ) == 5 && m("GenJet_pt",indGJ)>20 && fabs(m("GenJet_eta",indGJ)<2.4)) {//b GenJet criteria
               //Calculate dR GenJet, FatJetCand, see if it's less than 0.8
                double deltaR = AnalysisManager::EvalDeltaR( m("GenJet_eta",indGJ), m("GenJet_phi",indGJ), m("FatJetCand_eta"), m("FatJetCand_phi"));
                if(deltaR < 0.8) nGenBJetsInFatJetCand++;
            }
        }
        *in["nGenBJetsInFatJetCand"] = nGenBJetsInFatJetCand;
        
        //Gen AK8 Jets
        //Find closest Gen AK8 Jet to FatJetCand, and save its flavour
        double minDeltaR = 99.99;
        *in["GenFatJetCand_partonFlavour"]=-9999;//was f
        for(int indAK8Jet=0; indAK8Jet<mInt("nGenJetAK8"); indAK8Jet++){
            double deltaR = AnalysisManager::EvalDeltaR(m("GenJetAK8_eta",indAK8Jet), m("GenJetAK8_phi",indAK8Jet)  , m("FatJetCand_eta"), m("FatJetCand_phi"));
            if(deltaR < minDeltaR){
                minDeltaR = deltaR;
                *in["hGenFatJetCand"] = indAK8Jet;
                *in["GenFatJetCand_partonFlavour"]= mInt("GenJetAK8_partonFlavour",indAK8Jet);
            }
        }


        //first check for Higgs bosons:
        TLorentzVector GenBJ1, GenBJ2, GenBJJ;
        int mother_index=-1;
        int dau1_index=-1;
        int dau2_index=-1;
        int top_index_1=-1;
        int top_index_2=-1;
        *in["nGenTop_re"] = 0;
        *in["nGenVbosons_re"] = 0;
        *f["LeadGenVBoson_pt_re"] = 0.;
        *in["LeadGenVBoson_pdgId_re"] = 0;
        for(int indGP=0; indGP<mInt("nGenPart"); indGP++){
            //Check for Higgs boson and make sure it's the last copy -> 8192 = 2^13, 13th bit is IsLastCopy flag:
            if(fabs(mInt("GenPart_pdgId",indGP))==25 && (mInt("GenPart_statusFlags",indGP) & 8192)==8192 ) {
                mother_index=indGP;
            }
        }
        if(mother_index!=-1){
            for(int indGP=0; indGP<mInt("nGenPart"); indGP++){
                //Find b-quarks with a Higgs mother
                if(mInt("GenPart_genPartIdxMother",indGP)==mother_index && fabs(mInt("GenPart_pdgId",indGP))==5){
                    if(dau1_index>-1 && dau2_index>-1){
                        std::cout<<"This isn't supposed to happen!"<<std::endl;
                    }else if(dau1_index>-1){
                        dau2_index=indGP;
                    }else{
                        dau1_index=indGP;
                    }
                }
            }
        } else { //can also check for top quarks...
            for(int indGP=0; indGP<mInt("nGenPart"); indGP++){
                if(fabs(mInt("GenPart_pdgId",indGP))==6 && (mInt("GenPart_statusFlags",indGP) & 8192)==8192 ) {
                    *in["nGenTop_re"]=mInt("nGenTop_re")+1;
                    if(top_index_1>-1&&top_index_2>-1){
                        std::cout<<"This isn't supposed to happen!"<<std::endl;
                    } else if(top_index_1>-1){
                        top_index_2=indGP;
                    } else {
                        top_index_1=indGP;
                    }
                }
            }
            if(top_index_1!=-1){
                for(int indGP=0; indGP<mInt("nGenPart"); indGP++){
                    if(mInt("GenPart_genPartIdxMother",indGP)==top_index_1 && fabs(mInt("GenPart_pdgId",indGP))==5) dau1_index=indGP;
                }
            }
            if(top_index_2!=-1){
                for(int indGP=0; indGP<mInt("nGenPart"); indGP++){
                    if(mInt("GenPart_genPartIdxMother", indGP)==top_index_2 && fabs(mInt("GenPart_pdgId",indGP))==5){
                        if(dau1_index>-1) dau2_index=indGP;
                        else dau1_index=indGP;
                    }
                }
            }
        }
        if(dau1_index>-1){
            GenBJ1.SetPtEtaPhiM(m("GenPart_pt",dau1_index),m("GenPart_eta",dau1_index),m("GenPart_phi",dau1_index),4.2); //Mass only stored in nano if > 10 GeV
        }
        if(dau2_index>-1){
            GenBJ2.SetPtEtaPhiM(m("GenPart_pt",dau2_index),m("GenPart_eta",dau2_index),m("GenPart_phi",dau2_index),4.2);
        }

       *f["GenBJ1_pt_re"] = GenBJ1.Pt();
       *f["GenBJ1_eta_re"] = GenBJ1.Eta();
       *f["GenBJ1_phi_re"] = GenBJ1.Phi();
       *f["GenBJ1_mass_re"] = GenBJ1.M();
       *f["GenBJ2_pt_re"] = GenBJ2.Pt();
       *f["GenBJ2_eta_re"] = GenBJ2.Eta();
       *f["GenBJ2_phi_re"] = GenBJ2.Phi();
       *f["GenBJ2_mass_re"] = GenBJ2.M();

       GenBJJ = GenBJ1 + GenBJ2;
       *f["GenBJJ_pt_re"] = GenBJJ.Pt();
       *f["GenBJJ_eta_re"] = GenBJJ.Eta();
       *f["GenBJJ_phi_re"] = GenBJJ.Phi();
       *f["GenBJJ_mass_re"] = GenBJJ.M();
       *f["GenBJJ_dPhi_re"] = GenBJ2.DeltaPhi(GenBJ1);
       *f["GenBJJ_dR_re"] = GenBJ2.DeltaR(GenBJ1);
       *f["GenBJJ_dEta_re"] = fabs(GenBJ1.Eta() - GenBJ2.Eta());


       TLorentzVector GenLep1, GenLep2; // closest gen lep to either jet1 or jet2. Sometimes these could be the same lepton.
       double minDR1 = 999;
       double minDR2 = 999;
       double minSecDR2 = 999; //second closest gen lepton to second jet
       double part_mass=0.;
       int GenLepIndex1 = -1; // index of the lepton closest to jet 1
       int GenLepIndex2 = -1; // index of the lepton closest to jet 2
       int GenLepSecIndex2 = -1; // index of the lepton second closest to jet 2
       std::vector<int> VBosonIndices;
       VBosonIndices.clear();
       std::vector<int> LeptonIndices;
       LeptonIndices.clear();
       *ui["nW_re"]=0;
       *ui["nWlep_re"]=0;
       std::vector<unsigned int> WBosonIndices;
       std::vector<unsigned int> WBosonLeptonicIndices;
       WBosonIndices.clear();
       WBosonLeptonicIndices.clear();

       for (int i = 0; i< mInt("nGenPart"); i++){
           if((fabs(mInt("GenPart_pdgId",i))==24) && (mInt("GenPart_statusFlags",i) & 8192)==8192 ){
               VBosonIndices.push_back(i);
               WBosonIndices.push_back(i);
           }
           if((fabs(mInt("GenPart_pdgId",i))==23) && (mInt("GenPart_statusFlags",i) & 8192)==8192 ){
               VBosonIndices.push_back(i);
           }
       }
       *in["nGenVbosons_re"] = VBosonIndices.size();
       if(VBosonIndices.size()>0){
            for (int i = 0; i<mInt("nGenPart"); i++){
                if(fabs(mInt("GenPart_pdgId",i))>=11 && fabs(mInt("GenPart_pdgId",i))<=16 && mInt("GenPart_genPartIdxMother",i)>-1){
                    if(abs(mInt("GenPart_pdgId",mInt("GenPart_genPartIdxMother",i)))==24&& (mInt("GenPart_statusFlags",mInt("GenPart_genPartIdxMother",i)) & 8192)==8192 ){
                        if(std::find(WBosonLeptonicIndices.begin(), WBosonLeptonicIndices.end(), mInt("GenPart_genPartIdxMother",i)) == WBosonLeptonicIndices.end()){
                            WBosonLeptonicIndices.push_back(mInt("GenPart_genPartIdxMother",i));
                        }
                    }
                }
            }
       }
       //Check for W+Jet sample indices (40-99), DY->ellell sample indices (100-149) and Z->nunu sample indices (150-199)
       if(VBosonIndices.size()==0 && ( (mInt("sampleIndex")>=40 && mInt("sampleIndex")<=99) || (mInt("sampleIndex")>=100 && mInt("sampleIndex")<=149) || (mInt("sampleIndex")>=150 && mInt("sampleIndex")<=199) ) ){
           for (int i = 0; i<mInt("nGenPart"); i++){
               if((fabs(mInt("GenPart_pdgId",i))>=11 && fabs(mInt("GenPart_pdgId",i))<=16 && mInt("GenPart_status",i)==1 && (mInt("GenPart_statusFlags",i)&1)==1) || (fabs(mInt("GenPart_pdgId",i))==15 && (mInt("GenPart_statusFlags",i)&1)==1 && (mInt("GenPart_statusFlags",i)&2)==2) ){
                   LeptonIndices.push_back(i);
               }
           }
       }
       *ui["nW_re"]   =WBosonIndices.size();
       *ui["nWlep_re"]=WBosonLeptonicIndices.size();

       if (LeptonIndices.size()>1){
           std::sort(LeptonIndices.begin(),LeptonIndices.end(),[=](const int i1, const int i2){
                return f["GenPart_pt"][i1] > f["GenPart_pt"][i2];
           });
           TLorentzVector lep1,lep2;
           lep1.SetPtEtaPhiM(m("GenPart_pt",LeptonIndices.at(0)),m("GenPart_eta",LeptonIndices.at(0)),m("GenPart_phi",LeptonIndices.at(0)),m("GenPart_mass",LeptonIndices.at(0)));
           lep2.SetPtEtaPhiM(m("GenPart_pt",LeptonIndices.at(1)),m("GenPart_eta",LeptonIndices.at(1)),m("GenPart_phi",LeptonIndices.at(1)),m("GenPart_mass",LeptonIndices.at(1)));
           *f["LeadGenVBoson_pt_re"] = (lep1+lep2).Pt();
           *in["LeadGenVBoson_pdgId_re"] = 23;
           *in["nGenVbosons_re"] = 1;
       }
       if (VBosonIndices.size()>0){
           std::sort(VBosonIndices.begin(),VBosonIndices.end(),[=](const int i1, const int i2){
               return f["GenPart_pt"][i1] > f["GenPart_pt"][i2];
           });
           *f["LeadGenVBoson_pt_re"] = m("GenPart_pt",VBosonIndices.at(0));
           *in["LeadGenVBoson_pdgId_re"] = m("GenPart_pdgId",VBosonIndices.at(0));
           if(VBosonIndices.size()==1 && ( (mInt("sampleIndex")>=30 && mInt("sampleIndex")<=39) || (mInt("sampleIndex")>=340 && mInt("sampleIndex")<=350 ) ) ) *in["nGenVbosons_re"]=2; // Hack to make sure EWK/QCD reweighting doesn't get applied to diboson samples
       }
       for (int i = 0; i < mInt("nGenPart"); i++) {
           //Continue if not a final-state electron or muon which is either prompt or from a tau decay.
           if(!((fabs(mInt("GenPart_pdgId",i))==11||fabs(mInt("GenPart_pdgId",i))==13)&& mInt("GenPart_status",i)==1 && ((mInt("GenPart_statusFlags",i)&1)==1 || (mInt("GenPart_statusFlags",i)&32)==32  ))) continue;
           TLorentzVector gl;
           if(fabs(mInt("GenPart_pdgId",i))==11) part_mass=0.000511;
           if(fabs(mInt("GenPart_pdgId",i))==13) part_mass=0.105;
           gl.SetPtEtaPhiM(m("GenPart_pt",i), m("GenPart_eta",i), m("GenPart_phi",i),part_mass);
           double DR1 = gl.DeltaR(GenBJ1);
           double DR2 = gl.DeltaR(GenBJ2);

           if (DR1 <= minDR1) {
               minDR1 = DR1;
               GenLepIndex1 = i;
           }

           if (DR2 <= minDR2) {
               minSecDR2 = minDR2;
               GenLepSecIndex2 = GenLepIndex2;
               minDR2 = DR2;
               GenLepIndex2 = i;
           } else if (DR2 < minSecDR2) {
               minSecDR2 = DR2;
               GenLepSecIndex2 = i;
           }
       }

       if (GenLepIndex1 == GenLepIndex2) {
           // don't allow us to use the same lepton for each jet
            GenLepIndex2 = GenLepSecIndex2;
       }

       *in["GenLepIndex1_re"] = GenLepIndex1;
       *in["GenLepIndex2_re"] = GenLepIndex2;

        if (GenLepIndex1 != -1) {
            if(fabs(mInt("GenPart_pdgId",GenLepIndex1))==11) part_mass=0.000511;
            if(fabs(mInt("GenPart_pdgId",GenLepIndex1))==13) part_mass=0.105;
            GenLep1.SetPtEtaPhiM(m("GenPart_pt",GenLepIndex1), m("GenPart_eta",GenLepIndex1), m("GenPart_phi",GenLepIndex1),part_mass);
            *f["GenLep_GenBJ1_dR_re"] = GenLep1.DeltaR(GenBJ1);
            *f["GenLep_GenBJ1_dEta_re"] = fabs(GenLep1.Eta() - GenBJ1.Eta());
            *f["GenLep_GenBJ1_dPhi_re"] = GenLep1.DeltaPhi(GenBJ1);

            // try to reconstruct the top mass, although we've lost the neutrino so it will be shifted left
            TLorentzVector GenTop1 = GenLep1 + GenBJ1;
            *f["GenTop1_mass_re"] = GenTop1.M();
        } else {
            *f["GenLep_GenBJ1_dR_re"] = -999;
            *f["GenLep_GenBJ1_dEta_re"] = -999;
            *f["GenLep_GenBJ1_dPhi_re"] = -999;
            *f["GenTop1_mass_re"] = -999;
        }

        if (GenLepIndex2 != -1) {
            if(fabs(mInt("GenPart_pdgId",GenLepIndex2))==11) part_mass=0.000511;
            if(fabs(mInt("GenPart_pdgId",GenLepIndex2))==13) part_mass=0.105;
            GenLep2.SetPtEtaPhiM(m("GenPart_pt",GenLepIndex2), m("GenPart_eta",GenLepIndex2), m("GenPart_phi",GenLepIndex2),part_mass);
            *f["GenLep_GenBJ2_dR_re"] = GenLep2.DeltaR(GenBJ2);
            *f["GenLep_GenBJ2_dEta_re"] = (GenLep2.Eta(), GenBJ2.Eta());
            *f["GenLep_GenBJ2_dPhi_re"] = GenLep2.DeltaPhi(GenBJ2);

           // try to reconstruct the top mass, although we've lost the neutrino so it will be shifted left
            TLorentzVector GenTop2 = GenLep2 + GenBJ2;
            *f["GenTop2_mass_re"] = GenTop2.M();
        } else {
            *f["GenLep_GenBJ2_dR_re"] = -999;
            *f["GenLep_GenBJ2_dEta_re"] = -999;
            *f["GenLep_GenBJ2_dPhi_re"] = -999;
            *f["GenTop2_mass_re"] = -999;
        }


        // construct Gen W
        int w_index_1=-1;
        int w_index_2=-1;
        TLorentzVector GenW1, GenW2;
        for(int indGP=0; indGP<mInt("nGenPart"); indGP++){
            //Check for W boson and make sure it's the last copy -> 8192 = 2^13, 13th bit is IsLastCopy flag:
            if(fabs(mInt("GenPart_pdgId",indGP))==24 && (mInt("GenPart_statusFlags",indGP) & 8192)==8192 ) {
                if(w_index_1>-1&&w_index_2>-1){
                    std::cout<<"This isn't supposed to happen!"<<std::endl;
                } else if(w_index_1>-1){
                    w_index_2=indGP;
                } else w_index_1=indGP;
            }
        }
        if(w_index_1>-1){
            GenW1.SetPtEtaPhiM(m("GenPart_pt",w_index_1), m("GenPart_eta",w_index_1), m("GenPart_phi",w_index_1), m("GenPart_mass",w_index_1));
            *f["GenW_GenBJJ_dPhi_re"] = GenW1.DeltaPhi(GenBJJ);
            *f["GenW_GenBJJ_dEta_re"] = fabs(GenW1.Eta() - GenBJJ.Eta());
            // grab both W's in ttbar events
            if (w_index_2>-1){
                GenW2.SetPtEtaPhiM(m("GenPart_pt",w_index_2), m("GenPart_eta",w_index_2), m("GenPart_phi",w_index_2), m("GenPart_mass",w_index_2));
            }
        } else {
            *f["GenW_GenBJJ_dPhi_re"] = -999;
            *f["GenW_GenBJJ_dEta_re"] = -999;
        }

        std::vector<TLorentzVector> genWQuarks; // gen quarks from hadronic gen W decay
        if(w_index_1>-1||w_index_2>-1){
            for (int i = 0; i < mInt("nGenPart"); i++) {
                if( ( (w_index_1 > -1 && mInt("GenPart_genPartIdxMother",i)==w_index_1) ||(w_index_2>-1 && mInt("GenPart_genPartIdxMother",i)==w_index_2) ) && fabs(mInt("GenPart_pdgId",i))>=1&&fabs(mInt("GenPart_pdgId",i))<=6){
                    TLorentzVector v;
                    v.SetPtEtaPhiM(m("GenPart_pt",i), m("GenPart_eta",i), m("GenPart_phi",i), 0.); //FIXME could try to do something better here...
                    genWQuarks.push_back(v);
                }
            }
        }

        //int nSelectedJetsMatched = 0; // count the number (0, 1, 2) of selected jets matched to the real bottom quarks
        // Match Jets with Gen B Jets from Higgs/Tops
        for (int i = 0; i < mInt("nJet"); i++) {
            in["Jet_genJetMatchId_re"][i] = 0; // 0 if no gen match, 1 for pt-leading b-jet, 2 for pt sub-leading b-jet, 3 if matched to jet from hadronic W decay
            TLorentzVector Jet;
            Jet.SetPtEtaPhiM(m("Jet_bReg",i), m("Jet_eta",i), m("Jet_phi",i), m("Jet_mass",i) * (m("Jet_bReg",i) / m("Jet_Pt",i)));

            //double dR1 = Jet.DeltaR(GenHJ1);
            //double dR2 = Jet.DeltaR(GenHJ2);
            double dR1 = 999;
            if (GenBJ1.Pt() > 0) dR1 = Jet.DeltaR(GenBJ1);
            double dR2 = 999;
            if (GenBJ2.Pt() > 0) dR2 = Jet.DeltaR(GenBJ2);

            // try to match the jet to one of the jets from hadronic W decay
            double dR3 = 999;
            for (int j = 0; j < (int) genWQuarks.size(); j++) {
                double Jet_genWQuarkDR = Jet.DeltaR(genWQuarks[j]);
                if (Jet_genWQuarkDR < dR3) {
                    dR3 = Jet_genWQuarkDR;
                }
            }

            f["Jet_genWQuarkDR_re"][i] = dR3;
            if (dR3 < std::min(dR1, dR2) && dR3 < 0.5) {
                in["Jet_genJetMatchId_re"][i] = 3;
            } else if (dR1 <= dR2 && dR1 < 0.5) {
                in["Jet_genJetMatchId_re"][i] = 1;
            } else if (dR2 < 0.5) {
                in["Jet_genJetMatchId_re"][i] = 2;
            }

            f["Jet_genHJetMinDR_re"][i] = std::min(dR1, dR2);

            if (i == mInt("hJetInd1")) {
                *f["hJet1_matchedMinDR_re"] = m("Jet_genHJetMinDR_re",i);
            } else if (i == mInt("hJetInd2")) {
                *f["hJet2_matchedMinDR_re"] = m("Jet_genHJetMinDR_re",i);
            }
        }
    }
    std::vector<int> ptBoundL = {20, 30, 40, 60, 100};
    std::vector<int> ptBoundH = {30, 40, 60, 100, 1000};
    std::vector<double> etaBoundL = {0.0,0.8,1.6};
    std::vector<double> etaBoundH = {0.8,1.6,2.4};

    std::map<std::string,float> bTagWeightSys;

    *f["bTagWeight"] = 1.0;

    // b tag weights
    *f["bTagWeight_JESUp"] = 1.0;
    *f["bTagWeight_JESDown"] = 1.0;
    *f["bTagWeight_HFUp"] = 1.0;
    *f["bTagWeight_HFDown"] = 1.0;
    *f["bTagWeight_LFUp"] = 1.0;
    *f["bTagWeight_LFDown"] = 1.0;
    *f["bTagWeight_LFStats1Up"] = 1.0;
    *f["bTagWeight_LFStats1Down"] = 1.0;
    *f["bTagWeight_LFStats2Up"] = 1.0;
    *f["bTagWeight_LFStats2Down"] = 1.0;
    *f["bTagWeight_HFStats1Up"] = 1.0;
    *f["bTagWeight_HFStats1Down"] = 1.0;
    *f["bTagWeight_HFStats2Up"] = 1.0;
    *f["bTagWeight_HFStats2Down"] = 1.0;
    *f["bTagWeight_cErr1Up"] = 1.0;
    *f["bTagWeight_cErr1Down"] = 1.0;
    *f["bTagWeight_cErr2Up"] = 1.0;
    *f["bTagWeight_cErr2Down"] = 1.0;
    for (int i = 0; i < (int) ptBoundL.size(); i++) {
        for (int j = 0; j < (int) etaBoundL.size(); j++) {

             // b tag weights
             *f[Form("bTagWeight_JESUp_pt%i_eta%i",i,j)] = 1.0;
             *f[Form("bTagWeight_JESDown_pt%i_eta%i",i,j)] = 1.0;
             *f[Form("bTagWeight_HFUp_pt%i_eta%i",i,j)] = 1.0;
             *f[Form("bTagWeight_HFDown_pt%i_eta%i",i,j)] = 1.0;
             *f[Form("bTagWeight_LFUp_pt%i_eta%i",i,j)] = 1.0;
             *f[Form("bTagWeight_LFDown_pt%i_eta%i",i,j)] = 1.0;
             *f[Form("bTagWeight_LFStats1Up_pt%i_eta%i",i,j)] = 1.0;
             *f[Form("bTagWeight_LFStats1Down_pt%i_eta%i",i,j)] = 1.0;
             *f[Form("bTagWeight_LFStats2Up_pt%i_eta%i",i,j)] = 1.0;
             *f[Form("bTagWeight_LFStats2Down_pt%i_eta%i",i,j)] = 1.0;
             *f[Form("bTagWeight_HFStats1Up_pt%i_eta%i",i,j)] = 1.0;
             *f[Form("bTagWeight_HFStats1Down_pt%i_eta%i",i,j)] = 1.0;
             *f[Form("bTagWeight_HFStats2Up_pt%i_eta%i",i,j)] = 1.0;
             *f[Form("bTagWeight_HFStats2Down_pt%i_eta%i",i,j)] = 1.0;
             *f[Form("bTagWeight_cErr1Up_pt%i_eta%i",i,j)] = 1.0;
             *f[Form("bTagWeight_cErr1Down_pt%i_eta%i",i,j)] = 1.0;
             *f[Form("bTagWeight_cErr2Up_pt%i_eta%i",i,j)] = 1.0;
             *f[Form("bTagWeight_cErr2Down_pt%i_eta%i",i,j)] = 1.0;

             bTagWeightSys[Form("bTagWeight_JESUp_pt%i_eta%i",i,j)] = 1.0;
             bTagWeightSys[Form("bTagWeight_JESDown_pt%i_eta%i",i,j)] = 1.0;
             bTagWeightSys[Form("bTagWeight_HFUp_pt%i_eta%i",i,j)] = 1.0;
             bTagWeightSys[Form("bTagWeight_HFDown_pt%i_eta%i",i,j)] = 1.0;
             bTagWeightSys[Form("bTagWeight_LFUp_pt%i_eta%i",i,j)] = 1.0;
             bTagWeightSys[Form("bTagWeight_LFDown_pt%i_eta%i",i,j)] = 1.0;
             bTagWeightSys[Form("bTagWeight_LFStats1Up_pt%i_eta%i",i,j)] = 1.0;
             bTagWeightSys[Form("bTagWeight_LFStats1Down_pt%i_eta%i",i,j)] = 1.0;
             bTagWeightSys[Form("bTagWeight_LFStats2Up_pt%i_eta%i",i,j)] = 1.0;
             bTagWeightSys[Form("bTagWeight_LFStats2Down_pt%i_eta%i",i,j)] = 1.0;
             bTagWeightSys[Form("bTagWeight_HFStats1Up_pt%i_eta%i",i,j)] = 1.0;
             bTagWeightSys[Form("bTagWeight_HFStats1Down_pt%i_eta%i",i,j)] = 1.0;
             bTagWeightSys[Form("bTagWeight_HFStats2Up_pt%i_eta%i",i,j)] = 1.0;
             bTagWeightSys[Form("bTagWeight_HFStats2Down_pt%i_eta%i",i,j)] = 1.0;
             bTagWeightSys[Form("bTagWeight_cErr1Up_pt%i_eta%i",i,j)] = 1.0;
             bTagWeightSys[Form("bTagWeight_cErr1Down_pt%i_eta%i",i,j)] = 1.0;
             bTagWeightSys[Form("bTagWeight_cErr2Up_pt%i_eta%i",i,j)] = 1.0;
             bTagWeightSys[Form("bTagWeight_cErr2Down_pt%i_eta%i",i,j)] = 1.0;
        }
    }

    if (bTagCalibReader) {
        if (debug>101) std::cout<<"evaluating btag scale factors"<<std::endl;

        // TODO uncertainties
        float bTagWeight = 1.;

        std::string taggerForEvaluation;
        if(m("dataYear")==2016){
            taggerForEvaluation.append("Jet_btagCMVA");
        } else {
            taggerForEvaluation=taggerName;
        }


        for (int i=0; i<mInt("nJet"); i++)
        {
            if ( (m("Jet_puId", i) > 6 || m("Jet_Pt",i)>50)
                && m("Jet_lepFilter", i) > 0
                && m("Jet_Pt", i)>m("JetPtPresel")
                && fabs(m("Jet_eta", i))<=m("JetEtaCut"))
            {
                int hadron_flav = mInt("Jet_hadronFlavour", i);

                auto flav = (hadron_flav==5) ? BTagEntry::FLAV_B :
                            (hadron_flav==4) ? BTagEntry::FLAV_C :
                            BTagEntry::FLAV_UDSG;
                float bTagWeight_jet = bTagCalibReader->eval_auto_bounds(
                    "central",
                    flav,
                    fabs(m("Jet_eta", i)),
                    m("Jet_Pt", i),
                    m(taggerForEvaluation, i)
                );
                bTagWeight *= bTagWeight_jet;
               float bTagWeight_jes_up = 1.;
               float bTagWeight_jes_down = 1.;
               float bTagWeight_hf_up = 1.;
               float bTagWeight_hf_down = 1.;
               float bTagWeight_lf_up = 1.;
               float bTagWeight_lf_down = 1.;
               float bTagWeight_lfstats1_up = 1.;
               float bTagWeight_lfstats1_down = 1.;
               float bTagWeight_lfstats2_up = 1.;
               float bTagWeight_lfstats2_down = 1.;
               float bTagWeight_hfstats1_up = 1.;
               float bTagWeight_hfstats1_down = 1.;
               float bTagWeight_hfstats2_up = 1.;
               float bTagWeight_hfstats2_down = 1.;
               float bTagWeight_cferr1_up = 1.;
               float bTagWeight_cferr1_down = 1.;
               float bTagWeight_cferr2_up = 1.;
               float bTagWeight_cferr2_down = 1.;
               if(hadron_flav==5){
                   bTagWeight_jes_up = bTagCalibReader->eval_auto_bounds(
                       "up_jes", flav,fabs(m("Jet_eta", i)),m("Jet_Pt", i),m(taggerForEvaluation, i));
                   bTagWeight_jes_down = bTagCalibReader->eval_auto_bounds(
                       "down_jes", flav,fabs(m("Jet_eta", i)),m("Jet_Pt", i),m(taggerForEvaluation, i));
                   bTagWeight_lf_up = bTagCalibReader->eval_auto_bounds(
                       "up_lf", flav,fabs(m("Jet_eta", i)),m("Jet_Pt", i),m(taggerForEvaluation, i));
                   bTagWeight_lf_down = bTagCalibReader->eval_auto_bounds(
                       "down_lf", flav,fabs(m("Jet_eta", i)),m("Jet_Pt", i),m(taggerForEvaluation, i));
                   bTagWeight_hfstats1_up = bTagCalibReader->eval_auto_bounds(
                       "up_hfstats1", flav,fabs(m("Jet_eta", i)),m("Jet_Pt", i),m(taggerForEvaluation, i));
                   bTagWeight_hfstats1_down = bTagCalibReader->eval_auto_bounds(
                       "down_hfstats1", flav,fabs(m("Jet_eta", i)),m("Jet_Pt", i),m(taggerForEvaluation, i));
                   bTagWeight_hfstats2_up = bTagCalibReader->eval_auto_bounds(
                       "up_hfstats2", flav,fabs(m("Jet_eta", i)),m("Jet_Pt", i),m(taggerForEvaluation, i));
                   bTagWeight_hfstats2_down = bTagCalibReader->eval_auto_bounds(
                       "down_hfstats2", flav,fabs(m("Jet_eta", i)),m("Jet_Pt", i),m(taggerForEvaluation, i));
                   bTagWeight_cferr1_up = bTagCalibReader->eval_auto_bounds(
                       "central", flav,fabs(m("Jet_eta", i)),m("Jet_Pt", i),m(taggerForEvaluation, i));
                   bTagWeight_cferr1_down = bTagCalibReader->eval_auto_bounds(
                       "central", flav,fabs(m("Jet_eta", i)),m("Jet_Pt", i),m(taggerForEvaluation, i));
                   bTagWeight_cferr2_up = bTagCalibReader->eval_auto_bounds(
                       "central", flav,fabs(m("Jet_eta", i)),m("Jet_Pt", i),m(taggerForEvaluation, i));
                   bTagWeight_cferr2_down = bTagCalibReader->eval_auto_bounds(
                       "central", flav,fabs(m("Jet_eta", i)),m("Jet_Pt", i),m(taggerForEvaluation, i));
                   bTagWeight_hf_up = bTagCalibReader->eval_auto_bounds(
                       "central", flav,fabs(m("Jet_eta", i)),m("Jet_Pt", i),m(taggerForEvaluation, i));
                   bTagWeight_hf_down = bTagCalibReader->eval_auto_bounds(
                       "central", flav,fabs(m("Jet_eta", i)),m("Jet_Pt", i),m(taggerForEvaluation, i));
                   bTagWeight_lfstats1_up = bTagCalibReader->eval_auto_bounds(
                       "central", flav,fabs(m("Jet_eta", i)),m("Jet_Pt", i),m(taggerForEvaluation, i));
                   bTagWeight_lfstats1_down = bTagCalibReader->eval_auto_bounds(
                       "central", flav,fabs(m("Jet_eta", i)),m("Jet_Pt", i),m(taggerForEvaluation, i));
                   bTagWeight_lfstats2_up = bTagCalibReader->eval_auto_bounds(
                       "central", flav,fabs(m("Jet_eta", i)),m("Jet_Pt", i),m(taggerForEvaluation, i));
                   bTagWeight_lfstats2_down = bTagCalibReader->eval_auto_bounds(
                       "central", flav,fabs(m("Jet_eta", i)),m("Jet_Pt", i),m(taggerForEvaluation, i));
               } else if (hadron_flav==4){
                   bTagWeight_jes_up = bTagCalibReader->eval_auto_bounds(
                       "central", flav,fabs(m("Jet_eta", i)),m("Jet_Pt", i),m(taggerForEvaluation, i));
                   bTagWeight_jes_down = bTagCalibReader->eval_auto_bounds(
                       "central", flav,fabs(m("Jet_eta", i)),m("Jet_Pt", i),m(taggerForEvaluation, i));
                   bTagWeight_lf_up = bTagCalibReader->eval_auto_bounds(
                       "central", flav,fabs(m("Jet_eta", i)),m("Jet_Pt", i),m(taggerForEvaluation, i));
                   bTagWeight_lf_down = bTagCalibReader->eval_auto_bounds(
                       "central", flav,fabs(m("Jet_eta", i)),m("Jet_Pt", i),m(taggerForEvaluation, i));
                   bTagWeight_hfstats1_up = bTagCalibReader->eval_auto_bounds(
                       "central", flav,fabs(m("Jet_eta", i)),m("Jet_Pt", i),m(taggerForEvaluation, i));
                   bTagWeight_hfstats1_down = bTagCalibReader->eval_auto_bounds(
                       "central", flav,fabs(m("Jet_eta", i)),m("Jet_Pt", i),m(taggerForEvaluation, i));
                   bTagWeight_hfstats2_up = bTagCalibReader->eval_auto_bounds(
                       "central", flav,fabs(m("Jet_eta", i)),m("Jet_Pt", i),m(taggerForEvaluation, i));
                   bTagWeight_hfstats2_down = bTagCalibReader->eval_auto_bounds(
                       "central", flav,fabs(m("Jet_eta", i)),m("Jet_Pt", i),m(taggerForEvaluation, i));
                   bTagWeight_cferr1_up = bTagCalibReader->eval_auto_bounds(
                       "up_cferr1", flav,fabs(m("Jet_eta", i)),m("Jet_Pt", i),m(taggerForEvaluation, i));
                   bTagWeight_cferr1_down = bTagCalibReader->eval_auto_bounds(
                       "down_cferr1", flav,fabs(m("Jet_eta", i)),m("Jet_Pt", i),m(taggerForEvaluation, i));
                   bTagWeight_cferr2_up = bTagCalibReader->eval_auto_bounds(
                       "up_cferr2", flav,fabs(m("Jet_eta", i)),m("Jet_Pt", i),m(taggerForEvaluation, i));
                   bTagWeight_cferr2_down = bTagCalibReader->eval_auto_bounds(
                       "down_cferr2", flav,fabs(m("Jet_eta", i)),m("Jet_Pt", i),m(taggerForEvaluation, i));
                   bTagWeight_hf_up = bTagCalibReader->eval_auto_bounds(
                       "central", flav,fabs(m("Jet_eta", i)),m("Jet_Pt", i),m(taggerForEvaluation, i));
                   bTagWeight_hf_down = bTagCalibReader->eval_auto_bounds(
                       "central", flav,fabs(m("Jet_eta", i)),m("Jet_Pt", i),m(taggerForEvaluation, i));
                   bTagWeight_lfstats1_up = bTagCalibReader->eval_auto_bounds(
                       "central", flav,fabs(m("Jet_eta", i)),m("Jet_Pt", i),m(taggerForEvaluation, i));
                   bTagWeight_lfstats1_down = bTagCalibReader->eval_auto_bounds(
                       "central", flav,fabs(m("Jet_eta", i)),m("Jet_Pt", i),m(taggerForEvaluation, i));
                   bTagWeight_lfstats2_up = bTagCalibReader->eval_auto_bounds(
                       "central", flav,fabs(m("Jet_eta", i)),m("Jet_Pt", i),m(taggerForEvaluation, i));
                   bTagWeight_lfstats2_down = bTagCalibReader->eval_auto_bounds(
                       "central", flav,fabs(m("Jet_eta", i)),m("Jet_Pt", i),m(taggerForEvaluation, i));
               } else {
                   bTagWeight_jes_up = bTagCalibReader->eval_auto_bounds(
                       "up_jes", flav,fabs(m("Jet_eta", i)),m("Jet_Pt", i),m(taggerForEvaluation, i));
                   bTagWeight_jes_down = bTagCalibReader->eval_auto_bounds(
                       "down_jes", flav,fabs(m("Jet_eta", i)),m("Jet_Pt", i),m(taggerForEvaluation, i));
                   bTagWeight_hf_up = bTagCalibReader->eval_auto_bounds(
                       "up_hf", flav,fabs(m("Jet_eta", i)),m("Jet_Pt", i),m(taggerForEvaluation, i));
                   bTagWeight_hf_down = bTagCalibReader->eval_auto_bounds(
                       "down_hf", flav,fabs(m("Jet_eta", i)),m("Jet_Pt", i),m(taggerForEvaluation, i));
                   bTagWeight_lfstats1_up = bTagCalibReader->eval_auto_bounds(
                       "up_lfstats1", flav,fabs(m("Jet_eta", i)),m("Jet_Pt", i),m(taggerForEvaluation, i));
                   bTagWeight_lfstats1_down = bTagCalibReader->eval_auto_bounds(
                       "down_lfstats1", flav,fabs(m("Jet_eta", i)),m("Jet_Pt", i),m(taggerForEvaluation, i));
                   bTagWeight_lfstats2_up = bTagCalibReader->eval_auto_bounds(
                       "up_lfstats2", flav,fabs(m("Jet_eta", i)),m("Jet_Pt", i),m(taggerForEvaluation, i));
                   bTagWeight_lfstats2_down = bTagCalibReader->eval_auto_bounds(
                       "down_lfstats2", flav,fabs(m("Jet_eta", i)),m("Jet_Pt", i),m(taggerForEvaluation, i));
                   bTagWeight_lf_up = bTagCalibReader->eval_auto_bounds(
                       "central", flav,fabs(m("Jet_eta", i)),m("Jet_Pt", i),m(taggerForEvaluation, i));
                   bTagWeight_lf_down = bTagCalibReader->eval_auto_bounds(
                       "central", flav,fabs(m("Jet_eta", i)),m("Jet_Pt", i),m(taggerForEvaluation, i));
                   bTagWeight_hfstats1_up = bTagCalibReader->eval_auto_bounds(
                       "central", flav,fabs(m("Jet_eta", i)),m("Jet_Pt", i),m(taggerForEvaluation, i));
                   bTagWeight_hfstats1_down = bTagCalibReader->eval_auto_bounds(
                       "central", flav,fabs(m("Jet_eta", i)),m("Jet_Pt", i),m(taggerForEvaluation, i));
                   bTagWeight_hfstats2_up = bTagCalibReader->eval_auto_bounds(
                       "central", flav,fabs(m("Jet_eta", i)),m("Jet_Pt", i),m(taggerForEvaluation, i));
                   bTagWeight_hfstats2_down = bTagCalibReader->eval_auto_bounds(
                       "central", flav,fabs(m("Jet_eta", i)),m("Jet_Pt", i),m(taggerForEvaluation, i));
                   bTagWeight_cferr1_up = bTagCalibReader->eval_auto_bounds(
                       "central", flav,fabs(m("Jet_eta", i)),m("Jet_Pt", i),m(taggerForEvaluation, i));
                   bTagWeight_cferr1_down = bTagCalibReader->eval_auto_bounds(
                       "central", flav,fabs(m("Jet_eta", i)),m("Jet_Pt", i),m(taggerForEvaluation, i));
                   bTagWeight_cferr2_up = bTagCalibReader->eval_auto_bounds(
                       "central", flav,fabs(m("Jet_eta", i)),m("Jet_Pt", i),m(taggerForEvaluation, i));
                   bTagWeight_cferr2_down = bTagCalibReader->eval_auto_bounds(
                       "central", flav,fabs(m("Jet_eta", i)),m("Jet_Pt", i),m(taggerForEvaluation, i));
               }
                *f["bTagWeight_JESUp"]         *= bTagWeight_jes_up;
                *f["bTagWeight_JESDown"]       *= bTagWeight_jes_down;
                *f["bTagWeight_HFUp"]          *= bTagWeight_hf_up;
                *f["bTagWeight_HFDown"]        *= bTagWeight_hf_down;
                *f["bTagWeight_LFUp"]          *= bTagWeight_lf_up;
                *f["bTagWeight_LFDown"]        *= bTagWeight_lf_down;
                *f["bTagWeight_LFStats1Up"]    *= bTagWeight_lfstats1_up;
                *f["bTagWeight_LFStats1Down"]  *= bTagWeight_lfstats1_down;
                *f["bTagWeight_LFStats2Up"]    *= bTagWeight_lfstats2_up;
                *f["bTagWeight_LFStats2Down"]  *= bTagWeight_lfstats2_down;
                *f["bTagWeight_HFStats1Up"]    *= bTagWeight_hfstats1_up;
                *f["bTagWeight_HFStats1Down"]  *= bTagWeight_hfstats1_down;
                *f["bTagWeight_HFStats2Up"]    *= bTagWeight_hfstats2_up;
                *f["bTagWeight_HFStats2Down"]  *= bTagWeight_hfstats2_down;
                *f["bTagWeight_cErr1Up"]       *= bTagWeight_cferr1_up;
                *f["bTagWeight_cErr1Down"]     *= bTagWeight_cferr1_down;
                *f["bTagWeight_cErr2Up"]       *= bTagWeight_cferr2_up;
                *f["bTagWeight_cErr2Down"]     *= bTagWeight_cferr2_down;
                for (int j = 0; j < (int) ptBoundL.size(); j++) {
                    for (int k = 0; k < (int) etaBoundL.size(); k++) {
                       if (m("Jet_Pt", i) < ptBoundH[j] && m("Jet_Pt", i) >= ptBoundL[j] && fabs(m("Jet_eta", i)) > etaBoundL[k] && fabs(m("Jet_eta", i)) <= etaBoundH[k]) {
                           bTagWeightSys[Form("bTagWeight_JESUp_pt%i_eta%i",j,k)] *= bTagWeight_jes_up;
                           bTagWeightSys[Form("bTagWeight_JESDown_pt%i_eta%i",j,k)] *= bTagWeight_jes_down;
                           bTagWeightSys[Form("bTagWeight_HFUp_pt%i_eta%i",j,k)]  *= bTagWeight_hf_up;
                           bTagWeightSys[Form("bTagWeight_HFDown_pt%i_eta%i",j,k)] *= bTagWeight_hf_down;
                           bTagWeightSys[Form("bTagWeight_LFUp_pt%i_eta%i",j,k)] *= bTagWeight_lf_up;
                           bTagWeightSys[Form("bTagWeight_LFDown_pt%i_eta%i",j,k)] *= bTagWeight_lf_down;
                           bTagWeightSys[Form("bTagWeight_LFStats1Up_pt%i_eta%i",j,k)] *= bTagWeight_lfstats1_up;
                           bTagWeightSys[Form("bTagWeight_LFStats1Down_pt%i_eta%i",j,k)] *= bTagWeight_lfstats1_down;
                           bTagWeightSys[Form("bTagWeight_LFStats2Up_pt%i_eta%i",j,k)] *= bTagWeight_lfstats2_up;
                           bTagWeightSys[Form("bTagWeight_LFStats2Down_pt%i_eta%i",j,k)] *= bTagWeight_lfstats2_down;
                           bTagWeightSys[Form("bTagWeight_HFStats1Up_pt%i_eta%i",j,k)] *= bTagWeight_hfstats1_up;
                           bTagWeightSys[Form("bTagWeight_HFStats1Down_pt%i_eta%i",j,k)] *= bTagWeight_hfstats1_down;
                           bTagWeightSys[Form("bTagWeight_HFStats2Up_pt%i_eta%i",j,k)] *= bTagWeight_hfstats2_up;
                           bTagWeightSys[Form("bTagWeight_HFStats2Down_pt%i_eta%i",j,k)] *= bTagWeight_hfstats2_down;
                           bTagWeightSys[Form("bTagWeight_cErr1Up_pt%i_eta%i",j,k)] *= bTagWeight_cferr1_up;
                           bTagWeightSys[Form("bTagWeight_cErr1Down_pt%i_eta%i",j,k)] *= bTagWeight_cferr1_down;
                           bTagWeightSys[Form("bTagWeight_cErr2Up_pt%i_eta%i",j,k)] *= bTagWeight_cferr2_up;
                           bTagWeightSys[Form("bTagWeight_cErr2Down_pt%i_eta%i",j,k)] *= bTagWeight_cferr2_down;
                        }
                        else {
                           bTagWeightSys[Form("bTagWeight_JESUp_pt%i_eta%i",j,k)] *= bTagWeight_jet;
                           bTagWeightSys[Form("bTagWeight_JESDown_pt%i_eta%i",j,k)] *= bTagWeight_jet;
                           bTagWeightSys[Form("bTagWeight_HFUp_pt%i_eta%i",j,k)]  *= bTagWeight_jet;
                           bTagWeightSys[Form("bTagWeight_HFDown_pt%i_eta%i",j,k)] *= bTagWeight_jet;
                           bTagWeightSys[Form("bTagWeight_LFUp_pt%i_eta%i",j,k)] *= bTagWeight_jet;
                           bTagWeightSys[Form("bTagWeight_LFDown_pt%i_eta%i",j,k)] *= bTagWeight_jet;
                           bTagWeightSys[Form("bTagWeight_LFStats1Up_pt%i_eta%i",j,k)] *= bTagWeight_jet;
                           bTagWeightSys[Form("bTagWeight_LFStats1Down_pt%i_eta%i",j,k)] *= bTagWeight_jet;
                           bTagWeightSys[Form("bTagWeight_LFStats2Up_pt%i_eta%i",j,k)] *= bTagWeight_jet;
                           bTagWeightSys[Form("bTagWeight_LFStats2Down_pt%i_eta%i",j,k)] *= bTagWeight_jet;
                           bTagWeightSys[Form("bTagWeight_HFStats1Up_pt%i_eta%i",j,k)] *= bTagWeight_jet;
                           bTagWeightSys[Form("bTagWeight_HFStats1Down_pt%i_eta%i",j,k)] *= bTagWeight_jet;
                           bTagWeightSys[Form("bTagWeight_HFStats2Up_pt%i_eta%i",j,k)] *= bTagWeight_jet;
                           bTagWeightSys[Form("bTagWeight_HFStats2Down_pt%i_eta%i",j,k)] *= bTagWeight_jet;
                           bTagWeightSys[Form("bTagWeight_cErr1Up_pt%i_eta%i",j,k)] *= bTagWeight_jet;
                           bTagWeightSys[Form("bTagWeight_cErr1Down_pt%i_eta%i",j,k)] *= bTagWeight_jet;
                           bTagWeightSys[Form("bTagWeight_cErr2Up_pt%i_eta%i",j,k)] *= bTagWeight_jet;
                           bTagWeightSys[Form("bTagWeight_cErr2Down_pt%i_eta%i",j,k)] *= bTagWeight_jet;
                        }
                    }
                }
            }
        }
        *f["bTagWeight"] = bTagWeight;
        //if(m("dataYear")==2016){
        //    *f["bTagWeight"] = 1.0; //stealing CMVA systematics but not central value
       // }
        *f["bTagWeight_JESUp"]         /= bTagWeight;
        *f["bTagWeight_JESDown"]       /= bTagWeight;
        *f["bTagWeight_HFUp"]          /= bTagWeight;
        *f["bTagWeight_HFDown"]        /= bTagWeight;
        *f["bTagWeight_LFUp"]          /= bTagWeight;
        *f["bTagWeight_LFDown"]        /= bTagWeight;
        *f["bTagWeight_LFStats1Up"]    /= bTagWeight;
        *f["bTagWeight_LFStats1Down"]  /= bTagWeight;
        *f["bTagWeight_LFStats2Up"]    /= bTagWeight;
        *f["bTagWeight_LFStats2Down"]  /= bTagWeight;
        *f["bTagWeight_HFStats1Up"]    /= bTagWeight;
        *f["bTagWeight_HFStats1Down"]  /= bTagWeight;
        *f["bTagWeight_HFStats2Up"]    /= bTagWeight;
        *f["bTagWeight_HFStats2Down"]  /= bTagWeight;
        *f["bTagWeight_cErr1Up"]       /= bTagWeight;
        *f["bTagWeight_cErr1Down"]     /= bTagWeight;
        *f["bTagWeight_cErr2Up"]       /= bTagWeight;
        *f["bTagWeight_cErr2Down"]     /= bTagWeight;
        for (int i = 0; i < (int) ptBoundL.size(); i++) {
            for (int j = 0; j < (int) etaBoundL.size(); j++) {
                *f[Form("bTagWeight_JESUp_pt%i_eta%i",i,j)] = bTagWeightSys[Form("bTagWeight_JESUp_pt%i_eta%i",i,j)]/bTagWeight;
                *f[Form("bTagWeight_JESDown_pt%i_eta%i",i,j)] = bTagWeightSys[Form("bTagWeight_JESDown_pt%i_eta%i",i,j)]/bTagWeight;
                *f[Form("bTagWeight_HFUp_pt%i_eta%i",i,j)] = bTagWeightSys[Form("bTagWeight_HFUp_pt%i_eta%i",i,j)] /bTagWeight;
                *f[Form("bTagWeight_HFDown_pt%i_eta%i",i,j)] = bTagWeightSys[Form("bTagWeight_HFDown_pt%i_eta%i",i,j)]/bTagWeight;
                *f[Form("bTagWeight_LFUp_pt%i_eta%i",i,j)] = bTagWeightSys[Form("bTagWeight_LFUp_pt%i_eta%i",i,j)]/bTagWeight;
                *f[Form("bTagWeight_LFDown_pt%i_eta%i",i,j)] = bTagWeightSys[Form("bTagWeight_LFDown_pt%i_eta%i",i,j)]/bTagWeight;
                *f[Form("bTagWeight_LFStats1Up_pt%i_eta%i",i,j)] = bTagWeightSys[Form("bTagWeight_LFStats1Up_pt%i_eta%i",i,j)]/bTagWeight;
                *f[Form("bTagWeight_LFStats1Down_pt%i_eta%i",i,j)] = bTagWeightSys[Form("bTagWeight_LFStats1Down_pt%i_eta%i",i,j)] /bTagWeight;
                *f[Form("bTagWeight_LFStats2Up_pt%i_eta%i",i,j)] =  bTagWeightSys[Form("bTagWeight_LFStats2Up_pt%i_eta%i",i,j)]/bTagWeight;
                *f[Form("bTagWeight_LFStats2Down_pt%i_eta%i",i,j)] = bTagWeightSys[Form("bTagWeight_LFStats2Down_pt%i_eta%i",i,j)] /bTagWeight;
                *f[Form("bTagWeight_HFStats1Up_pt%i_eta%i",i,j)] = bTagWeightSys[Form("bTagWeight_HFStats1Up_pt%i_eta%i",i,j)]/bTagWeight;
                *f[Form("bTagWeight_HFStats1Down_pt%i_eta%i",i,j)] = bTagWeightSys[Form("bTagWeight_HFStats1Down_pt%i_eta%i",i,j)] /bTagWeight;
                *f[Form("bTagWeight_HFStats2Up_pt%i_eta%i",i,j)] = bTagWeightSys[Form("bTagWeight_HFStats2Up_pt%i_eta%i",i,j)]/bTagWeight;
                *f[Form("bTagWeight_HFStats2Down_pt%i_eta%i",i,j)] = bTagWeightSys[Form("bTagWeight_HFStats2Down_pt%i_eta%i",i,j)]/bTagWeight;
                *f[Form("bTagWeight_cErr1Up_pt%i_eta%i",i,j)] = bTagWeightSys[Form("bTagWeight_cErr1Up_pt%i_eta%i",i,j)]/bTagWeight;
                *f[Form("bTagWeight_cErr1Down_pt%i_eta%i",i,j)] = bTagWeightSys[Form("bTagWeight_cErr1Down_pt%i_eta%i",i,j)]/bTagWeight;
                *f[Form("bTagWeight_cErr2Up_pt%i_eta%i",i,j)] = bTagWeightSys[Form("bTagWeight_cErr2Up_pt%i_eta%i",i,j)]/bTagWeight;
                *f[Form("bTagWeight_cErr2Down_pt%i_eta%i",i,j)] = bTagWeightSys[Form("bTagWeight_cErr2Down_pt%i_eta%i",i,j)]/bTagWeight;
           }
        }
    }

    //MET trigger efficiency
    *f["weight_mettrigSF"] = 1.0;
    *f["weight_mettrigSF_up"] = 1.0;
    *f["weight_mettrigSF_down"] = 1.0;
    if ( mInt("Vtype") == 4 && mInt("sampleIndex")!=0) {
        double met_mht = std::min(m("MET_pt"),m("MHT_pt"));
        auto args = std::vector<double>{met_mht};
        double met_trigger_sf=1.0;
        double met_trigger_sf_up=1.0;
        double met_trigger_sf_down=1.0;
        if (m("dataYear")==2018){
            if(met_mht > 500) met_mht=500;
            if(met_mht >= 100) {
                met_trigger_sf = met_trigger_sf120_2018_func->eval(args.data());
                met_trigger_sf_up = (met_trigger_sf120_2018_func_up->eval(args.data()))/met_trigger_sf;
                met_trigger_sf_down = (met_trigger_sf120_2018_func_down->eval(args.data()))/met_trigger_sf;
            } else {
                 met_trigger_sf=0.;
                 met_trigger_sf_up=1.;
                 met_trigger_sf_down=1.;
            }
        }
        if (m("dataYear")==2017){
            if(met_mht > 500) met_mht=500;
            if(met_mht >= 100) {
                met_trigger_sf = met_trigger_sf120_2017_func->eval(args.data());
                met_trigger_sf_up = (met_trigger_sf120_2017_func_up->eval(args.data()))/met_trigger_sf;
                met_trigger_sf_down = (met_trigger_sf120_2017_func_down->eval(args.data()))/met_trigger_sf;
            } else { 
                met_trigger_sf=0.;
                met_trigger_sf_up=1.;
                met_trigger_sf_down=1.;
            }
        }
        if (m("dataYear")==2016){
            if(met_mht > 500) met_mht=500;
            if(met_mht >= 100){
                met_trigger_sf = met_trigger_sf110OR170_2016_func->eval(args.data());
                met_trigger_sf_up = (met_trigger_sf110OR170_2016_func_up->eval(args.data()))/met_trigger_sf;
                met_trigger_sf_down = (met_trigger_sf110OR170_2016_func_down->eval(args.data()))/met_trigger_sf;
            } else { 
                met_trigger_sf=0.;
                met_trigger_sf_up=1.;
                met_trigger_sf_down=1.;
            }
        }
        *f["weight_mettrigSF"] = met_trigger_sf;
        *f["weight_mettrigSF_up"] = met_trigger_sf_up;
        *f["weight_mettrigSF_down"] = met_trigger_sf_down;
    }
    //if (bool(*f["doCutFlow"])) {
    //    ofile->cd();
    //    outputTree->Fill();
    //    return;
    //}
    // General use variables

    if (cursample->sampleNum == 53 || cursample->sampleNum == 54) {
        *f["nProcEvents"] = cursample->CountWeighted->GetBinContent(1);
    }
    else {
        *f["nProcEvents"] = cursample->processedEvents;
    }
    if (*in["sampleIndex"] != 0) {
        TH1F *CountWeightedLHEWeightScale = cursample->CountWeightedLHEWeightScale;
        TH1F *CountWeightedLHEWeightPdf = cursample->CountWeightedLHEWeightPdf;
        //Need to scale CountWeightedLHEWeightScale and CountWeightedLHEWeightPdf:
        CountWeightedLHEWeightScale->Scale(1./round(CountWeightedLHEWeightScale->GetBinContent(5)));
        CountWeightedLHEWeightPdf->Scale(1./round(CountWeightedLHEWeightScale->GetBinContent(5)));
        for (int i=0; i<mInt("nLHEScaleWeight");i++) {
            f["LHE_weights_scale_normwgt"][i] = 1.0 / CountWeightedLHEWeightScale->GetBinContent(CountWeightedLHEWeightScale->FindBin(i));
            f["LHEScaleWeight"][i] = m("LHEScaleWeight",i)*m("LHE_weights_scale_normwgt",i);
        }
        *f["LHE_weights_scale_muRUp"] = m("LHEScaleWeight",7);
        *f["LHE_weights_scale_muFUp"] = m("LHEScaleWeight",5);
        *f["LHE_weights_scale_muRDown"] = m("LHEScaleWeight",1);
        *f["LHE_weights_scale_muFDown"] = m("LHEScaleWeight",3);
        for (int i=0; i<mInt("nLHEPdfWeight");i++) {
            if (CountWeightedLHEWeightPdf->GetBinContent(CountWeightedLHEWeightPdf->FindBin(i)) != 0) {
                f["LHE_weights_pdf_normwgt"][i] = 1.0 / CountWeightedLHEWeightPdf->GetBinContent(CountWeightedLHEWeightPdf->FindBin(i));
                f["LHEPdfWeight"][i] = m("LHEPdfWeight",i)*m("LHE_weights_pdf_normwgt",i);
            }
            else {
                f["LHE_weights_pdf_normwgt"][i] = 1.0;
            }
        }
    }

    if(mInt("sampleIndex")!=0){
        //if (m("genWeight") > 0) {
        //    *f["weight"] = cursample->intWeight;
        //} else {
        //    *f["weight"] = cursample->intWeight * -1;
       // }
        //if (cursample->sampleNum == 49 || cursample->sampleNum == 491) {
        //    // special prescription for WJets_BGenFilter sample
        //    *f["weight"] = m("weight")*fabs(m("genWeight"));
        //}
        // with NanoAOD we use genEventSumw to normalize which includes absolute value
        // of genWeight
        *f["intWeight"] = cursample->intWeight;
        if(m("dataYear")!=2018){
            *f["weight"] = cursample->intWeight * m("genWeight") * m("PrefireWeight");
            *f["PrefireWeight_Up"] = m("PrefireWeight_Up")/m("PrefireWeight");
            *f["PrefireWeight_Down"] = m("PrefireWeight_Down")/m("PrefireWeight");
         } else {
            *f["weight"] = cursample->intWeight * m("genWeight");
         }
        //std::cout<"intWeight = "<<std::endl;
    }
    else {
        *f["weight"] = 1.0;
        *f["intWeight"] = 1.0;
    }


    *f["weight_PU_2016to2017"] = 1.0;
    if(mInt("sampleIndex")!=0){
        *f["weight_PU_2016to2017"] = puWeight_2016to2017(m("Pileup_nTrueInt"));
        if (m("dataYear")==2016){
            if (m("doICHEP") != 1) {
                *f["weight_PU"] = m("puWeight");
                //TEMPORARY SOLUTION: REMEMBER TO FIX IT BACK ONCE weightUP/DOWN WILL BE IN NANOAOD
                  *f["weight_PUUp"] = puWeight_2016Up(m("Pileup_nTrueInt"));
                  *f["weight_PUDown"] = puWeight_2016Down(m("Pileup_nTrueInt"));
                  //*f["weight_PUUp"] = m("puWeightUp") / m("puWeight");
                //*f["weight_PUUp"] = m("puWeight");
                  //*f["weight_PUDown"] = m("puWeightDown") / m("puWeight");
                //*f["weight_PUDown"] = m("puWeight");
            } else {
                //*f["weight_PU"] = *f["puWeight"];
                //*f["weight_PU"]=ReWeightMC(int(*f["nTrueInt"])+0.5); // it is now a float (continuous distribution?) so we round to the nearest int
                *f["weight_PU"]=puWeight_ichep(int(m("nTrueInt"))); // it is now a float (continuous distribution?) so we round to the nearest int
                *f["weight_PUUp"]=(puWeight_ichep_up(int(m("nTrueInt")))) / m("weight_PU"); // it is now a float (continuous distribution?) so we round to the nearest int
                *f["weight_PUDown"]=(puWeight_ichep_down(int(m("nTrueInt")))) / m("weight_PU"); // it is now a float (continuous distribution?) so we round to the nearest int
            }
        }
        if (m("dataYear")==2017){
            //*f["weight_PU"] = m("puWeight");
            *f["weight_PU"] = GetPUWeight(m("Pileup_nTrueInt"));
            *f["weight_PUUp"] = GetPUWeight(m("Pileup_nTrueInt"),1);
            *f["weight_PUDown"] = GetPUWeight(m("Pileup_nTrueInt"),-1);
        }
        if (m("dataYear")==2018){
            //*f["weight_PU"] = m("puWeight");
	    *f["weight_PU"] = GetPUWeight(m("Pileup_nTrueInt"));
	    *f["weight_PUUp"] = GetPUWeight(m("Pileup_nTrueInt"),1);
	    *f["weight_PUDown"] = GetPUWeight(m("Pileup_nTrueInt"),-1);
        }
        /*if (mInt("nGenTop")==0 && mInt("nGenVbosons")>0) {
            // only apply to Z/W+jet samples
            *f["weight_ptQCD"]=ptWeightQCD(mInt("nGenVbosons"), m("LHE_HT"), mInt("LeadGenVBoson_pdgId"));
        }*/
    } else {
        *f["weight_PU"]=1;
        *f["weight_PUUp"]=1;
        *f["weight_PUUp"]=1;
        *f["puWeight"]=1;
    }

    // From Silvio
    // https://github.com/silviodonato/Xbb/blob/V21/python/ZvvHbb13TeVconfig/samples_nosplit.ini
    // calculated here: https://github.com/silviodonato/Xbb/blob/V21/python/getWeights.py

    float WBjets_ptVMin = 40.;
    float WBjets_ptVMax = 10000000.;
    float WjetsBgen_ptVMin  = 40.;
    float WjetsBgen_ptVMax  = 10000000.;

    float weightWBjetsHT100,weightWBjetsHT200,weightWBjetsHT400,weightWBjetsHT600,weightWBjetsHT800,weightWBjetsHT1200,weightWBjetsHT2500;
    float weightWjetsBgenHT100,weightWjetsBgenHT200,weightWjetsBgenHT400,weightWjetsBgenHT600,weightWjetsBgenHT800,weightWjetsBgenHT1200,weightWjetsBgenHT2500;

    if (int(m("do2015")) == 1) {
        // weights for V21 ntuples (2015 analysis)
        weightWBjetsHT100=    0.22;
        weightWBjetsHT200=    0.34;
        weightWBjetsHT400=    0.62;
        weightWBjetsHT600=    0.73;
        weightWBjetsHT800=    0.73;
        weightWBjetsHT1200=   0.73;
        weightWBjetsHT2500=    0.73;

        weightWjetsBgenHT100=    0.24;
        weightWjetsBgenHT200=    0.39;
        weightWjetsBgenHT400=    0.69;
        weightWjetsBgenHT600=    0.80;
        weightWjetsBgenHT800=    0.80;
        weightWjetsBgenHT1200=    0.80;
        weightWjetsBgenHT2500=    0.80;
    }

    else if (int(m("doICHEP")) == 1) {
        // weights for V24 ntuples (2016 analysis, 22/fb)
        weightWBjetsHT100=    0.50;
        weightWBjetsHT200=    0.67;
        weightWBjetsHT400=    0.86;
        weightWBjetsHT600=    0.99;
        weightWBjetsHT800=    0.94;
        weightWBjetsHT1200=    1.0;
        weightWBjetsHT2500=    1.0;

        weightWjetsBgenHT100=    0.59;
        weightWjetsBgenHT200=    0.75;
        weightWjetsBgenHT400=    0.90;
        weightWjetsBgenHT600=    0.99;
        weightWjetsBgenHT800=    0.95;
        weightWjetsBgenHT1200= 1.0;
        weightWjetsBgenHT2500= 1.0;

        /*// weights for V22 ntuples (2016 analysis)
        weightWBjetsHT100=    0.42;
        weightWBjetsHT200=    0.67;
        weightWBjetsHT400=    0.62;
        weightWBjetsHT600=    0.94;
        weightWBjetsHT800=    0.99;
        weightWBjetsHT1200=    1.0;
        weightWBjetsHT2500=    1.0;
        weightWjetsBgenHT100=    0.51;
        weightWjetsBgenHT200=    0.75;
        weightWjetsBgenHT400=    0.71;
        weightWjetsBgenHT600=    0.95;
        weightWjetsBgenHT800=    0.99;
        weightWjetsBgenHT1200= 1.0;
        weightWjetsBgenHT2500= 1.0;*/
    }

    else {
        WBjets_ptVMin = 100.;
        WBjets_ptVMax = 200.;
        WjetsBgen_ptVMin = 200.;
        WjetsBgen_ptVMax = 10000000.;

        weightWBjetsHT100 = 0.12;
        weightWBjetsHT200 = 0.10;
        weightWBjetsHT400 = 0.13;
        weightWBjetsHT600 = 0.77;
        weightWBjetsHT800 = 0.76;
        weightWBjetsHT1200 = 0.91;
        weightWBjetsHT2500 = 0.99;

        weightWjetsBgenHT100 = 0.0;
        weightWjetsBgenHT200 = 0.035;
        weightWjetsBgenHT400 = 0.052;
        weightWjetsBgenHT600 = 0.55;
        weightWjetsBgenHT800 = 0.54;
        weightWjetsBgenHT1200 = 0.80;
        weightWjetsBgenHT2500 = 0.99;

    }

    // stitch together W b-enriched samples with HT-binned samples in order to maximize statistical power
    float WJetStitchWeight = 1.0;
    int nbHad = 0;
    if (int(m("reRunGenInfo"))==1){
        nbHad = mInt("nGenStatus2bHad_re");
    } else {
        nbHad = mInt("nGenStatus2bHad");
    }
    if ( m("dataYear")==2016 && ((cursample->sampleNum>=40 && cursample->sampleNum<=47) || (cursample->sampleNum>=50 && cursample->sampleNum<=54)) ) {
        if (m("LHE_HT")>100 && m("LHE_HT")<200) {
            if (m("LHE_Vpt") > WBjets_ptVMin && m("LHE_Vpt") < WBjets_ptVMax && mInt("LHE_Nb") > 0) {
                if (cursample->sampleNum == 50) {
                    WJetStitchWeight = (1 - weightWBjetsHT100);
                }
                else {
                    WJetStitchWeight = weightWBjetsHT100;
                }
            }
            else if (m("LHE_Vpt") > WjetsBgen_ptVMin && m("LHE_Vpt") < WjetsBgen_ptVMax && mInt("LHE_Nb") == 0 && nbHad > 0) {
                if (cursample->sampleNum == 53) {
                    WJetStitchWeight = (1 - weightWjetsBgenHT100);
                }
                else {
                    WJetStitchWeight = weightWjetsBgenHT100;
                }
            }
        }
        if (m("LHE_HT")>200 && m("LHE_HT")<400) {
            if (m("LHE_Vpt") > WBjets_ptVMin && m("LHE_Vpt") < WBjets_ptVMax && mInt("LHE_Nb") > 0) {
                if (cursample->sampleNum == 50) {
                    WJetStitchWeight = (1 - weightWBjetsHT200);
                }
                else {
                    WJetStitchWeight = weightWBjetsHT200;
                }
            }
            else if (m("LHE_Vpt") > WjetsBgen_ptVMin && m("LHE_Vpt") < WjetsBgen_ptVMax && mInt("LHE_Nb") == 0 && nbHad > 0) {
                if (cursample->sampleNum == 53) {
                    WJetStitchWeight = (1 - weightWjetsBgenHT200);
                }
                else {
                    WJetStitchWeight = weightWjetsBgenHT200;
                }
            }
        }
        if (m("LHE_HT")>400 && m("LHE_HT")<600) {
            if (m("LHE_Vpt") > WBjets_ptVMin && m("LHE_Vpt") < WBjets_ptVMax && mInt("LHE_Nb") > 0) {
                if (cursample->sampleNum == 50) {
                    WJetStitchWeight = (1 - weightWBjetsHT400);
                }
                else {
                    WJetStitchWeight = weightWBjetsHT400;
                }
            }
            else if (m("LHE_Vpt") > WjetsBgen_ptVMin && m("LHE_Vpt") < WjetsBgen_ptVMax && mInt("LHE_Nb") == 0 && nbHad > 0) {
                if (cursample->sampleNum == 53) {
                    WJetStitchWeight = (1 - weightWjetsBgenHT400);
                }
                else {
                    WJetStitchWeight = weightWjetsBgenHT400;
                }
            }
        }
        if (m("LHE_HT")>600 && m("LHE_HT")<800) {
            if (m("LHE_Vpt")  > WBjets_ptVMin && m("LHE_Vpt") < WBjets_ptVMax && mInt("LHE_Nb") > 0) {
                if (cursample->sampleNum == 50) {
                    WJetStitchWeight = (1 - weightWBjetsHT600);
                }
                else {
                    WJetStitchWeight = weightWBjetsHT600;
                }
            }
        else if (m("LHE_Vpt") > WjetsBgen_ptVMin && m("LHE_Vpt") < WjetsBgen_ptVMax && mInt("LHE_Nb") == 0 && nbHad > 0) {
                if (cursample->sampleNum == 53) {
                    WJetStitchWeight = (1 - weightWjetsBgenHT600);
                }
                else {
                    WJetStitchWeight = weightWjetsBgenHT600;
                }
            }
        }
        if (m("LHE_HT")>800 && m("LHE_HT")<1200) {
            if (m("LHE_Vpt") > WBjets_ptVMin && m("LHE_Vpt") < WBjets_ptVMax && mInt("LHE_Nb") > 0) {
                if (cursample->sampleNum == 50) {
                    WJetStitchWeight = (1 - weightWBjetsHT800);
                }
                else {
                    WJetStitchWeight = weightWBjetsHT800;
                }
            }
            else if (m("LHE_Vpt") > WjetsBgen_ptVMin && m("LHE_Vpt") < WjetsBgen_ptVMax && mInt("LHE_Nb") == 0 && nbHad > 0) {
                if (cursample->sampleNum == 53) {
                    WJetStitchWeight = (1 - weightWjetsBgenHT800);
                }
                else {
                    WJetStitchWeight = weightWjetsBgenHT800;
                }
            }
        }
        if (m("LHE_HT")>1200 && m("LHE_HT")<2500) {
            if (m("LHE_Vpt") > WBjets_ptVMin && m("LHE_Vpt") < WBjets_ptVMax && mInt("LHE_Nb") > 0) {
                if (cursample->sampleNum == 50) {
                    WJetStitchWeight = (1 - weightWBjetsHT1200);
                }
                else {
                    WJetStitchWeight = weightWBjetsHT1200;
                }
            }
            else if (m("LHE_Vpt") > WjetsBgen_ptVMin && m("LHE_Vpt") < WjetsBgen_ptVMax && mInt("LHE_Nb") == 0 && nbHad > 0) {
                if (cursample->sampleNum == 53) {
                    WJetStitchWeight = (1 - weightWjetsBgenHT1200);
                }
                else {
                    WJetStitchWeight = weightWjetsBgenHT1200;
                }
            }
        }
        if (m("LHE_HT")>2500) {
            if (m("LHE_Vpt") > WBjets_ptVMin && m("LHE_Vpt") < WBjets_ptVMax && mInt("LHE_Nb") > 0) {
                if (cursample->sampleNum == 50) {
                    WJetStitchWeight = (1 - weightWBjetsHT2500);
                }
                else {
                    WJetStitchWeight = weightWBjetsHT2500;
                }
            }
            else if (m("LHE_Vpt") > WjetsBgen_ptVMin && m("LHE_Vpt") < WjetsBgen_ptVMax && mInt("LHE_Nb") == 0 && nbHad > 0) {
                if (cursample->sampleNum == 53) {
                    WJetStitchWeight = (1 - weightWjetsBgenHT2500);
                }
                else {
                    WJetStitchWeight = weightWjetsBgenHT2500;
                }
            }
        }
        //*f["weight"] = *f["weight"] * WJetStitchWeight;
    }

    // for now we don't do stitching since the W+jets b-enriched statistics are very high
    *f["WJetStitchWeight"] = WJetStitchWeight;

    *b["usingBEnriched"] =  true; // if using b-enriched need to stitch properly
    if (cursample->sampleNum >= 40 && cursample->sampleNum <=47) {
        if (m("LHE_Vpt") > 100) {
	     if (mInt("LHE_Nb") != 0 || nbHad != 0) *b["usingBEnriched"]=false;
        }
    } else if (cursample->sampleNum == 50) {
        if (m("LHE_Vpt") < 100 || m("LHE_Vpt") > 200 || mInt("LHE_Nb") == 0) *b["usingBEnriched"]=false;
    } else if (cursample->sampleNum == 51) {
        if (m("LHE_Vpt") < 200 || mInt("LHE_Nb") == 0) *b["usingBEnriched"]=false;
    } else if (cursample->sampleNum == 53) {
        //if (m("LHE_Vpt") < 100 || m("LHE_Vpt") > 200 ) *b["usingBEnriched"]=false;
        //if (m("LHE_Vpt") < 100 || m("LHE_Vpt") > 200 || nbHad == 0) *b["usingBEnriched"]=false;
        if (m("LHE_Vpt") < 100 || m("LHE_Vpt") > 200 || nbHad == 0 || mInt("LHE_Nb") != 0) *b["usingBEnriched"]=false;
    } else if (cursample->sampleNum == 54) {
        //if (m("LHE_Vpt") < 200) *b["usingBEnriched"]=false;
        //if (m("LHE_Vpt") < 200 || nbHad == 0) *b["usingBEnriched"]=false;
        if (m("LHE_Vpt") < 200 || nbHad == 0 || mInt("LHE_Nb") != 0) *b["usingBEnriched"]=false;
    }
    if (cursample->sampleNum >= 110 && cursample->sampleNum<=117){
        if (m("LHE_Vpt") > 100) {
	     if (mInt("LHE_Nb") != 0 || nbHad != 0) *b["usingBEnriched"]=false;
        }
    } else if (cursample->sampleNum == 121){
        if (m("LHE_Vpt") < 100 || m("LHE_Vpt") > 200 || mInt("LHE_Nb") == 0) *b["usingBEnriched"]=false;
    } else if (cursample->sampleNum == 122){
        if (m("LHE_Vpt") < 200 || mInt("LHE_Nb") == 0) *b["usingBEnriched"]=false;
    } else if (cursample->sampleNum == 141){
        //if (m("LHE_Vpt") < 100 || m("LHE_Vpt") > 200 || nbHad == 0) *b["usingBEnriched"]=false;
        if (m("LHE_Vpt") < 100 || m("LHE_Vpt") > 200 || nbHad == 0 || mInt("LHE_Nb") != 0) *b["usingBEnriched"]=false;
    } else if (cursample->sampleNum == 142){
        //if (m("LHE_Vpt") < 200 || nbHad == 0) *b["usingBEnriched"]=false;
        if (m("LHE_Vpt") < 200 || nbHad == 0 || mInt("LHE_Nb") != 0) *b["usingBEnriched"]=false;
    }
    if (cursample->sampleNum >= 150 && cursample->sampleNum<=156){
        if (m("LHE_Vpt") > 100) {
	     if (mInt("LHE_Nb") != 0 || nbHad != 0) *b["usingBEnriched"]=false;
        }
    } else if (cursample->sampleNum == 160){
        if (m("LHE_Vpt") < 100 || m("LHE_Vpt") > 200 || mInt("LHE_Nb") == 0) *b["usingBEnriched"]=false;
    } else if (cursample->sampleNum == 161){
        if (m("LHE_Vpt") < 200 || mInt("LHE_Nb") == 0) *b["usingBEnriched"]=false;
    } else if (cursample->sampleNum == 162){
        //if (m("LHE_Vpt") < 100 || m("LHE_Vpt") > 200 || nbHad == 0) *b["usingBEnriched"]=false;
        if (m("LHE_Vpt") < 100 || m("LHE_Vpt") > 200 || nbHad == 0 || mInt("LHE_Nb") != 0) *b["usingBEnriched"]=false;
    } else if (cursample->sampleNum == 163){
        //if (m("LHE_Vpt") < 200 || nbHad == 0) *b["usingBEnriched"]=false;
        if (m("LHE_Vpt") < 200 || nbHad == 0 || mInt("LHE_Nb") != 0) *b["usingBEnriched"]=false;
    }


    *b["useLOVV"] = true; 
    // FIXME I think this needs to be adjusted to really included all the NLO VV, at the time of writing this it is not the full phase space in samples_2016.txt I think
    if (cursample->sampleNum==31 || cursample->sampleNum==32 || cursample->sampleNum==34 || cursample->sampleNum==36 || cursample->sampleNum==37) {
        *b["useLOVV"] = false;
    }
    *b["useNLOVV"] = true; 
    // FIXME I think this needs to be adjusted to really included all the NLO VV, at the time of writing this it is not the full phase space in samples_2016.txt I think
    if (cursample->sampleNum==30 || cursample->sampleNum==33 || cursample->sampleNum==35) {
        *b["useNLOVV"] = false;
    }

    // we need to just save the bTagWeight since we only want to apply it
    // for the nominal shape
    //if (cursyst->name == "nominal") {
    //    *f["weight"] = *f["weight"] * *f["bTagWeight"];
    //}*/


    float VBenrichReweight = 1.0;
    
    if (m("LHE_Vpt") > 100 && mInt("LHE_Nb") == 0 && mInt("nGenStatus2bHad") > 0){
      //ZJets_BGenFilter, WJets_BGenFilter, DYJets_BGenFilter                                                                                                   

      if(m("dataYear") == 2016){
	if (cursample->sampleNum == 162 || cursample->sampleNum == 163){
	  VBenrichReweight = 3*0.92;
	}else if (cursample->sampleNum == 53 || cursample->sampleNum == 54){
	  VBenrichReweight = 0.99;
	}else if (cursample->sampleNum == 141 || cursample->sampleNum == 142){
	  VBenrichReweight = 0.94;
	}
      }else if(m("dataYear") == 2017){
	if (cursample->sampleNum == 162 || cursample->sampleNum == 163){
	  VBenrichReweight = 3.725*(7.7e-01 + 1.184e-03*std::min(m("LHE_Vpt"),800.) - 9.181e-07*TMath::Power(std::min(m("LHE_Vpt"),800.),2));
	}else if (cursample->sampleNum == 53 || cursample->sampleNum == 54){
	  VBenrichReweight = 1.248*(8.325e-01 + 1.054e-03*std::min(m("LHE_Vpt"),800.) - 1.067e-06*TMath::Power(std::min(m("LHE_Vpt"),800.),2));
	}else if (cursample->sampleNum == 141 || cursample->sampleNum == 142){
	  VBenrichReweight = 1.171*(7.825e-01 + 1.529e-03*std::min(m("LHE_Vpt"),800.) - 9.667e-07*TMath::Power(std::min(m("LHE_Vpt"),800.),2));
	}
      }else if(m("dataYear") == 2018){
	if (cursample->sampleNum == 162 || cursample->sampleNum == 163){
	  VBenrichReweight = 3.768*(3.17457e-00 + 2.671e-03*std::min(m("LHE_Vpt"),800.) - 7.573e-07*TMath::Power(std::min(m("LHE_Vpt"),800.),2));
	}else if (cursample->sampleNum == 53 || cursample->sampleNum == 54){
	  VBenrichReweight = 1.2656*(1.05375 + 1.44192e-03*std::min(m("LHE_Vpt"),800.) - 1.80578e-06*TMath::Power(std::min(m("LHE_Vpt"),800.),2));
	}else if (cursample->sampleNum == 141 || cursample->sampleNum == 142){
	  VBenrichReweight = 1.1796*(8.4438e-01 + 1.02536e-03*std::min(m("LHE_Vpt"),800.) - 2.76698e-07*TMath::Power(std::min(m("LHE_Vpt"),800.),2));
	}
      }
    }else if (m("LHE_Vpt") > 100 && mInt("LHE_Nb") > 0){
      //ZBJets, WBJets, DYBJets                                                                                                                                            
      if(m("dataYear") == 2016){
	if (cursample->sampleNum == 160 || cursample->sampleNum == 161){
	  VBenrichReweight = 1.02;
	}else if (cursample->sampleNum == 50 || cursample->sampleNum == 51){
	  VBenrichReweight = 1.04;
	}else if (cursample->sampleNum == 121 || cursample->sampleNum == 122){
	  VBenrichReweight = 1.04;
	}
      }else if(m("dataYear") == 2017){
	if (cursample->sampleNum == 160 || cursample->sampleNum == 161){
	  VBenrichReweight = 1.332*(6.968e-01 + 1.764e-03*std::min(m("LHE_Vpt"),800.) - 1.526e-06*TMath::Power(std::min(m("LHE_Vpt"),800.),2));
	}else if (cursample->sampleNum == 50 || cursample->sampleNum == 51){
	  VBenrichReweight = 0.977*(1.005e-00 + 5.043e-04*std::min(m("LHE_Vpt"),800.) - 3.894e-07*TMath::Power(std::min(m("LHE_Vpt"),800.),2));
	}else if (cursample->sampleNum == 121 || cursample->sampleNum == 122){
	  VBenrichReweight = 1.259*(7.519e-01 + 1.975e-03*std::min(m("LHE_Vpt"),800.) - 1.836e-06*TMath::Power(std::min(m("LHE_Vpt"),800.),2));
	}
      }else if(m("dataYear") == 2018){
	if (cursample->sampleNum == 160 || cursample->sampleNum == 161){
	  VBenrichReweight = 1.337*(9.955e-01 + 1.842e-03*std::min(m("LHE_Vpt"),800.) - 1.215e-06*TMath::Power(std::min(m("LHE_Vpt"),800.),2));
	}else if (cursample->sampleNum == 50 || cursample->sampleNum == 51){
	  VBenrichReweight = 1.0079*(7.97747e-01 + 1.51461e-03*std::min(m("LHE_Vpt"),800.) - 9.57889e-07*TMath::Power(std::min(m("LHE_Vpt"),800.),2));
	}else if (cursample->sampleNum == 121 || cursample->sampleNum == 122){
	  VBenrichReweight = 1.2474*(7.97747e-01 + 1.51461e-03*std::min(m("LHE_Vpt"),800.) - 9.57889e-07*TMath::Power(std::min(m("LHE_Vpt"),800.),2));
	}
      }
    }
    
    *f["VBenrichReweight"] = VBenrichReweight;
    *f["weight"] = *f["weight"] * *f["VBenrichReweight"];
    




    // Split WJets and ZJets samples by jet parton flavor
    *in["bMCFlavorSum"] = 0;
    *in["bMCFlavorSumSelected"] = 0;
    *in["bGenJetSum"] = 0;
    *in["bGenJetBSum"] = 0;
    *in["bGenBJetSum"] = 0.;

    *in["sampleIndex_sel"] =
    *in["sampleIndex_GenJetSum"] =
    *in["sampleIndex_GenJetSumNB"] =
    *in["sampleIndex_GenBJetSum"] =
    *in["sampleIndex"];

    if (cursample->doJetFlavorSplit) {
        *in["sampleIndex"] = mInt("sampleIndex")*100;
        *in["sampleIndex_sel"] = mInt("sampleIndex_sel")*100;
        *in["sampleIndex_GenJetSum"] = mInt("sampleIndex_GenJetSum")*100;
        *in["sampleIndex_GenBJetSum"] = mInt("sampleIndex_GenBJetSum")*100;
        *in["sampleIndex_GenJetSumNB"] = mInt("sampleIndex_GenJetSumNB")*100;

        if (fabs(mInt("Jet_hadronFlavour",mInt("hJetInd1"))) == 5)  *in["bMCFlavorSumSelected"]=mInt("bMCFlavorSumSelected")+1;
        if (fabs(mInt("Jet_hadronFlavour",mInt("hJetInd2"))) == 5)  *in["bMCFlavorSumSelected"]=mInt("bMCFlavorSumSelected")+1;

        for(int iJet=0;iJet<mInt("nJet");iJet++){
            if(fabs(mInt("Jet_hadronFlavour",iJet))==5){
               *in["bMCFlavorSum"]=mInt("bMCFlavorSum")+1;
            }
        }


        for(int indGJ=0; indGJ<mInt("nGenJet"); indGJ++){
            *in["bGenJetSum"]=mInt("bGenJetSum")+1;
        }

        /*FIXME we might want to add this back in - do the pruned gen parts store enough information to do this?,
        for(int indGJ=0; indGJ<mInt("nGenJet"); indGJ++){
            *in["bGenJetBSum"]=mInt("bGenJetBSum")+mInt("GenJet_numBHadrons",indGJ);
        }*/
        for(int indGJ=0; indGJ<mInt("nGenJet"); indGJ++){
            if (mInt("GenJet_hadronFlavour",indGJ) == 5 && m("GenJet_pt",indGJ)>20 && fabs(m("GenJet_eta",indGJ)<2.4)) {
                *in["bGenBJetSum"]=mInt("bGenBJetSum") + 1;
            }
        }

        //if(*in["bMCFlavorSum"]==1){
        //    *in["sampleIndex"] = *in["sampleIndex"] + 1;
        //}else if(*in["bMCFlavorSum"]>1){
        //    *in["sampleIndex"] = *in["sampleIndex"] + 2;
        //}

        if(mInt("bMCFlavorSumSelected")==1){
            *in["sampleIndex_sel"] = mInt("sampleIndex_sel") + 1;
        }else if(mInt("bMCFlavorSumSelected")>1){
            *in["sampleIndex_sel"] = mInt("sampleIndex_sel") + 2;
        }

        if(mInt("bGenJetBSum")==1){
            *in["sampleIndex_GenJetSum"] = mInt("sampleIndex_GenJetSum") + 1;
        }else if(mInt("bGenJetBSum")>1){
            *in["sampleIndex_GenJetSum"] = mInt("sampleIndex_GenJetSum") + 2;
        }

        if(mInt("bGenBJetSum")==1){
            *in["sampleIndex_GenBJetSum"] = mInt("sampleIndex_GenBJetSum") + 1;
        }else if(mInt("bGenBJetSum")>1){
            *in["sampleIndex_GenBJetSum"] = mInt("sampleIndex_GenBJetSum") + 2;
        }

        if(mInt("bGenJetSum")==1){
            *in["sampleIndex_GenJetSumNB"] = mInt("sampleIndex_GenJetSumNB") + 1;
        }else if(mInt("bGenJetSum")>1){
            *in["sampleIndex_GenJetSumNB"] = mInt("sampleIndex_GenJetSumNB") + 2;
        }

        // default now the number of b's counted
        //*in["sampleIndex"]=*in["sampleIndex_GenJetSum"];
        *in["sampleIndex"]=mInt("sampleIndex_GenBJetSum");
    }

    // Split by boost category
    if(mInt("isWenu")) {
        if(m("V_pt") >= 100 && m("V_pt") < 150) *in["eventClass"] += 2;
        else if(m("V_pt") >= 150) *in["eventClass"] += 1;
        else *in["eventClass"] += 3;
    }
    else if(mInt("isWmunu")) {
        if(m("V_pt") >= 100 && m("V_pt") < 130) *in["eventClass"] += 3;
        else if(m("V_pt") >= 130 && m("V_pt") < 180) *in["eventClass"] += 2;
        else if(m("V_pt") >= 180) *in["eventClass"] += 1;
        else *in["eventClass"] += 4;
    }

    //if (mInt("isWmunu")) {
    //    *f["selLeptons_relIso_0"] = m("selLeptons_relIso04",mInt("lepInd1"));
    //}
    //else if (mInt("isWenu")) {
    //    *f["selLeptons_relIso_0"] = m("selLeptons_relIso03",mInt("lepInd1"));
    //}
    //else if (mInt("isZmm")) {
    //    *f["selLeptons_relIso_0"] = m("selLeptons_relIso04",mInt("lepInd1"));
    //    *f["selLeptons_relIso_1"] = m("selLeptons_relIso04",mInt("lepInd2"));
    //}
    //else if (mInt("isZee")) {
    //    *f["selLeptons_relIso_0"] = m("selLeptons_relIso03",mInt("lepInd1"));
    //    *f["selLeptons_relIso_1"] = m("selLeptons_relIso03",mInt("lepInd2"));
    //}


    //if(mInt("sampleIndex")!=0){
    //    if (mInt("nGenLep") == 0) {
    //        // gen lep is originally a tau?
    //        *f["selLeptons_genLepDR_0"] = -1;
    //    } else {
    //        TLorentzVector GenLep, El;
    //        GenLep.SetPtEtaPhiM(m("GenLep_pt",0),m("GenLep_eta",0),m("GenLep_phi",0),m("GenLep_mass",0));
    //        El.SetPtEtaPhiM(m("selLeptons_pt",mInt("lepInd1")), m("selLeptons_eta",mInt("lepInd1")), m("selLeptons_phi",mInt("lepInd1")), m("selLeptons_mass",mInt("lepInd1")));
    //        *f["selLeptons_genLepDR_0"] = El.DeltaR(GenLep);
    //    }
    //}

    //// Reconstruct Higgs and W and recalculate variables ourselves
    //if(debug>1000) std::cout<<"Making composite candidates"<<std::endl;
    //TLorentzVector MET,Lep,W,HJ1,HJ2,Hbb;
    //MET.SetPtEtaPhiM(m("MET_pt"), 0., m("MET_phi"), 0.); // Eta/M don't affect calculation of W.pt and W.phi
    //Lep.SetPtEtaPhiM(m("selLeptons_pt",mInt("lepInd1")), m("selLeptons_eta",mInt("lepInd1")), m("selLeptons_phi",mInt("lepInd1")), m("selLeptons_mass",mInt("lepInd1")));
    //W = MET + Lep;

    //HJ1.SetPtEtaPhiM(m("Jet_bReg",mInt("hJetInd1")), m("Jet_eta",mInt("hJetInd1")), m("Jet_phi",mInt("hJetInd1")), m("Jet_mass",mInt("hJetInd1")) * (m("Jet_bReg",mInt("hJetInd1")) / m("Jet_pt",mInt("hJetInd1")) ) );
    //HJ2.SetPtEtaPhiM(m("Jet_bReg",mInt("hJetInd2")), m("Jet_eta",mInt("hJetInd2")), m("Jet_phi",mInt("hJetInd2")), m("Jet_mass",mInt("hJetInd2")) * (m("Jet_bReg",mInt("hJetInd2")) / m("Jet_pt",mInt("hJetInd2")) ) );
    //Hbb = HJ1 + HJ2;

    // We already calculate these in Analyze()
    //*d["H_mass"] = Hbb.M();
    //*d["H_pt"] = Hbb.Pt();
    //*d["V_pt"] = W.Pt();
    //*d["HVdPhi"] = Hbb.DeltaPhi(W);

    // Set variables used by the BDT
    //*f["H_mass_f"] = (float) *f["H_mass"];
    //*f["H_pt_f"] = (float) *f["H_pt"];
    //*f["V_pt_f"] = (float) *f["V_pt"];

    // General BDT inputs

    *f["hJets_btagged_0"] = (float) m(taggerName,mInt("hJetInd1"));
    *f["hJets_btagged_1"] = (float) m(taggerName,mInt("hJetInd2"));

    *f["hJets_btagWP_0"] = (float) BtagWPForJet(mInt("hJetInd1"));
    *f["hJets_btagWP_1"] = (float) BtagWPForJet(mInt("hJetInd2"));

    //*f["hJets_mt_0"] = HJ1.Mt();
    //*f["hJets_mt_1"] = HJ2.Mt();
    //*f["H_dR"] = (float) HJ1.DeltaR(HJ2);
    //*f["absDeltaPullAngle"] = 0.; //FIXME what is this in the new ntuples??
    *f["hJets_pt_0"] = (float) m("Jet_bReg",mInt("hJetInd1"));
    *f["hJets_pt_1"] = (float) m("Jet_bReg",mInt("hJetInd2"));

    //redefine the leading and subleading jet pt in case of FSR
    if(!m("recoFSR")){
      *f["hJets_leadingPt_noFSR"]    = std::max(*f["hJets_pt_0"],*f["hJets_pt_1"]);
      *f["hJets_subleadingPt_noFSR"] = std::min(*f["hJets_pt_0"],*f["hJets_pt_1"]);
    }else{
      *f["hJets_leadingPt"]    = std::max(*f["HJ1_pt"],*f["HJ2_pt"]);
      *f["hJets_subleadingPt"] = std::min(*f["HJ1_pt"],*f["HJ2_pt"]);
    }

    // Channel specific BDT inputs
   

    std::string bdt_branch_label;
    if(mInt("isZnn")) {
        if(debug>10000) {
            std::cout<<"setting up bdt inputs for isZnn"<<std::endl;
        }
        std::vector<std::string> bdtNames;
        bdtNames.clear();

        if(m("twoResolvedJets")){
            thisBDTInfo = bdtInfos.find("bdt_0lep");
            if(thisBDTInfo != bdtInfos.end()){
                bdtNames.push_back("bdt_0lep");
            }
            thisBDTInfo = bdtInfos.find("bdt_0lep_vzbb");
            if(thisBDTInfo != bdtInfos.end()){
                bdtNames.push_back("bdt_0lep_vzbb");
            }
            thisBDTInfo = bdtInfos.find("bdt_0lep_massless");
            if(thisBDTInfo != bdtInfos.end()){
                float logMjj = log(m("H_mass"));
                *f["H_pt_overLogM"] = m("H_pt")/logMjj;
                *f["hJets_leadingPt_overLogM"]    = m("hJets_leadingPt")/logMjj;
                *f["hJets_subleadingPt_overLogM"] = m("hJets_subleadingPt")/logMjj;
                bdtNames.push_back("bdt_0lep_massless");
            }
        }

        if(*b["oneMergedJet"]){
            thisBDTInfo = bdtInfos.find("BDT_0lep_boosted");
            if(thisBDTInfo != bdtInfos.end()){
                bdtNames.push_back("BDT_0lep_boosted");
            }
            thisBDTInfo = bdtInfos.find("BDT_1lep_boosted");
            if(thisBDTInfo != bdtInfos.end()){
                bdtNames.push_back("BDT_1lep_boosted");
            }
            thisBDTInfo = bdtInfos.find("BDT_2lep_boosted");
            if(thisBDTInfo != bdtInfos.end()){
                bdtNames.push_back("BDT_2lep_boosted");
            }
        }

        if(bdtNames.size()>0){
            //loop through jets for best of other variables
            *f["otherJetsBestBtag"]    = -99;
            *f["otherJetsHighestPt"]   = -99;
            *f["minDPhiFromOtherJets"] = 99;
            for(int iJet=0; iJet<mInt("nJet"); iJet++){
                if(iJet==mInt("hJetInd1")) continue;
                if(iJet==mInt("hJetInd2")) continue;
                if(mInt("Jet_lepFilter",iJet)==0) continue;
                if( (mInt("Jet_puId",iJet) > 6 || m("Jet_Pt",iJet)>50) && m("Jet_bReg",iJet)>25){
                    if(*f["otherJetsBestBtag"]< m(taggerName,iJet)){
                        *f["otherJetsBestBtag"]=m(taggerName,iJet);
                    }
                    if(*f["otherJetsHighestPt"]< m("Jet_bReg",iJet)){
                        *f["otherJetsHighestPt"]=m("Jet_bReg",iJet);
                    }
                }
                if( (mInt("Jet_puId",iJet)>6 || m("Jet_Pt",iJet)>50) && m("Jet_Pt",iJet)>30){
                    if(*f["minDPhiFromOtherJets"]>fabs(EvalDeltaPhi(m("MET_Phi"), m("Jet_phi",iJet)))){
                        *f["minDPhiFromOtherJets"]=fabs(EvalDeltaPhi(m("MET_Phi"), m("Jet_phi",iJet)));
                    }
                }
            }

            for(unsigned int iBDT=0; iBDT<bdtNames.size(); iBDT++){
                std::string bdtname(bdtNames[iBDT]);
                if(debug>5000) {
                    std::cout<<"Evaluating BDT... "<<bdtNames[iBDT]<<std::endl;
                    PrintBDTInfoValues(bdtInfos[bdtNames[iBDT]]);
                    std::cout<<"BDT evaluates to: "<<EvaluateMVA(bdtInfos[bdtNames[iBDT]])<<std::endl;
                }
                bdt_branch_label = bdtInfos[bdtNames[iBDT]]->bdtname;
                if (cursyst->name != "nominal") {
                    bdt_branch_label.append("_");
                    bdt_branch_label.append(cursyst->name);
                }
                *f[bdt_branch_label] = EvaluateMVA(bdtInfos[bdtNames[iBDT]]);
            }
        }

    } else if(mInt("isWenu") || mInt("isWmunu")) {
      //        *f["nAddJets302p5_puid_f"] = (float) mInt("nAddJets302p5_puid");
        *f["nAddJet_f"] = (float) mInt("nAddJets252p9_puid");
        *f["nAddLep_f"] = (float) mInt("nAddLeptons");
        *f["isWenu_f"] = (float) mInt("isWenu");
        *f["isWmunu_f"] = (float) mInt("isWmunu");

        std::vector<std::string> bdtNames;
        bdtNames.clear();
        thisBDTInfo = bdtInfos.find("bdt_1lep");

        if(m("twoResolvedJets")){
            if(thisBDTInfo != bdtInfos.end()){
                bdtNames.push_back("bdt_1lep");
            }
            thisBDTInfo = bdtInfos.find("bdt_1lep_vzbb");
            if(thisBDTInfo != bdtInfos.end()){
                bdtNames.push_back("bdt_1lep_vzbb");
            }
        }

        if(*b["oneMergedJet"]){
            thisBDTInfo = bdtInfos.find("bdt_boosted_wdB");
            if(thisBDTInfo != bdtInfos.end()){
                bdtNames.push_back("bdt_boosted_wdB");
            }
            thisBDTInfo = bdtInfos.find("bdt_boosted_nodB");
            if(thisBDTInfo != bdtInfos.end()){
                bdtNames.push_back("bdt_boosted_nodB");
            }
        }

        if(bdtNames.size()>0){
            for(unsigned int iBDT=0; iBDT<bdtNames.size(); iBDT++){
                std::string bdtname(bdtNames[iBDT]);
                if(debug>5000) {
                    std::cout<<"Evaluating BDT... "<<bdtNames[iBDT]<<std::endl;
                    PrintBDTInfoValues(bdtInfos[bdtNames[iBDT]]);
                    std::cout<<"BDT evaluates to: "<<EvaluateMVA(bdtInfos[bdtNames[iBDT]])<<std::endl;
                }
                bdt_branch_label = bdtInfos[bdtNames[iBDT]]->bdtname;
                if (cursyst->name != "nominal") {
                    bdt_branch_label.append("_");
                    bdt_branch_label.append(cursyst->name);
                }
                *f[bdt_branch_label] = EvaluateMVA(bdtInfos[bdtNames[iBDT]]);
            }
        }
    } else if(mInt("isZee") || mInt("isZmm")) {

        std::vector<std::string> bdtNames;
        bdtNames.clear();
        thisBDTInfo = bdtInfos.find("bdt_2lep_highPt");

        if(m("twoResolvedJets")){
            if(thisBDTInfo != bdtInfos.end()){
                if(m("V_pt")>=150){
                    bdtNames.push_back("bdt_2lep_highPt");
                }
            }
            thisBDTInfo = bdtInfos.find("bdt_2lep_lowPt");
            if(thisBDTInfo != bdtInfos.end()){
                if(m("V_pt")<150){
                    bdtNames.push_back("bdt_2lep_lowPt");
                }
            }
            thisBDTInfo = bdtInfos.find("bdt_2lep_vzbb");
            if(thisBDTInfo != bdtInfos.end()){
                bdtNames.push_back("bdt_2lep_vzbb");
            }
        }

        if(*b["oneMergedJet"]){
            thisBDTInfo = bdtInfos.find("bdt_boosted_wdB");
            if(thisBDTInfo != bdtInfos.end()){
                bdtNames.push_back("bdt_boosted_wdB");
            }
            thisBDTInfo = bdtInfos.find("bdt_boosted_nodB");
            if(thisBDTInfo != bdtInfos.end()){
                bdtNames.push_back("bdt_boosted_nodB");
            }
        }

        if(bdtNames.size()>0){
            for(unsigned int iBDT=0; iBDT<bdtNames.size(); iBDT++){
                std::string bdtname(bdtNames[iBDT]);
                if(debug>5000) {
                    std::cout<<"Evaluating BDT... "<<bdtNames[iBDT]<<std::endl;
                    PrintBDTInfoValues(bdtInfos[bdtNames[iBDT]]);
                    std::cout<<"BDT evaluates to: "<<EvaluateMVA(bdtInfos[bdtNames[iBDT]])<<std::endl;
                }
                bdt_branch_label = bdtInfos[bdtNames[iBDT]]->bdtname;
                if (cursyst->name != "nominal") {
                    bdt_branch_label.append("_");
                    bdt_branch_label.append(cursyst->name);
                }
                *f[bdt_branch_label] = EvaluateMVA(bdtInfos[bdtNames[iBDT]]);
            }
        }
    } else {
        std::cout << mInt("Vtype") << " "
                  << mInt("isZnn") << " "
                  << mInt("isWmunu") << " "
                  << mInt("isWenu") << " "
                  << mInt("isZmm") << " "
                  << mInt("isZee") << " "
                  << mInt("controlSample")
                  << std::endl;
        std::cout<<"ALL CHANNELS APPEAR TO BE FALSE"<<std::endl;
    }

    // FIXME For the code to be meaningful it should go far earlier.
    if (f.count("bdt_bjetreg")>0 && m("twoResolvedJets")) {
        if(debug>10000) {
            std::cout<<"Evaluating the Jet Energy Regression..."<<std::endl;
            PrintBDTInfoValues(bdtInfos["bdt_bjetreg"]);
        }
        double r1Pt = evaluateRegression(mInt("hJetInd1"));
        double r2Pt = evaluateRegression(mInt("hJetInd2"));

        *f["Jet1_regWeight"] = r1Pt/(m("hJets_pt_0"));
        *f["Jet2_regWeight"] = r2Pt/(m("hJets_pt_1"));

        TLorentzVector hJ1_reg = TLorentzVector();
        TLorentzVector hJ2_reg = TLorentzVector();

        *f["Jet1_pt_reg"] = r1Pt;
        *f["Jet2_pt_reg"] = r2Pt;
        hJ1_reg.SetPtEtaPhiM(r1Pt, m("Jet_eta",mInt("hJetInd1")), m("Jet_phi",mInt("hJetInd1")), m("Jet_mass",mInt("hJetInd1")) * (r1Pt/m("Jet_Pt",mInt("hJetInd1"))) );
        hJ2_reg.SetPtEtaPhiM(r2Pt, m("Jet_eta",mInt("hJetInd2")), m("Jet_phi",mInt("hJetInd2")), m("Jet_mass",mInt("hJetInd2")) * (r2Pt/m("Jet_Pt",mInt("hJetInd2"))) );

        //*f["hJets_mt_0"] = hJ1_reg.Mt();
        //*f["hJets_mt_1"] = hJ2_reg.Mt();

        TLorentzVector H_vec_reg = hJ1_reg + hJ2_reg;
        *f["H_mass_reg"] = H_vec_reg.M();
        *f["H_pt_reg"] = H_vec_reg.Pt();
        *f["H_eta_reg"] = H_vec_reg.Eta();
        *f["H_phi_reg"] = H_vec_reg.Phi();
        *f["H_dR_reg"] = hJ1_reg.DeltaR(hJ2_reg);
        *f["H_dPhi_reg"] = hJ1_reg.DeltaPhi(hJ2_reg);
        *f["H_dEta_reg"] = fabs(hJ1_reg.Eta() - hJ2_reg.Eta() );
    }

    // add control sample fitted scale factor (already computed)
    // FIXME beautify later.  maybe put into scale factor container?
    *f["CS_SF"] = 1.0;
    if (mInt("sampleIndex")==4000 || mInt("sampleIndex")==4100 || mInt("sampleIndex")==4200 || mInt("sampleIndex")==4300 || mInt("sampleIndex")==4400 || mInt("sampleIndex")==4500 || mInt("sampleIndex")==4600 || mInt("sampleIndex") == 4700 || mInt("sampleIndex") == 5000 || mInt("sampleIndex")==5100 || mInt("sampleIndex")==5300 || mInt("sampleIndex")==5400) {
        *f["CS_SF"] = m("SF_Wj0b");
    }
    /*else if (*in["sampleIndex"]==2201 || *in["sampleIndex"]==4401 || *in["sampleIndex"]==4501 || *in["sampleIndex"]==4601 || *in["sampleIndex"]==4701 || *in["sampleIndex"]==4801 || *in["sampleIndex"]==4901 || *in["sampleIndex"]==2202 || *in["sampleIndex"]==4402 || *in["sampleIndex"]==4502 || *in["sampleIndex"]==4602 || *in["sampleIndex"]==4702 || *in["sampleIndex"]==4802 || *in["sampleIndex"]==4902) {
        *f["CS_SF"] = *f["SF_WHF"];
    }*/
    else if (mInt("sampleIndex")==4001 || mInt("sampleIndex")==4101 || mInt("sampleIndex")==4201 || mInt("sampleIndex")==4301 || mInt("sampleIndex")==4401 || mInt("sampleIndex")==4501 || mInt("sampleIndex")==4601 || mInt("sampleIndex") == 4701 || mInt("sampleIndex") == 5001 || mInt("sampleIndex")==5101 || mInt("sampleIndex")==5301 || mInt("sampleIndex")==5401) {
        *f["CS_SF"] = m("SF_Wj1b");
    }
    else if (mInt("sampleIndex")==4002 || mInt("sampleIndex")==4102 || mInt("sampleIndex")==4202 || mInt("sampleIndex")==4302 || mInt("sampleIndex")==4402 || mInt("sampleIndex")==4502 || mInt("sampleIndex")==4602 || mInt("sampleIndex") == 4702 || mInt("sampleIndex") == 5002 || mInt("sampleIndex")==5102 || mInt("sampleIndex")==5302 || mInt("sampleIndex")==5402) {
        *f["CS_SF"] = m("SF_Wj2b");
    }
    else if (mInt("sampleIndex")>=200 && mInt("sampleIndex")<=210 ) {
        *f["CS_SF"] = m("SF_TT");
    }

    if (mInt("sampleIndex")!=0) {
        *f["Lep_SF"] = 1.0;
        if (mInt("isZmm") == 1) {
            if(m("dataYear") == 2017) {
                *f["Lep_SF"] = m("SF_DoubleMuId",mInt("lepInd1")) * m("SF_DoubleMuIso",mInt("lepInd1")) *  m("SF_DoubleMuId",mInt("lepInd2")) * m("SF_DoubleMuIso",mInt("lepInd2"));
                //m("SF_DoubleMuTriggerLeg1",mInt("lepInd1")) * m("SF_DoubleMuTriggerLeg2",mInt("lepInd2"));
            }
            if(m("dataYear") == 2018) {
                //TODO: Add the DZ and Mass cut eff
	      *f["Lep_SF"] = m("SF_DoubleMu_ID_AD2018",mInt("lepInd1")) * m("SF_DoubleMu_ISO_AD2018",mInt("lepInd1")) *  m("SF_DoubleMu_ID_AD2018",mInt("lepInd2")) * m("SF_DoubleMu_ISO_AD2018",mInt("lepInd2")) * computeEventSFForDoubleLeptonTrig("SF_Mu8Leg_Data", "SF_Mu17Leg_Data", "SF_Mu8Leg_MC", "SF_Mu17Leg_MC");
            }
	    if(m("dataYear") == 2016) {
                //TODO: Add DZ and Mass cut eff
	      *f["Lep_SF"] = ((20.1/36.4) *m("SF_DoubleMu_ID_BF2016",mInt("lepInd1")) * m("SF_DoubleMu_ISO_BF2016",mInt("lepInd1")) *  m("SF_DoubleMu_ID_BF2016",mInt("lepInd2"))* m("SF_DoubleMu_ISO_BF2016",mInt("lepInd2")) * computeEventSFForDoubleLeptonTrig("SF_Mu8Leg_BF_Data", "SF_Mu17Leg_BF_Data", "SF_Mu8Leg_BF_MC", "SF_Mu17Leg_BF_MC") + (16.3/36.4) * m("SF_DoubleMu_ID_GH2016",mInt("lepInd1")) * m("SF_DoubleMu_ISO_GH2016",mInt("lepInd1")) * m("SF_DoubleMu_ID_GH2016",mInt("lepInd2"))* m("SF_DoubleMu_ISO_GH2016",mInt("lepInd2")) *computeEventSFForDoubleLeptonTrig("SF_Mu8Leg_GH_Data", "SF_Mu17Leg_GH_Data", "SF_Mu8Leg_GH_MC", "SF_Mu17Leg_GH_MC")  ) ;
            }
        }else
        if (mInt("isZee") == 1) {
            if(m("dataYear") == 2017) {
                *f["Lep_SF"] = m("SF_DoubleElIdIso",mInt("lepInd1")) * m("SF_DoubleElIdIso",mInt("lepInd2")) * m("SF_DoubleElTriggerLeg1",mInt("lepInd1")) * m("SF_DoubleElTriggerLeg2",mInt("lepInd2"));
                *f["Lep_SFUp"] = (m("SF_DoubleElIdIso",mInt("lepInd1"))+m("SF_DoubleElIdIso_err",mInt("lepInd1"))) * (m("SF_DoubleElIdIso",mInt("lepInd2"))+m("SF_DoubleElIdIso_err",mInt("lepInd2"))) * (m("SF_DoubleElTriggerLeg1",mInt("lepInd1"))+m("SF_DoubleElTriggerLeg1_err",mInt("lepInd1"))) * (m("SF_DoubleElTriggerLeg2",mInt("lepInd2"))+m("SF_DoubleElTriggerLeg2_err",mInt("lepInd2")));
                *f["Lep_SFDown"] = (m("SF_DoubleElIdIso",mInt("lepInd1"))-m("SF_DoubleElIdIso_err",mInt("lepInd1"))) * (m("SF_DoubleElIdIso",mInt("lepInd2"))-m("SF_DoubleElIdIso_err",mInt("lepInd2"))) * (m("SF_DoubleElTriggerLeg1",mInt("lepInd1"))-m("SF_DoubleElTriggerLeg1_err",mInt("lepInd1"))) * (m("SF_DoubleElTriggerLeg2",mInt("lepInd2"))-m("SF_DoubleElTriggerLeg2_err",mInt("lepInd2")));
            }
        }else if (mInt("isWmunu") == 1) {
            if (m("do2015") == 1) {
                // used for 2015 analysis
                *f["Lep_SF"] = m("selLeptons_SF_IsoTight",mInt("lepInd1")) * m("selLeptons_SF_IdCutTight",mInt("lepInd1")) * m("selLeptons_SF_HLT_RunD4p3",mInt("lepInd1"));
            } else if(m("dataYear") == 2016) {
                *f["Lep_SF"] = ( (20.1/36.4) * m("SF_MuIDTightBCDEF",mInt("lepInd1")) + (16.3/36.4) * m("SF_MuIDTightGH",mInt("lepInd1"))) * ( (20.1/36.4) * m("SF_MuIsoTightBCDEF",mInt("lepInd1")) + (16.3/36.4) * m("SF_MuIsoTightGH",mInt("lepInd1")) ) *  ( (20.1/36.4) * m("SF_MuTriggerBCDEF",mInt("lepInd1")) + (16.3/36.4) * m("SF_MuTriggerGH",mInt("lepInd1")));
                *f["Lep_SFUp"] = ( (20.1/36.4) * (m("SF_MuIDTightBCDEF",mInt("lepInd1")) + m("SF_MuIDTightBCDEF_err",mInt("lepInd1"))) + (16.3/36.4) * (m("SF_MuIDTightGH",mInt("lepInd1")) + m("SF_MuIDTightGH_err",mInt("lepInd1"))) ) * ( (20.1/36.4) * (m("SF_MuIsoTightBCDEF",mInt("lepInd1")) + m("SF_MuIsoTightBCDEF_err",mInt("lepInd1")) ) + (16.3/36.4) * (m("SF_MuIsoTightGH",mInt("lepInd1")) + m("SF_MuIsoTightGH_err",mInt("lepInd1"))) ) *  ( (20.1/36.4) * (m("SF_MuTriggerBCDEF",mInt("lepInd1")) + m("SF_MuTriggerBCDEF_err",mInt("lepInd1")) )+ (16.3/36.4) * (m("SF_MuTriggerGH",mInt("lepInd1"))) + m("SF_MuTriggerGH_err",mInt("lepInd1"))) ;
                *f["Lep_SFUp"] = m("Lep_SFUp")/ m("Lep_SF");
                *f["Lep_SFDown"] = ( (20.1/36.4) * (m("SF_MuIDTightBCDEF",mInt("lepInd1")) - m("SF_MuIDTightBCDEF_err",mInt("lepInd1"))) + (16.3/36.4) * (m("SF_MuIDTightGH",mInt("lepInd1")) - m("SF_MuIDTightGH_err",mInt("lepInd1"))) ) * ( (20.1/36.4) * (m("SF_MuIsoTightBCDEF",mInt("lepInd1")) - m("SF_MuIsoTightBCDEF_err",mInt("lepInd1")) ) + (16.3/36.4) * (m("SF_MuIsoTightGH",mInt("lepInd1")) - m("SF_MuIsoTightGH_err",mInt("lepInd1"))) ) *  ( (20.1/36.4) * (m("SF_MuTriggerBCDEF",mInt("lepInd1")) - m("SF_MuTriggerBCDEF_err",mInt("lepInd1")) )+ (16.3/36.4) * (m("SF_MuTriggerGH",mInt("lepInd1"))) - m("SF_MuTriggerGH_err",mInt("lepInd1")))  ;
                *f["Lep_SFDown"] = m("Lep_SFDown") / m("Lep_SF");
            } else if(m("dataYear") == 2017) {
                // keeping 2016 tracking efficiency but not its error
                *f["Lep_SF"] = m("SF_SingleMuTrigger",mInt("lepInd1")) * m("SF_SingleMuIso",mInt("lepInd1")) * m("SF_SingleMuId",mInt("lepInd1")) * ( (20.1/36.4) * m("SF_MuTrackerBCDEF",mInt("lepInd1")) + (16.3/36.4) * m("SF_MuTrackerGH",mInt("lepInd1")));

                *f["Lep_SFUp"] = (m("SF_SingleMuTrigger",mInt("lepInd1")) + m("SF_SingleMuTrigger_err",mInt("lepInd1")) )* (m("SF_SingleMuIso",mInt("lepInd1")) + m("SF_SingleMuIso_err",mInt("lepInd1")) )* (m("SF_SingleMuId",mInt("lepInd1")) + m("SF_SingleMuId_err",mInt("lepInd1")) );
                *f["Lep_SFUp"] = m("Lep_SFUp") / m("Lep_SF");
                *f["Lep_SFDown"] = (m("SF_SingleMuTrigger",mInt("lepInd1")) - m("SF_SingleMuTrigger_err",mInt("lepInd1")) )* (m("SF_SingleMuIso",mInt("lepInd1")) - m("SF_SingleMuIso_err",mInt("lepInd1")) ) * (m("SF_SingleMuId",mInt("lepInd1")) - m("SF_SingleMuId_err",mInt("lepInd1")) ) ;
                *f["Lep_SFDown"] = m("Lep_SFDown") / m("Lep_SF");
            } else if(m("dataYear") == 2018) {
                *f["Lep_SF"] = m("SF_SingleMuTrigger",mInt("lepInd1")) * m("SF_SingleMuIso",mInt("lepInd1")) * m("SF_SingleMuId",mInt("lepInd1"));

                *f["Lep_SFUp"] = (m("SF_SingleMuTrigger",mInt("lepInd1")) + m("SF_SingleMuTrigger_err",mInt("lepInd1")) )* (m("SF_SingleMuIso",mInt("lepInd1")) + m("SF_SingleMuIso_err",mInt("lepInd1")) )* (m("SF_SingleMuId",mInt("lepInd1")) + m("SF_SingleMuId_err",mInt("lepInd1")) );
                *f["Lep_SFUp"] = m("Lep_SFUp") / m("Lep_SF");
                *f["Lep_SFDown"] = (m("SF_SingleMuTrigger",mInt("lepInd1")) - m("SF_SingleMuTrigger_err",mInt("lepInd1")) )* (m("SF_SingleMuIso",mInt("lepInd1")) - m("SF_SingleMuIso_err",mInt("lepInd1")) ) * (m("SF_SingleMuId",mInt("lepInd1")) - m("SF_SingleMuId_err",mInt("lepInd1")) ) ;
                *f["Lep_SFDown"] = m("Lep_SFDown") / m("Lep_SF");
            }
        } else if (mInt("isWenu") == 1) {
           if (m("do2015") == 1) {
               // used for 2015 analysis
               *f["Lep_SF"] = m("selLeptons_SF_IsoTight",mInt("lepInd1")) * m("selLeptons_SF_IdMVATight",mInt("lepInd1")) * m("SF_HLT_Ele23_WPLoose",mInt("lepInd1"));
           } else if(m("dataYear") == 2016 || m("dataYear") == 2017){
               *f["Lep_SF"] = m("SF_SingleElTrigger",mInt("lepInd1")) * m("SF_SingleElIdIso",mInt("lepInd1")) *  m("SF_egammaEffi_tracker",mInt("lepInd1"));

               *f["Lep_SFUp"] = (m("SF_SingleElTrigger",mInt("lepInd1")) + m("SF_SingleElTrigger_err",mInt("lepInd1")) )* (m("SF_SingleElIdIso",mInt("lepInd1")) + m("SF_SingleElIdIso_err",mInt("lepInd1")) ) *  (m("SF_egammaEffi_tracker",mInt("lepInd1")) + m("SF_egammaEffi_tracker_err",mInt("lepInd1")) );
               *f["Lep_SFUp"] = m("Lep_SFUp") / m("Lep_SF");
               *f["Lep_SFDown"] = (m("SF_SingleElTrigger",mInt("lepInd1")) - m("SF_SingleElTrigger_err",mInt("lepInd1")) )* (m("SF_SingleElIdIso",mInt("lepInd1")) - m("SF_SingleElIdIso_err",mInt("lepInd1")) ) *  (m("SF_egammaEffi_tracker",mInt("lepInd1")) - m("SF_egammaEffi_tracker_err",mInt("lepInd1")) );
                *f["Lep_SFDown"] = m("Lep_SFDown") / m("Lep_SF");

               // failing SF to correct for data/MC differences in lepton veto yields
               *f["LepFail_SF"] = 1.0;
               *f["LepFail_SFUp"] = 1.0;
               *f["LepFail_SFDown"] = 1.0;
               if (mInt("nElectron") > 1) {
                   for (int i=0; i < mInt("nElectron"); i++) {
                       if (i == mInt("lepInd1")) continue;
                       *f["LepFail_SF"] = m("SF_EleVeto",i);
                       *f["LepFail_SFUp"] = m("SF_EleVeto",i) + m("SF_EleVeto_err",i);
                       *f["LepFail_SFUp"] = m("LepFail_SFUp") / m("LepFail_SF");
                       *f["LepFail_SFDown"] = m("SF_EleVeto",i) - m("SF_EleVeto_err",i);
                       *f["LepFail_SFDown"] = m("LepFail_SFDown") / m("LepFail_SF");
                       break; // apply for highest pT lepton besides selected one
                   }
               }
           }
        }

        /* GET ELECTROWEAK CORRECTION */
        *f["weight_ptEWK"] = 1.0;
        int nGenT=0;
        int nGenV=0;
        double GenVpt=0;
        int Genpdgid=0;
        if(mInt("sampleIndex")!=0){
            if(int(m("reRunGenInfo"))==1){
                nGenT = mInt("nGenTop_re");
                nGenV = mInt("nGenVbosons_re");
                GenVpt = m("LeadGenVBoson_pt_re");
                Genpdgid = m("LeadGenVBoson_pdgId_re");
            } else {
                nGenT = mInt("nGenTop");
                nGenV = mInt("nGenVBosons");
                GenVpt = m("LeadGenVBoson_pt");
                Genpdgid = m("LeadGenVBoson_pdgId");
            }
        }
        if(mInt("sampleIndex")<0){
            if( nGenV == 1 ){
                TH1D* thisHist = NULL;
                if( mInt("sampleIndex") == -12500 ){
                    thisHist=ewkCorrHist_wp ;
                } else if( mInt("sampleIndex") == -12501 ){
                    thisHist=ewkCorrHist_wm ;
                } else if( mInt("sampleIndex") == -12502 ){
                    thisHist=ewkCorrHist_zll ;
                } else if( mInt("sampleIndex") == -12504 ){
                    thisHist=ewkCorrHist_znn ;
                }
                if(thisHist!=NULL){
                    *f["weight_ptEWK"] = GetVHEWKCorrFactor( m("V_pt"), thisHist );
                }
                if(debug>1000) std::cout<<"weight_ptEWK V_pt "<< m("V_pt")<<" "<<*f["weight_ptEWK"]<<std::endl;
            }
        } else if(mInt("sampleIndex")>0 && nGenT==0 && nGenV>0){
            *f["weight_ptEWK"]=ptWeightEWK(nGenV, GenVpt, m("Vtype"), Genpdgid);
        }

        /* ELECTROWEAK CORRECTION APPLIED*/

        if (m("doCutFlow")>0 && m("cutFlow")<5) {
            // lepton scale factor calculation will break for some events in the cutflow before lepton selection
            *f["weight"] = m("weight") * m("weight_PU") * m("bTagWeight") * m("weight_ptEWK") * m("weight_mettrigSF");
        } else {
            *f["weight"] = m("weight") * m("weight_PU") * m("bTagWeight") * m("weight_ptEWK") * m("Lep_SF") * m("weight_mettrigSF");
        }

        // 2016 pT(W) re-weighting from fit to data
        *f["recoWReWeight"] = 1.0;
        *f["recoWReWeightUp"] = 1.0;
        *f["recoWReWeightDown"] = 1.0;
        *f["recoWReWeight"] = getVPtCorrFactor(m("V_pt"), m("sampleIndex"), 0);
        *f["recoWReWeightUp"] = getVPtCorrFactor(m("V_pt"), m("sampleIndex"), 1);
        *f["recoWReWeightDown"] = getVPtCorrFactor(m("V_pt"), m("sampleIndex"), -1);
        *f["weight"] = m("weight") * m("recoWReWeight");

        // Add NLO to LO W+jet re-weighting from Z(ll)H(bb)
        // FIXME need a flag for LO vs NLO
        int sampleIndex = mInt("sampleIndex");
        float WJetNLOWeight = 1.0;
        if ((sampleIndex<4000 || sampleIndex>4702) && (sampleIndex<5000||sampleIndex>5402) && (sampleIndex<11000 || sampleIndex>11702) && (sampleIndex<15000 || sampleIndex>15602)) { WJetNLOWeight = 1.0; }
        else if(m("twoResolvedJets")) {
            float deta_bb = fabs(m("Jet_eta",mInt("hJetInd1")) - m("Jet_eta",mInt("hJetInd2")));
            if(sampleIndex==4000 || sampleIndex==4100 || sampleIndex==4200 || sampleIndex==4300 || sampleIndex==4400 || sampleIndex==4500 || sampleIndex==4600 || sampleIndex==4700 || sampleIndex==5000 || sampleIndex==5100 || sampleIndex==5300 || sampleIndex==5400
               || sampleIndex==11000 || sampleIndex==11100 || sampleIndex == 11200 || sampleIndex == 11300 ||sampleIndex==11400 || sampleIndex==11500 || sampleIndex==11600 ||sampleIndex==11700
               || sampleIndex==15000 || sampleIndex==15100 || sampleIndex == 15200 || sampleIndex == 15300 || sampleIndex == 15400 || sampleIndex==15500|| sampleIndex==15600) {
                WJetNLOWeight = LOtoNLOWeightBjetSplitEtabb(deta_bb, 0);
                WJetNLOWeight = WJetNLOWeight*1.153;
            }
            else if(sampleIndex==4001 || sampleIndex==4101 || sampleIndex==4201 || sampleIndex==4301 || sampleIndex==4401 || sampleIndex==4501 || sampleIndex==4601 || sampleIndex==4701 || sampleIndex==5001 || sampleIndex==5101 || sampleIndex==5301 || sampleIndex==5401
               || sampleIndex==11001 || sampleIndex==11101 || sampleIndex==11201 || sampleIndex==11301 || sampleIndex==11401 || sampleIndex==11501 || sampleIndex==11601 ||sampleIndex==11701
               || sampleIndex==15001 || sampleIndex==15101 || sampleIndex==15201 || sampleIndex==15301 || sampleIndex==15401 || sampleIndex==15501|| sampleIndex==15601) {
                WJetNLOWeight = LOtoNLOWeightBjetSplitEtabb(deta_bb, 1);
                WJetNLOWeight = WJetNLOWeight*1.153;
            }
            else if(sampleIndex==4002 || sampleIndex==4102 || sampleIndex==4202 || sampleIndex==4302 || sampleIndex==4402 || sampleIndex==4502 || sampleIndex==4602 || sampleIndex==4702 || sampleIndex==5002 || sampleIndex==5102 || sampleIndex==5302 || sampleIndex==5402
               || sampleIndex==11002 || sampleIndex==11102 || sampleIndex==11202 || sampleIndex==11302 || sampleIndex==11402 || sampleIndex==11502 || sampleIndex==11602 ||sampleIndex==11702
               || sampleIndex==15002 || sampleIndex==15102 || sampleIndex==15202 || sampleIndex==15302 || sampleIndex==15402 || sampleIndex==15502|| sampleIndex==15602) {
                WJetNLOWeight = LOtoNLOWeightBjetSplitEtabb(deta_bb, 2);
                WJetNLOWeight = WJetNLOWeight*1.153;
            }
        }
        *f["weight"] = m("weight") * WJetNLOWeight;
        *f["WJetNLOWeight"] = WJetNLOWeight;
        if (cursyst->name != "nominal") {
            *f[Form("weight_%s", cursyst->name.c_str())] = m("weight");
        }
    }

    if (cursyst->name != "nominal") {
        if(debug>100) std::cout<<"setting up branches for systematic"<<cursyst->name<<std::endl;
        *in[Form("controlSample_%s", cursyst->name.c_str())]                     = mInt("controlSample");

        *f[ Form("H_mass_%s",                            cursyst->name.c_str())] = m("H_mass");
        *f[ Form("H_pt_%s",                              cursyst->name.c_str())] = m("H_pt");
        *f[ Form("V_pt_%s",                              cursyst->name.c_str())] = m("V_pt");
        *in[Form("nAddJets252p9_puid_%s",                cursyst->name.c_str())] = mInt("nAddJets252p9_puid");

        *f[ Form("HJ1_HJ2_dEta_%s",                      cursyst->name.c_str())] = m("HJ1_HJ2_dEta");
        *f[ Form("HJ1_HJ2_dPhi_%s",                      cursyst->name.c_str())] = m("HJ1_HJ2_dPhi");
        *f[ Form("HJ1_HJ2_dR_%s",                        cursyst->name.c_str())] = m("HJ1_HJ2_dR");
        *f[ Form("HVdPhi_%s",                            cursyst->name.c_str())] = m("HVdPhi");
        *f[ Form("HVdR_%s",                              cursyst->name.c_str())] = m("HVdR");
        *f[ Form("MET_pt_%s",                            cursyst->name.c_str())] = m("MET_pt");
        *f[ Form("Top1_mass_fromLepton_regPT_w4MET_%s",  cursyst->name.c_str())] = m("Top1_mass_fromLepton_regPT_w4MET");
        *f[ Form("V_mt_%s",                              cursyst->name.c_str())] = m("V_mt");
        *f[ Form("hJets_btagged_0_%s",                   cursyst->name.c_str())] = m("hJets_btagged_0");
        *f[ Form("hJets_btagged_1_%s",                   cursyst->name.c_str())] = m("hJets_btagged_1");
        *f[ Form("hJets_btagWP_0_%s",                    cursyst->name.c_str())] = m("hJets_btagWP_0");
        *f[ Form("hJets_btagWP_1_%s",                    cursyst->name.c_str())] = m("hJets_btagWP_1");
        *f[ Form("hJets_leadingPt_%s",                   cursyst->name.c_str())] = m("hJets_leadingPt");
        *f[ Form("hJets_pt_0_%s",                        cursyst->name.c_str())] = m("hJets_pt_0");
        *f[ Form("hJets_pt_1_%s",                        cursyst->name.c_str())] = m("hJets_pt_1");
        *f[ Form("hJets_subleadingPt_%s",                cursyst->name.c_str())] = m("hJets_subleadingPt");
        *f[ Form("jjVPtRatio_%s",                        cursyst->name.c_str())] = m("jjVPtRatio");
        *f[ Form("lepMetDPhi_%s",                        cursyst->name.c_str())] = m("lepMetDPhi");
        *f[ Form("minDPhiFromOtherJets_%s",              cursyst->name.c_str())] = m("minDPhiFromOtherJets");
        *f[ Form("nAddJet_f_%s",                         cursyst->name.c_str())] = m("nAddJet_f");
        *f[ Form("nAddJets302p5_puid_%s",                cursyst->name.c_str())] = m("nAddJets302p5_puid");
        *f[ Form("nAddJets252p5_puid_%s",                cursyst->name.c_str())] = m("nAddJets252p5_puid");
        *f[ Form("nAddJets_2lep_%s",                     cursyst->name.c_str())] = m("nAddJets_2lep");
        *f[ Form("nJets25_dR06_%s",                      cursyst->name.c_str())] = m("nJets25_dR06");
        *f[ Form("nJets30_0lep_%s",                      cursyst->name.c_str())] = m("nJets30_0lep");
        *f[ Form("nJets30_2lep_%s",                      cursyst->name.c_str())] = m("nJets30_2lep");
        *f[ Form("nLooseBtagsDR0p8_%s",                  cursyst->name.c_str())] = m("nLooseBtagsDR0p8");
        *f[ Form("nLooseBtagsDR1p0_%s",                  cursyst->name.c_str())] = m("nLooseBtagsDR1p0");
        *f[ Form("otherJetsBestBtag_%s",                 cursyst->name.c_str())] = m("otherJetsBestBtag");
        // *f[ Form("otherJetsHighestP_%s",                 cursyst->name.c_str())] = m("otherJetsHighestP");
        *in[Form("hJetInd1_%s",                          cursyst->name.c_str())] = mInt("hJetInd1");
        *in[Form("hJetInd2_%s",                          cursyst->name.c_str())] = mInt("hJetInd2");

        // *f[ Form("FatJetCand_Msoftdrop_corrected_%s",    cursyst->name.c_str())] = m("FatJetCand_Msoftdrop_corrected");
        // *f[ Form("FatJetCand_doubleB_%s",                cursyst->name.c_str())] = m("FatJetCand_doubleB");
        // *f[ Form("FatJetCand_pt_%s",                     cursyst->name.c_str())] = m("FatJetCand_pt");
        // *f[ Form("FatJetCand_tau21_%s",                  cursyst->name.c_str())] = m("FatJetCand_tau21");
        // *f[ Form("FatJetCand_tau32_%s",                  cursyst->name.c_str())] = m("FatJetCand_tau32");
        // *f[ Form("FatJetVPtRatio_%s",                    cursyst->name.c_str())] = m("FatJetVPtRatio");
        // *f[ Form("FatJetVdPhi_%s",                       cursyst->name.c_str())] = m("FatJetVdPhi");
        // *in[Form("nFatJets200_%s",                       cursyst->name.c_str())] = mInt("nFatJets200");
    }


    // FIXME nominal must be last
    if(cursyst->name=="nominal"){
        tempfile->cd();
        //ofile->cd();
        if (debug>10000) std::cout<<"filling output tree"<<std::endl;
        outputTree->Fill();
    }

    return;
}

void VHbbAnalysis::TermAnalysis(bool skim){
    if(debug>10) std::cout<<"START TermAnalysis()"<<std::endl;
    ofile->cd();
    if(m("slimmedOutput")>0.5 && !skim){
        outputTree->SetBranchStatus("*",0);
        for(std::map<std::string,BranchInfo*>::iterator ibranch=branchInfos.begin();
            ibranch!=branchInfos.end(); ++ibranch){
            if( !(ibranch->second->dropBranch) && ibranch->second->prov!="settings"){
                outputTree->SetBranchStatus((ibranch->first).c_str(),1);
            }
        }
    }

    outputTreeSlim = outputTree->CloneTree();
    outputTreeSlim->Write();
    ofile->Close();
    if(debug>10) std::cout<<"DONE TermAnalysis()"<<std::endl;
    return;
}

//      _|                        _|
//      _|    _|_|    _|_|_|    _|_|_|_|    _|_|    _|_|_|      _|_|_|
//      _|  _|_|_|_|  _|    _|    _|      _|    _|  _|    _|  _|_|
//      _|  _|        _|    _|    _|      _|    _|  _|    _|      _|_|
//      _|    _|_|_|  _|_|_|        _|_|    _|_|    _|    _|  _|_|_|
//                    _|
//                    _|

bool VHbbAnalysis::SelectLeptonChannel(){
    bool pass=false;
    
    // leptons have to be selected before bjets, since we need to know which channel we're in
    *in["lepInd1"] = *in["muInd1"] = *in["elInd1"] = -1;
    *in["lepInd2"] = *in["muInd2"] = *in["elInd2"] = -1;
    *in["isZnn"] = 0;
    *in["isWmunu"] = 0;
    *in["isWenu"] = 0;
    *in["isZmm"] = 0;
    *in["isZee"] = 0;

    if (debug > 1000) std::cout << "Selecting leptons" << std::endl;

    if (mInt("Vtype") == 0 && (cursample->lepFlav == -1 || cursample->lepFlav == 0)) {
        std::pair<int,int> good_muons_2lep = HighestPtGoodMuonsOppCharge(    m("muptcut_2lepchan"), m("murelisocut_2lepchan"));
        if (good_muons_2lep.second > -1) {
            *in["isZmm"] = 1;
            *in["lepInd1"] = *in["muInd1"] = good_muons_2lep.first;
            *in["lepInd2"] = *in["muInd2"] = good_muons_2lep.second;
            pass=true;
        } else {
            *in["controlSample"] = -1;
        }
    } else if (mInt("Vtype") == 1 && (cursample->lepFlav == -1 || cursample->lepFlav == 1)) {
        std::pair<int,int> good_elecs_2lep = HighestPtGoodElectronsOppCharge(m("elptcut_2lepchan"), m("elrelisocut_2lepchan"), m("elidcut_2lepchan"),0);
        if (good_elecs_2lep.second > -1) {
            *in["isZee"] = 1;
            *in["lepInd1"] = *in["elInd1"] = good_elecs_2lep.first;
            *in["lepInd2"] = *in["elInd2"] = good_elecs_2lep.second;
            pass=true;
        } else {
            *in["controlSample"] = -1;
        }
    } else if (mInt("Vtype") == 2 && (cursample->lepFlav == -1 || cursample->lepFlav == 2)) {
        std::pair<int,int> good_muons_1lep = HighestPtGoodMuonsOppCharge(    m("muptcut_1lepchan"), m("murelisocut_1lepchan"), true);
        if (good_muons_1lep.first > -1) {
            *in["isWmunu"] = 1;
            *in["eventClass"] = 0;      // TODO is this still needed????
            *in["lepInd1"] = *in["muInd1"] = good_muons_1lep.first;
            pass=true;
        } else {
            *in["controlSample"] = -1;
        }
    } else if (mInt("Vtype") == 3 && (cursample->lepFlav == -1 || cursample->lepFlav == 3)) {
        std::pair<int,int> good_elecs_1lep = HighestPtGoodElectronsOppCharge(m("elptcut_1lepchan"), m("elrelisocut_1lepchan"), m("elidcut_1lepchan"),1);
        if (good_elecs_1lep.first > -1) {
            *in["isWenu"] = 1;
            *in["eventClass"] = 1000;   // TODO is this still needed????
            *in["lepInd1"] = *in["elInd1"] = good_elecs_1lep.first;
            pass=true;
        } else {
            *in["controlSample"] = -1;
        }
    } else if (mInt("Vtype") == 4 && (cursample->lepFlav == -1 || cursample->lepFlav == 4)) {
        bool passMetFilters = false;
        if (m("dataYear") == 2016) {
            passMetFilters = (
                m("Flag_goodVertices")
                && m("Flag_globalTightHalo2016Filter")
                && m("Flag_HBHENoiseFilter")
                && m("Flag_HBHENoiseIsoFilter")
                && m("Flag_EcalDeadCellTriggerPrimitiveFilter")
            );
        } else if (m("dataYear") == 2017) {
            passMetFilters = (
                m("Flag_goodVertices")
                && m("Flag_globalTightHalo2016Filter")
                && m("Flag_HBHENoiseFilter")
                && m("Flag_HBHENoiseIsoFilter")
                && m("Flag_EcalDeadCellTriggerPrimitiveFilter")
                && m("Flag_BadPFMuonFilter")
                && m("Flag_BadChargedCandidateFilter")
                && m("Flag_ecalBadCalibFilter")
            );
        } else if (m("dataYear") == 2018){
            passMetFilters = (
                m("Flag_goodVertices")
                && m("Flag_globalSuperTightHalo2016Filter")
                && m("Flag_HBHENoiseFilter")
                && m("Flag_HBHENoiseIsoFilter")
                && m("Flag_EcalDeadCellTriggerPrimitiveFilter")
                && m("Flag_BadPFMuonFilter")
                && m("Flag_BadChargedCandidateFilter")
                && m("Flag_ecalBadCalibFilter")
            );
        }
        if (mInt("sampleIndex") == 0) { // Data
            passMetFilters = passMetFilters && m("Flag_eeBadScFilter");
        }
        if ((m("MET_Pt") > m("metcut_0lepchan")) && passMetFilters) {
            *in["isZnn"] = 1;
            pass=true;
        } else {
            *in["controlSample"] = -1;
        }
    } else {
        *in["controlSample"] = -1;
    }

    if (debug > 1000) {
        std::cout << "Vtype isZnn isWmunu isWenu isZmm isZee controlSample ";
        std::cout << mInt("Vtype") << " "
                  << mInt("isZnn") << " "
                  << mInt("isWmunu") << " "
                  << mInt("isWenu") << " "
                  << mInt("isZmm") << " "
                  << mInt("isZee") << " "
                  << mInt("controlSample")
                  << std::endl;
    }

    return pass;
}


//        _|              _|
//              _|_|    _|_|_|_|    _|_|_|
//        _|  _|_|_|_|    _|      _|_|
//        _|  _|          _|          _|_|
//        _|    _|_|_|      _|_|  _|_|_|
//        _|
//      _|


bool VHbbAnalysis::SelectJets(){
    bool pass=true;

    if (debug > 1000) std::cout << "Looping over jets" << std::endl;
    UpdateJetPts();
    for (int i = 0; i < mInt("nJet"); i++) {
        if (int(m("doCMVA")) != 0) {
           // use CMVAv2 discriminator instead of CSV
           f["Jet_btagCSVV2"][i] = m("Jet_btagCMVA",i);
        }
    }

    // leptons had to be selected before bjets, since we need to know which channel we're in
    // keep track of each jet selection method separately
    // FIXME this could be set outside of the event loop
    float j1ptCut, j2ptCut;
    int j1BtagWP=-1;
    if (mInt("isZnn")) {
        j1ptCut = m("j2ptCut_0lepchan");
        j2ptCut = m("j2ptCut_0lepchan");
        j1BtagWP = m("j1BtagWP_0lepchan");
    } else if (mInt("isWmunu") || mInt("isWenu")) {
        j1ptCut = m("j1ptCut_1lepchan");
        j2ptCut = m("j2ptCut_1lepchan");
        j1BtagWP = m("j1BtagWP_1lepchan");
    } else if (mInt("isZmm") || mInt("isZee")) {
        j1ptCut = m("j1ptCut_2lepchan");
        j2ptCut = m("j2ptCut_2lepchan");
        j1BtagWP = m("j1BtagWP_2lepchan");
    }

    if(j1BtagWP==2){
        *f["j1Btag"]=m("tagWPT");
    } else if(j1BtagWP==1){
        *f["j1Btag"]=m("tagWPM");
    } else if(j1BtagWP==0){
        *f["j1Btag"]=m("tagWPL");
    } else {
        std::cout<<"leading bTag WP not defined; setting to -1."<<std::endl;
        *f["j1Btag"]=-1;
    }

    if(m("j2BtagWP")==2){
        *f["j2Btag"]=m("tagWPT");
    } else if(m("j2BtagWP")==1){
        *f["j2Btag"]=m("tagWPM");
    } else if(m("j2BtagWP")==0){
        *f["j2Btag"]=m("tagWPL");
    } else {
        std::cout<<"subleading bTag WP not defined; setting to -1."<<std::endl;
        *f["j2Btag"]=-1;
    }

    // for the moment keep also the old function to fill the branches declared in newbranches.txt
    std::pair<int,int> bjets_bestDeepCSV = HighestDeepCSVBJets(j1ptCut, j2ptCut);
    *in["hJetInd1_bestDeepCSV"] = bjets_bestDeepCSV.first;
    *in["hJetInd2_bestDeepCSV"] = bjets_bestDeepCSV.second;
    std::pair<int,int> bjets_bestCMVA = HighestCMVABJets(j1ptCut, j2ptCut);
    *in["hJetInd1_bestCMVA"] = bjets_bestCMVA.first;
    *in["hJetInd2_bestCMVA"] = bjets_bestCMVA.second;
    std::pair<int,int> bjets_bestCSV = HighestCSVBJets(j1ptCut, j2ptCut);
    *in["hJetInd1_bestCSV"] = bjets_bestCSV.first;
    *in["hJetInd2_bestCSV"] = bjets_bestCSV.second;
    std::pair<int,int> bjets_highestPt = HighestPtBJets();
    *in["hJetInd1_highestPt"] = bjets_highestPt.first;
    *in["hJetInd2_highestPt"] = bjets_highestPt.second;
    std::pair<int,int> bjets_highestPtJJ = HighestPtJJBJets();
    *in["hJetInd1_highestPtJJ"] = bjets_highestPtJJ.first;
    *in["hJetInd2_highestPtJJ"] = bjets_highestPtJJ.second;

    // the jet selection algorithm we actually use for the rest of the analysis chain
    std::pair<int, int> bjets_bestTagger =  HighestTaggerValueBJets(j1ptCut, j2ptCut, taggerName);

   // put B-Tagger cuts out of selection functions
    if (bjets_bestTagger.first != -1 && bjets_bestTagger.second != -1) {
        *b["twoResolvedJets"]=true; // jets don't need to be b-tagged to be "resolved"
        *in["hJetInd1"] = bjets_bestTagger.first;
        *in["hJetInd2"] = bjets_bestTagger.second;
        if (mInt("isZnn")) { // Higgs jet pT cuts are applied to the pair sorted by pT instead of b-tag value
            float j1pt = m("Jet_bReg",bjets_bestTagger.first);
            float j2pt = m("Jet_bReg",bjets_bestTagger.second);
            if (std::max(j1pt, j2pt) < m("j1ptCut_0lepchan")) {
                *in["controlSample"] = -1;
            }
        }
    }

    if (int(m("doBoost")) != 0) {
        // need to do some filtering on b-tagging for events that we won't use in the boosted analysis
        // (i.e. there is no fat jet) or file size for W+jets is ridiculous for boosted analysis
        FatJetSelection();

        if(mInt("FatJetCand_index")>-1){
            *b["oneMergedJet"]=true;
        }
    } 
    
    if (!m("twoResolvedJets") && !m("oneMergedJet")) {
        pass=false;
    }
    return pass;
}




//      _|        _|
//      _|_|_|          _|_|_|    _|_|_|    _|_|_|      _|  _|_|    _|_|      _|_|_|    _|_|
//      _|    _|  _|  _|    _|  _|    _|  _|_|          _|_|      _|_|_|_|  _|        _|    _|
//      _|    _|  _|  _|    _|  _|    _|      _|_|      _|        _|        _|        _|    _|
//      _|    _|  _|    _|_|_|    _|_|_|  _|_|_|        _|          _|_|_|    _|_|_|    _|_|
//                          _|        _|
//                      _|_|      _|_|


bool VHbbAnalysis::ReconstructHiggsCand(){
    bool pass=true;
    if(m("twoResolvedJets")){
        if (debug > 1000) {
            std::cout << "nJet = " << mInt("nJet") << std::endl;
            std::cout << "hJetInd1 = " << mInt("hJetInd1") << std::endl;
            std::cout << "hJetInd2 = " << mInt("hJetInd2") << std::endl;
            std::cout << "found two bjets with pt and "<< taggerName << ": "
                      << m("Jet_bReg",mInt("hJetInd1")) << " "
                      << m(taggerName,mInt("hJetInd1")) << " "
                      << m("Jet_bReg",mInt("hJetInd2")) << " "
                         << m(taggerName,mInt("hJetInd2")) << " "
                      << std::endl;
        }

        // Reconstruct Higgs
        HJ1.SetPtEtaPhiM(m("Jet_bReg",mInt("hJetInd1")), m("Jet_eta",mInt("hJetInd1")), m("Jet_phi",mInt("hJetInd1")), m("Jet_mass",mInt("hJetInd1")) * (m("Jet_bReg",mInt("hJetInd1")) / m("Jet_Pt",mInt("hJetInd1"))));
        HJ2.SetPtEtaPhiM(m("Jet_bReg",mInt("hJetInd2")), m("Jet_eta",mInt("hJetInd2")), m("Jet_phi",mInt("hJetInd2")), m("Jet_mass",mInt("hJetInd2")) * (m("Jet_bReg",mInt("hJetInd2")) / m("Jet_Pt",mInt("hJetInd2"))));
        Hbb = HJ1 + HJ2;

        HJ1_noreg.SetPtEtaPhiM(m("Jet_Pt",mInt("hJetInd1")), m("Jet_eta",mInt("hJetInd1")), m("Jet_phi",mInt("hJetInd1")), m("Jet_mass",mInt("hJetInd1")));
        HJ2_noreg.SetPtEtaPhiM(m("Jet_Pt",mInt("hJetInd2")), m("Jet_eta",mInt("hJetInd2")), m("Jet_phi",mInt("hJetInd2")), m("Jet_mass",mInt("hJetInd2")));
        Hbb_noreg = HJ1_noreg + HJ2_noreg;


        // CP ADD BACK
        if(debug>1000) std::cout<<"doOnlySignalRegion controlSample "<<m("doOnlySignalRegion")<<" "<<mInt("controlSample")<<std::endl;
        if (m("doOnlySignalRegion")>0 && mInt("controlSample") < 0) {
            pass=false;
        }

        //*f["H_mass_step2"] = m("H_mass"); // I don't think we need this.
        *f["H_mass_noreg"] = Hbb_noreg.M();

        //Recovering ISR/FSR
        bool enableFSRRecovery = int(m("enableFSRRecovery")) > 0;
        *b["recoFSR"] = false;
        int nFSRJets = 0;

        if( enableFSRRecovery ){
            *b["recoFSR"] = true;
            *f["H_mass_noFSR"] = Hbb.M();
            *f["H_pt_noFSR"] = Hbb.Pt();
            // di-jet kinematics

            *f["HJ1_HJ2_dPhi_noFSR"] = HJ1.DeltaPhi(HJ2);
            *f["HJ1_HJ2_dEta_noFSR"] = fabs(HJ1.Eta() - HJ2.Eta());
            *f["HJ1_HJ2_dR_noFSR"] = HJ1.DeltaR(HJ2);
            *f["JJEtaBal_noFSR"] = (fabs(m("Jet_eta",mInt("hJetInd1")) + m("Jet_eta",mInt("hJetInd2")))) / (fabs(m("Jet_eta",mInt("hJetInd1")) - m("Jet_eta",mInt("hJetInd2"))));

            for(int ijet = 0; ijet < mInt("nJet"); ijet++) {
                if( ijet == mInt("hJetInd1") || ijet == mInt("hJetInd2") ) continue;

                //Select FSR jets
                if( m("Jet_Pt", ijet) > 20 && abs(m("Jet_eta", ijet)) < 3.0 && (m("Jet_puId", ijet) > 6 || m("Jet_Pt",ijet)>50)  && m("Jet_lepFilter", ijet) > 0 ){

                    TLorentzVector Jet_FSR;
                    Jet_FSR.SetPtEtaPhiM( m("Jet_bReg",ijet), m("Jet_eta",ijet), m("Jet_phi",ijet), m("Jet_mass",ijet) * (m("Jet_bReg",ijet) / m("Jet_Pt",ijet)) );

                    if( std::min( Jet_FSR.DeltaR( HJ1 ), Jet_FSR.DeltaR( HJ2 ) ) < 0.8){
                        if(Jet_FSR.DeltaR( HJ1 )<Jet_FSR.DeltaR( HJ2 )){
                            //add FSR to HJ1
                            nFSRJets++;
                            HJ1 += Jet_FSR;
                        }else if(Jet_FSR.DeltaR( HJ2 )<Jet_FSR.DeltaR( HJ1 )){
                            //add FSR to HJ2
                            nFSRJets++;
                            HJ2 += Jet_FSR;
                        }
                    }
                }
            }

            //recompute the total Higgs 4-Momentum
            Hbb = ( HJ1 + HJ2 );

        }// end of FSR recovering

        // di-jet kinematics
        *f["HJ1_pt"] = HJ1.Pt();
        *f["HJ2_pt"] = HJ2.Pt();
        *f["HJ1_HJ2_dPhi"] = HJ1.DeltaPhi(HJ2);
        *f["HJ1_HJ2_dEta"] = fabs(HJ1.Eta() - HJ2.Eta());
        *f["HJ1_HJ2_dR"] = HJ1.DeltaR(HJ2);
        *f["JJEtaBal"] = ( fabs( HJ1.Eta() + HJ2.Eta() ) / fabs( HJ1.Eta() - HJ2.Eta() ) );

        // Now we can calculate whatever we want (transverse) with V and H four-vectors
        *in["nFSRJets"] = nFSRJets;
        *f["H_mass"] = Hbb.M(); // mass window cut? regression applied in FinishEvent
        *f["H_pt"] = Hbb.Pt();
        
        if (cursyst->name != "nominal") {
            *f[Form("H_mass_%s", cursyst->name.c_str())] = Hbb.M();
            *f[Form("H_pt_%s", cursyst->name.c_str())] = Hbb.Pt();
        }
        if (mInt("isZnn") && m("H_pt") < m("hptcut_0lepchan")) {
            if (int(m("doBoost")) == 0 || mInt("FatJetCand_index")==-1) {//FIXME
                *in["controlSample"] = -1;
            }
        } else if ((m("isWmunu") || m("isWenu")) && m("H_pt") < m("hptcut_1lepchan")) {
            if (int(m("doBoost")) == 0 || mInt("FatJetCand_index")==-1) {//FIXME
                *in["controlSample"] = -1;
            }
        } else if (mInt("controlSample") > -1) {
            *in["cutFlow"] += 1; // pT(jj) cut
        }
    }
    if(*b["oneMergedJet"]){
        //no construction necessary... the fat jet is already selected.
        pass=true;
    }

    return pass;
}


//      _|      _|
//      _|      _|      _|  _|_|    _|_|      _|_|_|    _|_|
//      _|      _|      _|_|      _|_|_|_|  _|        _|    _|
//        _|  _|        _|        _|        _|        _|    _|
//          _|          _|          _|_|_|    _|_|_|    _|_|

bool VHbbAnalysis::ReconstructVCand(){
    bool pass=true;

    // channel specific selection and vector boson reconstruction
    *f["Top1_mass_fromLepton_regPT"] = -999;
    *f["Top1_mass_fromLepton_regPT_wMET"] = -999;
    *f["Top1_mass_fromLepton_regPT_w4MET"] = -999;
    *f["Top1_mass_fromLepton_smearedPT"] = -999;
    *f["Top1_mass_fromLepton_smearedPT_wMET"] = -999;
    *f["Top1_mass_fromLepton_smearedPT_w4MET"] = -999;
    *f["Top1_mass_fromLepton_unsmearedPT"] = -999;
    *f["Top1_mass_fromLepton_unsmearedPT_wMET"] = -999;
    *f["Top1_mass_fromLepton_unsmearedPT_w4MET"] =  -999;


    if (mInt("isZee") == 1 || mInt("isZmm") == 1) {
        TLorentzVector lep1, lep2;
        if (mInt("isZee") == 1){
            lep1.SetPtEtaPhiM(m("Electron_pt",mInt("lepInd1")), m("Electron_eta",mInt("lepInd1")), m("Electron_phi",mInt("lepInd1")), m("Electron_mass",mInt("lepInd1")));
            lep2.SetPtEtaPhiM(m("Electron_pt",mInt("lepInd2")), m("Electron_eta",mInt("lepInd2")), m("Electron_phi",mInt("lepInd2")), m("Electron_mass",mInt("lepInd2")));
        } else if (mInt("isZmm") == 1){
            lep1.SetPtEtaPhiM(m("Muon_pt",mInt("lepInd1")), m("Muon_eta",mInt("lepInd1")), m("Muon_phi",mInt("lepInd1")), m("Muon_mass",mInt("lepInd1")));
            lep2.SetPtEtaPhiM(m("Muon_pt",mInt("lepInd2")), m("Muon_eta",mInt("lepInd2")), m("Muon_phi",mInt("lepInd2")), m("Muon_mass",mInt("lepInd2")));
        } else {
            std::cout<<"This shouldn't happen.  isZee  isZmm "<<mInt("isZee")<<" "<<mInt("isZmm")<<std::endl;
        }
        V = lep1 + lep2;
    } else if (mInt("isWenu") == 1 || mInt("isWmunu") == 1) {
        if (mInt("isWenu") == 1){
            *f["selLeptons_pt_0"] = m("Electron_pt",mInt("lepInd1"));
            *f["selLeptons_eta_0"] = m("Electron_eta",mInt("lepInd1"));
            *f["selLeptons_phi_0"] = m("Electron_phi",mInt("lepInd1"));
            *f["selLeptons_mass_0"] = m("Electron_mass",mInt("lepInd1"));
            *in["selLeptons_pdgId_0"] = mInt("Electron_charge",mInt("lepInd1"))*11;
        } else if (mInt("isWmunu") == 1){
            *f["selLeptons_pt_0"] = m("Muon_pt",mInt("lepInd1"));
            *f["selLeptons_eta_0"] = m("Muon_eta",mInt("lepInd1"));
            *f["selLeptons_phi_0"] = m("Muon_phi",mInt("lepInd1"));
            *f["selLeptons_mass_0"] = m("Muon_mass",mInt("lepInd1"));
            *in["selLeptons_pdgId_0"] = m("Muon_charge",mInt("lepInd1"))*13;
        } else {
            std::cout<<"This shouldn't happen.  isWenu  isWmunu "<<mInt("isWenu")<<" "<<mInt("isWmunu")<<std::endl;
        }

        //FIXME why do we do this?  Cut-flow?  CP
        if (mInt("lepInd1") == -1) {
            // not Wenu or Wmunu, use preselected lepton
            *in["lepInd1"] = 0;
        }

        if (debug > 1000) {
            std::cout << "cutting on dphi(lep, met)" << std::endl;
            std::cout << mInt("lepInd1") << " "
            //          << f["selLeptons_phi"][*in["lepInd1"]] << " "
            //          << f["selLeptons_eta"][*in["lepInd1"]] << " "
            //          << f["selLeptons_mass"][*in["lepInd1"]]
                      << std::endl;
        }

        *f["lepMetDPhi"] = fabs(EvalDeltaPhi(m("selLeptons_phi_0"), m("MET_Phi")));

        if (mInt("isWmunu") && m("lepMetDPhi") > m("muMetDPhiCut")) {
            *in["controlSample"] = -1;
        } else if (mInt("isWenu") && m("lepMetDPhi") > m("elMetDPhiCut")) {
            *in["controlSample"] = -1;
        }

        // CP ADD BACK
        if (m("doOnlySignalRegion")>0 && mInt("controlSample") < 0) {
            pass=false;
        }

        if (mInt("controlSample") > -1) {
            *in["cutFlow"] += 1;
        }

        TLorentzVector Lep, MET;
        // Reconstruct W
        if (debug > 1000) {
            std::cout << "met " << m("MET_Pt") << " " << m("MET_Phi") << std::endl;
        }

        MET.SetPtEtaPhiM(m("MET_Pt"), 0., m("MET_Phi"), 0.); // Eta/M don't affect calculation of W.pt and W.phi
        Lep.SetPtEtaPhiM(m("selLeptons_pt_0"), m("selLeptons_eta_0"), m("selLeptons_phi_0"), m("selLeptons_mass_0"));
        double cosPhi12 = (Lep.Px()*MET.Px() + Lep.Py()*MET.Py()) / (Lep.Pt() * MET.Pt()); // cos of the angle between the lepton and the missing energy
        *f["V_mt"] = TMath::Sqrt(2*Lep.Pt()*MET.Pt() * (1 - cosPhi12));
        if(mInt("hJetInd1")>-1 && mInt("hJetInd2")>-1) {
            *f["Lep_HJ1_dPhi"] = Lep.DeltaPhi(HJ1);
            *f["Lep_HJ2_dPhi"] = Lep.DeltaPhi(HJ2);
        }
        V = MET + Lep;

        TLorentzVector neutrino = getNu4Momentum(Lep, MET);
        W_withNuFromMWCon = neutrino + Lep;

        //FIXME we should reduce this code to keep only what we use

        *f["Top1_mass_fromLepton_regPT"] = GetRecoTopMass(Lep, false, 0, true, true); // construct top mass from closest jet to lepton
        *f["Top1_mass_fromLepton_regPT_wMET"] = GetRecoTopMass(Lep, false, 1, true, true); // construct top mass from closest jet to lepton
        *f["Top1_mass_fromLepton_regPT_w4MET"] = GetRecoTopMass(Lep, false, 2, true, true); // construct top mass from closest jet to lepton
        *f["Top1_mass_fromLepton_smearedPT"] = GetRecoTopMass(Lep, false, 0, false, true); // construct top mass from closest jet to lepton
        *f["Top1_mass_fromLepton_smearedPT_wMET"] = GetRecoTopMass(Lep, false, 1, false, true); // construct top mass from closest jet to lepton
        *f["Top1_mass_fromLepton_smearedPT_w4MET"] = GetRecoTopMass(Lep, false, 2, false, true); // construct top mass from closest jet to lepton
        *f["Top1_mass_fromLepton_unsmearedPT"] = GetRecoTopMass(Lep, false, 0, false, false); // construct top mass from closest jet to lepton
        *f["Top1_mass_fromLepton_unsmearedPT_wMET"] = GetRecoTopMass(Lep, false, 1, false, false); // construct top mass from closest jet to lepton
        *f["Top1_mass_fromLepton_unsmearedPT_w4MET"] = GetRecoTopMass(Lep, false, 2, false, false); // construct top mass from closest jet to lepton


    } else if (mInt("isZnn") == 1) {
        TLorentzVector MET;
        MET.SetPtEtaPhiM(m("MET_Pt"), 0., m("MET_Phi"), 0.); // Eta/M don't affect calculation of V.pt and V.phi
        V = MET;
    } else {
        std::cerr << "not any selection... can I go now?" << std::endl;  // haha, this is a good one.
    }

    *f["V_pt"] = V.Pt();
    if (cursyst->name != "nominal") {
        *f[Form("V_pt_%s", cursyst->name.c_str())] = V.Pt();
    }

    if (m("V_pt") < 50) {                               // V_pt is lower in the 2-lepton channel
        *in["controlSample"] = -1;
    } else if (!(mInt("isZee") == 1 || mInt("isZmm") == 1) && m("V_pt") < m("vptcut")) {
        *in["controlSample"] = -1;
    } else if (mInt("controlSample") > -1) {
        *in["cutFlow"] += 1; // pT(W) cut
    }

    if (mInt("isZee") == 1 || mInt("isZmm") == 1) {
        *f["V_mass"] = V.M();
        if (m("V_mass") < m("zmasslow")|| m("V_mass") > m("zmasshigh")) {
            *in["controlSample"] = -1;
        }
    }

    return pass;
}


//      _|      _|      _|      _|    _|      _|        _|                                                  _|      _|
//      _|      _|      _|      _|    _|      _|  _|        _|_|_|      _|_|    _|_|_|  _|_|      _|_|_|  _|_|_|_|        _|_|_|    _|_|_|
//      _|      _|  _|_|_|_|_|  _|_|_|_|      _|_|      _|  _|    _|  _|_|_|_|  _|    _|    _|  _|    _|    _|      _|  _|        _|_|
//        _|  _|        _|      _|    _|      _|  _|    _|  _|    _|  _|        _|    _|    _|  _|    _|    _|      _|  _|            _|_|
//          _|          _|      _|    _|      _|    _|  _|  _|    _|    _|_|_|  _|    _|    _|    _|_|_|      _|_|  _|    _|_|_|  _|_|_|

void VHbbAnalysis::ComputeVHKinematics(){

    // Calculate V+H kinematics

   // CP ADD BACK
   if (m("doOnlySignalRegion")>0 && mInt("controlSample") < 0) {
       return;
   }

    if(m("twoResolvedJets")){

        *f["jjVPtRatio"] = m("H_pt") / m("V_pt");
        *f["HVdPhi"] = fabs(Hbb.DeltaPhi(V));
        *f["HVdEta"] = fabs(Hbb.Eta() - V.Eta());
        *f["HVdR"]   = Hbb.DeltaR(V);

        if( enableFSRRecovery ){
            *f["jjVPtRatio_noFSR"] = m("H_pt_noFSR") / m("V_pt");
            *f["HVdPhi_noFSR"] = fabs(Hbb_noFSR.DeltaPhi(V));
            *f["HVdEta_noFSR"] = fabs(Hbb_noFSR.Eta() - V.Eta());
            *f["HVdR_noFSR"]   = Hbb_noFSR.DeltaR(V);
        }


        if (mInt("isWmunu") == 1 || mInt("isWenu") == 1) {
            *f["HVdEta_4MET"] = fabs(Hbb.Eta() -  W_withNuFromMWCon.Eta());
            if( enableFSRRecovery ){
                *f["HVdEta_4MET_noFSR"] = fabs(Hbb_noFSR.Eta() -  W_withNuFromMWCon.Eta());
            }
        }

    }// End of two resolved jets computation

    if (cursyst->name != "nominal") {
        *in[Form("nAddJets252p9_puid_%s", cursyst->name.c_str())] = mInt("nAddJets252p9_puid");
    }

    if(m("oneMergedJet")==1){
        *f["FatJetVPtRatio"] = fatJetCand.Pt() / m("V_pt");
        *f["FatJetVdPhi"] = fabs(fatJetCand.DeltaPhi(V));
        *f["FatJetVdEta"] = fabs(fatJetCand.Eta() - V.Eta());
        *f["FatJetVdR"]   = fatJetCand.DeltaR(V);

        BoostedSelection();
    }
}

//                      _|        _|  _|    _|      _|                                _|        _|            _|
//        _|_|_|    _|_|_|    _|_|_|      _|_|_|_|        _|_|    _|_|_|      _|_|_|  _|                      _|
//      _|    _|  _|    _|  _|    _|  _|    _|      _|  _|    _|  _|    _|  _|    _|  _|        _|            _|
//      _|    _|  _|    _|  _|    _|  _|    _|      _|  _|    _|  _|    _|  _|    _|  _|        _|            _|
//        _|_|_|    _|_|_|    _|_|_|  _|      _|_|  _|    _|_|    _|    _|    _|_|_|  _|        _|    _|      _|
//                                                                                              _|  _|
//                                                                                            _|

void VHbbAnalysis::ComputeOtherEventKinematics(){

    // count the number of additional leptons and jets, then cut on this number
    int nAddLep = 0;


    if (debug > 1000) std::cout << "counting additional jets and leptons" << std::endl;
    // 15 to 30 by 1 GeV, 1.5 to 3 w/ 0.1 in eta
    //std::vector<int> ptCuts = {15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35};
    //std::vector<double> etaCuts = {1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.1,3.2,3.3,3.4,3.5};
    std::vector<int> ptCuts = {25, 30};
    std::vector<double> etaCuts = {2.9, 2.5};

    for (int i = 0; i < (int) ptCuts.size(); i++) {
        for (int j = 0; j < (int) etaCuts.size(); j++) {
            float maxPt = 0; // max. pT of additional jets
            int nAddJet_tmp = 0;
            std::string eta_cut = Form("%.1f", etaCuts[j]); // convert, say, 2.5 to '2p5'
            std::replace(eta_cut.begin(), eta_cut.end(), '.', 'p');
            for (int k = 0; k < mInt("nJet"); k++) {
                if (k == mInt("hJetInd1") || k == mInt("hJetInd2")) continue;
                if (m("Jet_lepFilter",k) && m("Jet_Pt",k) > ptCuts[i] && fabs(m("Jet_eta",k)) < etaCuts[j] && (m("Jet_puId",k) > 6 || m("Jet_Pt",k)>50)) {
                    nAddJet_tmp++;
                    if (m("Jet_Pt",k) > maxPt) {
                        maxPt = m("Jet_Pt",k);
                        *f[Form("AddJets%i%s_puid_leadJet_pt", ptCuts[i], eta_cut.c_str())] = maxPt;
                        *f[Form("AddJets%i%s_puid_leadJet_eta", ptCuts[i], eta_cut.c_str())] = m("Jet_eta",k);
                        *f[Form("AddJets%i%s_puid_leadJet_phi", ptCuts[i], eta_cut.c_str())] = m("Jet_phi",k);
                        *f[Form("AddJets%i%s_puid_leadJet_btagged", ptCuts[i], eta_cut.c_str())] = m(taggerName,k);
                    }
                }
            }
            if( i==1 && j==1){
                *f[Form("nAddJets%i%s_puid", ptCuts[i], eta_cut.c_str())] = nAddJet_tmp;
            }else{
                *in[Form("nAddJets%i%s_puid", ptCuts[i], eta_cut.c_str())] = nAddJet_tmp;
            }
        }
    }

    *f["nAddJets_2lep"] = 0;
    for(int iJet=0; iJet<mInt("nJet"); iJet++){
        if(
               iJet!=mInt("hJetInd1")
            && iJet!=mInt("hJetInd2")
            && (mInt("Jet_puId",iJet)>6 || m("Jet_Pt",iJet)>6)
            && m("Jet_Pt",iJet)>30
            && mInt("Jet_lepFilter",iJet)>0
            && abs(m("Jet_eta",iJet))<2.4
        ) {
            *f["nAddJets_2lep"]=m("nAddJets_2lep")+1;
        }
    }

    // count additional leptons (check both collections, which are exclusive)
    for (int i = 0; i < mInt("nMuon"); i++) {
        if (i == mInt("muInd1")) continue; // don't look at the lepton we've selected from the W
        if (i == mInt("muInd2")) continue; // don't look at the lepton we've selected from the Z
        if (m("Muon_pt",i) > 15 && fabs(m("Muon_eta",i)) < 2.5 && m("Muon_pfRelIso04_all",i) < 0.1) {
            nAddLep++;
        }
    }

    for (int i = 0; i < mInt("nElectron"); i++) {
        if (i == mInt("elInd1")) continue; // don't look at the lepton we've selected from the W
        if (i == mInt("elInd2")) continue; // don't look at the lepton we've selected from the Z
        if (m("Electron_pt",i) > 15 && fabs(m("Electron_eta",i)) < 2.5 && m("Electron_pfRelIso03_all",i) < 0.1) {
            nAddLep++;
        }
    }


    *in["nAddLeptons"] = nAddLep;
    if ((mInt("isZnn") || mInt("isWmunu") || mInt("isWenu")) && nAddLep >= m("nAddLeptonsCut")) {
        *in["controlSample"] = -1;
    } else if (mInt("controlSample") > -1) {
        *in["cutFlow"] += 1; // additional lepton veto
    }


    if ((mInt("isZnn") || mInt("isWmunu") || mInt("isWenu")) && m("nAddJets302p5_puid") >= m("nAddJetsCut")) {
        *in["controlSample"] = -1;
    } else if (mInt("controlSample") > -1) {
        *in["cutFlow"] += 1; // additional jet veto
    }

    if (mInt("isZnn") && m("HVdPhi") < m("HVDPhiCut_0lepchan")) {
        *in["controlSample"] = -1;
    } else if ((m("isWmunu") || m("isWenu")) && m("HVdPhi") < m("HVDPhiCut_1lepchan")) {
        *in["controlSample"] = -1;
    } else if ((mInt("isZmm") || mInt("isZee")) && m("HVdPhi") < m("HVDPhiCut_2lepchan")) {
        *in["controlSample"] = -1;
    } else if ( mInt("controlSample") > -1) {
        *in["cutFlow"] += 1; // dPhi(jj,W) cut
    }

    if (m("doOnlySignalRegion")>0 && mInt("controlSample") < 0) {
        return;
    }

    // Iterate over the jet collection and count the number of central jets and
    // the number of jets azimuthally within 0.5 radians of the reconstructed
    // Z boson (MET). The latter is used to reject QCD events.
    int nJetsCloseToMET = 0;
    for (int i = 0; i < mInt("nJet"); i++) {
        if (m("Jet_Pt",i) > 30 && (mInt("Jet_puId",i) > 6 || m("Jet_Pt",i)>50) && m("Jet_lepFilter",i)) {
            float dPhi_Jet_MET = fabs(EvalDeltaPhi(m("Jet_phi",i), m("MET_Phi")));
            if (dPhi_Jet_MET < 0.5) {
                nJetsCloseToMET += 1;
            }
        }
    }
    *in["nJetsCloseToMET"] = nJetsCloseToMET;
    if (mInt("isZnn") && mInt("nJetsCloseToMET") != 0) {
        *in["controlSample"] = -1;
    }

    float dPhi_hJet1_MET = -99;
    float dPhi_hJet2_MET = -99;
    float min_dPhi_hJet_MET = -99;
    float dPhi_MET_TkMET = -99;
    if (m("twoResolvedJets")) {
        dPhi_hJet1_MET = fabs(EvalDeltaPhi(m("Jet_phi",mInt("hJetInd1")), m("MET_Phi")));
        dPhi_hJet2_MET = fabs(EvalDeltaPhi(m("Jet_phi",mInt("hJetInd2")), m("MET_Phi")));
        min_dPhi_hJet_MET = std::min(dPhi_hJet1_MET, dPhi_hJet2_MET);
        dPhi_MET_TkMET = fabs(EvalDeltaPhi(m("MET_Phi"), m("TkMET_phi")));
        *f["min_dPhi_hJet_MET"] = min_dPhi_hJet_MET;
        *f["dPhi_MET_TkMET"] = dPhi_MET_TkMET;
        if (mInt("isZnn") && m("dPhi_MET_TkMET") > m("dPhiMETTkMETCut")) {
            *in["controlSample"] = -1;
        }
    }
}


//                                      _|                          _|                                                    _|
//        _|_|_|    _|_|    _|_|_|    _|_|_|_|  _|  _|_|    _|_|    _|        _|_|_|    _|_|_|  _|_|_|  _|_|    _|_|_|    _|    _|_|      _|_|_|
//      _|        _|    _|  _|    _|    _|      _|_|      _|    _|  _|      _|_|      _|    _|  _|    _|    _|  _|    _|  _|  _|_|_|_|  _|_|
//      _|        _|    _|  _|    _|    _|      _|        _|    _|  _|          _|_|  _|    _|  _|    _|    _|  _|    _|  _|  _|            _|_|
//        _|_|_|    _|_|    _|    _|      _|_|  _|          _|_|    _|      _|_|_|      _|_|_|  _|    _|    _|  _|_|_|    _|    _|_|_|  _|_|_|
//                                                                                                              _|
//                                                                                                              _|

void VHbbAnalysis::ControlSampleSelection(){

    if (debug>1000) std::cout<<"about to apply control region selection"<<std::endl;

    // Control Samples

    if (m("twoResolvedJets")) {

        // Set aliases and precalculate quantities of interest.
        float H_mass = m("H_mass");
        float H_pt = m("H_pt");
        float V_mass = m("V_mass");
        float HVdPhi = m("HVdPhi");

        float hJet1_btag = m(taggerName,mInt("hJetInd1"));
        float hJet2_btag = m(taggerName,mInt("hJetInd2"));
        float min_hJet_btag = std::min(hJet1_btag, hJet2_btag);
        float max_hJet_btag = std::max(hJet1_btag, hJet2_btag);

        float hJet1_bReg = m("Jet_bReg",mInt("hJetInd1"));
        float hJet2_bReg = m("Jet_bReg",mInt("hJetInd2"));
        float min_hJet_bReg = std::min(hJet1_bReg, hJet2_bReg);
        float max_hJet_bReg = std::max(hJet1_bReg, hJet2_bReg);
        
        if(mInt("controlSample")!=-1){
            if (max_hJet_btag > m("j1Btag") && min_hJet_btag > m("j2Btag")) {
                *in["controlSample"] = 0;
            } else {
                *in["controlSample"] = -1;
            }
        }
            

        // 0-lepton
        bool base0LepCSSelection = (
            // Vector Boson Cuts
            (mInt("isWmunu") || mInt("isWenu") || mInt("isZnn"))
            && mInt("cutFlow") >= 2
            // Higgs Boson Cuts
            && H_mass < 500
            && H_pt > 120
            // Higgs Jet Cuts
            && max_hJet_bReg > 60
            && min_hJet_bReg > 35
            && hJet2_btag > m("tagWPL")
            // MET Filters
            && m("Flag_goodVertices")
            && m("Flag_globalTightHalo2016Filter")
            && m("Flag_HBHENoiseFilter")
            && m("Flag_HBHENoiseIsoFilter")
            && m("Flag_EcalDeadCellTriggerPrimitiveFilter")
        );

        if (m("dataYear") == 2017) {
            base0LepCSSelection = (
                base0LepCSSelection
                && m("Flag_BadPFMuonFilter")
                && m("Flag_BadChargedCandidateFilter")
                && m("Flag_ecalBadCalibFilter")
            );
        }

        if (mInt("sampleIndex") == 0) { // Data
            base0LepCSSelection = base0LepCSSelection && m("Flag_eeBadScFilter");
        }

        if (base0LepCSSelection) {
            if ((mInt("isWmunu") || mInt("isWenu")) && m("MET_Pt") > m("metcut_0lepchan")) {
                if ( *f["min_dPhi_hJet_MET"] < 1.57 && hJet1_btag > m("tagWPM") && m("nAddJets302p5_puid") >= 2 && HVdPhi > 2) {
                    *in["controlSample"] = 1; // TTbar Control Sample Index
                }
            } else if (mInt("isZnn")) {
                if (*in["nJetsCloseToMET"] == 0 && m("dPhi_MET_TkMET") < 0.5 && HVdPhi > 2) {
                if (hJet1_btag < m("tagWPM") && m("nAddJets302p5_puid") < 2) {
                        *in["controlSample"] = 2; // Z+Light Control Sample Index
                    } else if (hJet1_btag > m("tagWPT") && m("nAddJets302p5_puid") < 1 && (H_mass<60 || H_mass>160)) { //FIXME change tight->medium?
                        *in["controlSample"] = 3; // Z+bb Control Sample Index
                    }
                }
            }
        }

        // 1-lepton
        // if cutflow==0 then jets are bogus
        // jet pt
        // lepton pt, iso, id
        // MET_pt > threshold
        // deltaPhi(met,lep) -- NO
                //&& ((*in["isWmunu"] && *f["lepMetDPhi"] < *f["muMetDPhiCut"])
                //  || (*in["isWenu"] && *f["lepMetDPhi"] < *f["elMetDPhiCut"]))
        // V_pt > 100 and bb_mass<250
        bool base1LepCSSelection = (
            mInt("cutFlow") >= 2
            && (mInt("isWmunu") || mInt("isWenu"))
            && m("Jet_bReg",mInt("hJetInd1")) > m("j1ptCut_1lepchan")
            && m("Jet_bReg",mInt("hJetInd2")) > m("j2ptCut_1lepchan")
            && m("MET_Pt") > m("metcut_1lepchan")
            && mInt("nAddLeptons") == 0
            && m("V_pt") > m("vptcut")
            && H_mass < 250
            && H_pt > m("hptcut_1lepchan")
        );


        // htJet30 is not currenty in nanoAOD, need to calculate it
        *f["htJet30"] = 0.;
        for (int i=0; i<mInt("nJet");i++) {
            if (m("Jet_Pt",i)>30. && m("Jet_lepFilter",i) && (m("Jet_puId",i)>6 || m("Jet_Pt",i)>50)) {
                *f["htJet30"] = m("htJet30") + m("Jet_Pt",i);
            }
        }
        //for (int i=0; i<mInt("nMuon"); i++) {
        //    if (m("Muon_pt",i)>5 && m("Muon_pfRelIso04_all",i)<0.4) {
        //        *f["htJet30"] = m("htJet30") + m("Muon_pt",i);
        //    }
       // }
        //for (int i=0; i<mInt("nElectron"); i++) {
        //    if (m("Electron_pt",i)>5 && m("Electron_pfRelIso03_all",i)<0.4) {
        //        *f["htJet30"] = m("htJet30") + m("Electron_pt",i);
        //    }
        //}

        if (base1LepCSSelection) {
            if (max_hJet_btag > m("tagWPT")){ //ttbar or W+HF #FIXME change tight-> medium?
                if (m("nAddJets302p5_puid") > 1.5 && m("MET_Pt") < m("metcut_0lepchan")) { //ttbar, avoid overlap with Z(vv) TT CR
                    *in["controlSample"] = 11;
                //} else if (mInt("nAddJets252p9_puid") < 0.5 && m("MET_Pt")/sqrt(m("htJet30")) > 2.) { //W+HF // remove mass window so we can use the same ntuple for VV, just be careful that we always avoid overlap with SR
                } else if (m("nAddJets302p5_puid") < 1.5 && m("MET_Pt")/sqrt(m("htJet30")) > 2. && (H_mass<90 || H_mass>150)) {
                    *in["controlSample"] = 13;
                }
            }else if (max_hJet_btag > m("tagWPL") && max_hJet_btag < m("tagWPM") && m("MET_Pt")/sqrt(m("htJet30")) > 2.) { //W+LF
                *in["controlSample"] = 12;
            }

            if (mInt("sampleIndex") == 0 && debug>10) {
                std::cout << "data CS event " << mInt("controlSample")
                      << " max_hJet_btag " << max_hJet_btag
                          << " nAddJets302p5_puid " << m("nAddJets302p5_puid")
                          << " MET_Pt " << m("MET_Pt")
                          << " MET_sumEt " << m("MET_sumEt")
                          << " H_mass " << m("H_mass")
                          << std::endl;
            }
        }

        // 2-lepton
        bool base2LepCSSelection = (    // implementing table 15 of AN2015_168_v12, page 75
            mInt("cutFlow") >= 2         // require Vtype and trigger
            && (mInt("isZmm") || mInt("isZee"))
            && m("V_pt") > 50
        );

        if (base2LepCSSelection) {
            if (////////////////////////// ttbar
                max_hJet_btag > m("tagWPT") //fixme change tight->medium?
                && min_hJet_btag > m("tagWPL")
                && (V_mass > 10 && (V_mass < 75 || V_mass > 120))
            ) {
                *in["controlSample"] = 21;
            } else if (/////////////////// Z + LF
            max_hJet_btag < m("tagWPL") // Change loose to medium?
                && min_hJet_btag < m("tagWPL")
                && HVdPhi > 2.5
                && (V_mass > 75 && V_mass < 105)
            ) {
                if (H_mass > 90 && H_mass <= 150) {
                    *in["controlSample"] = 22;  // normal Z + LF CR
                } else if (m("V_pt") > 150) {
                    *in["controlSample"] = 24;  // keep these event for kinematic fit
                }
            } else if (/////////////////// Z + HF
            max_hJet_btag > m("tagWPT") //Change tight -> medium?
                && min_hJet_btag > m("tagWPL")
                && m("MET_Pt") < 60
                && HVdPhi > 2.5
                && (V_mass > 85 && V_mass <= 97)
            ) {
                if (H_mass <= 90 || H_mass > 150) {
                    *in["controlSample"] = 23;  // normal Z + HF CR
                } else if (
                    *in["controlSample"] != 0 // don't steal SR events
                    && m("V_pt") > 150
                ) {
                    *in["controlSample"] = 25;  // keep these event for kinematic fit
                }
            }
        }
        // end of 2-lepton
    }
    
    //BOOSTED CONTROL SAMPLE SELECTION
    *in["boostedControlSample"] = -1;
    if(int(m("doBoost"))==1 && m("oneMergedJet")){
        if(mInt("FatJetCand_index")>-1
            && m("V_pt")>250
            && mInt("nAddLeptons") == 0
            && m("lepMetDPhi") < 2){
            
            //std::cout << "\t\tDoing boosted control sample" << std::endl;
            //only single lep for now
            if(mInt("isWmunu")||mInt("isWenu")){    
                //Boosted Signal
                if((m("FatJetCand_Msoftdrop_corrected") > 90)
                    && (m("FatJetCand_Msoftdrop_corrected") < 150)
                    && (m("FatJetCand_doubleB") > 0.9 )
                    && (mInt("nBJetsOutsideFatJet") == 0)){
                        *in["boostedControlSample"] = 0;
                }
                //ttbar CR=1
                if((m("FatJetCand_Msoftdrop_corrected") > 50)
                    && (m("FatJetCand_doubleB") > 0.9)
                    && (mInt("nBJetsOutsideFatJet") > 0)){
                        *in["boostedControlSample"] = 11;
                }
                //V+light CR=2
                if((m("FatJetCand_Msoftdrop_corrected") > 50)
                    && (m("FatJetCand_doubleB") < 0.9)
                    && (mInt("nBJetsOutsideFatJet") == 0)){
                        *in["boostedControlSample"] = 12;
                }
                //V+hf CR=3
                if( ((m("FatJetCand_Msoftdrop_corrected") > 50 && m("FatJetCand_Msoftdrop_corrected") < 90) || (m("FatJetCand_Msoftdrop_corrected") > 150 && m("FatJetCand_Msoftdrop_corrected") < 200))
                    && m("FatJetCand_doubleB") > 0.9
                    && mInt("nBJetsOutsideFatJet") == 0){
                        *in["boostedControlSample"] = 13;
                }
            }
        }
    }
}

void VHbbAnalysis::FatJetSelection(){
    if(atLeastOnePreselFatJet){
        float fatJetPtCut, fatJetBBTaggerCut, tau2OverTau1Cut;
        if (mInt("isZnn")) {
            fatJetPtCut       = m("fatJetPtCut_0lepchan");
            fatJetBBTaggerCut = m("fatJetBBTaggerCut_0lepchan");
            tau2OverTau1Cut   = m("tau2OverTau1Cut_0lepchan");
        } else if (mInt("isWmunu") || mInt("isWenu")) {
            fatJetPtCut       = m("fatJetPtCut_1lepchan");
            fatJetBBTaggerCut = m("fatJetBBTaggerCut_1lepchan");
            tau2OverTau1Cut   = m("tau2OverTau1Cut_1lepchan");
        } else if (mInt("isZmm") || mInt("isZee")) {
            fatJetPtCut       = m("fatJetPtCut_2lepchan");
            fatJetBBTaggerCut = m("fatJetBBTaggerCut_2lepchan");
            tau2OverTau1Cut   = m("tau2OverTau1Cut_2lepchan");
        }

        //float highestPt=0;  // I'm assuming highest fat jet pt is best... maybe not?
        float tau2OverTau1=100;
        float highestDoubleB=-2;
        for(int iFatJet=0; iFatJet<mInt("nFatJet"); iFatJet++){
            tau2OverTau1=m("FatJet_tau2",iFatJet)/m("FatJet_tau1",iFatJet);
            if(m("FatJet_pt",iFatJet)>fatJetPtCut
                    && m("FatJet_btagHbb",iFatJet)>fatJetBBTaggerCut
                    && m("FatJet_mass", iFatJet)>50
                    && fabs(m("FatJet_eta",iFatJet))<2.5
                    && tau2OverTau1<tau2OverTau1Cut){//turned off for boosted
                //if(m("FatJet_pt",iFatJet)>highestPt)
                if(m("FatJet_btagHbb",iFatJet)>highestDoubleB){
                    *in["FatJetCand_index"]=iFatJet;
                    highestDoubleB=m("FatJet_btagHbb",iFatJet);
                    /*if(m("FatJet_pt",iFatJet)<highestPt){
                        std::cout<<"this pt is smaller than the last candidate"<<std::endl;
                    }*/
                    //highestPt=m("FatJet_pt",iFatJet);
                }
            }
        }
        if(*in["FatJetCand_index"]>-1){
            fatJetCand.SetPtEtaPhiM(m("FatJet_pt",mInt("FatJetCand_index")),m("FatJet_eta",mInt("FatJetCand_index")),m("FatJet_phi",mInt("FatJetCand_index")),m("FatJet_mass",mInt("FatJetCand_index")));
        }
    }
}

void VHbbAnalysis::BoostedSelection(){
    //bb-tag, tau2/1 and FatJetPT applied in FatJetSelection... is that ok?
    if(mInt("FatJetCand_index")>-1){
        if(m("isZnn")==1){
            *f["nLooseBtagsDR1p5"]=0;
            *f["nLooseBtagsDR1p0"]=0;
            *f["nLooseBtagsDR0p8"]=0;

            for(int iJet=0;iJet<mInt("nJet");iJet++){
                if(m(taggerName,iJet)>m("tagWPL")){
                    TLorentzVector looseBJet;
                    looseBJet.SetPtEtaPhiM(m("Jet_bReg",iJet),m("Jet_eta",iJet),m("Jet_phi",iJet),m("Jet_mass",iJet) * (m("Jet_bReg",iJet) / m("Jet_Pt",iJet) ) );
                    float dR = fatJetCand.DeltaR(looseBJet);
                    if(dR>1.5) *f["nLooseBtagsDR1p5"]=m("nLooseBtagsDR1p5")+1;
                    if(dR>1.0) *f["nLooseBtagsDR1p0"]=m("nLooseBtagsDR1p0")+1;
                    if(dR>0.8) *f["nLooseBtagsDR0p8"]=m("nLooseBtagsDR0p8")+1;
                }
            }

            if(m("nLooseBtagsDR1p5")>0){
                *in["FatJetCand_index"]=-1;
                *b["oneMergedJet"]=false;
            }
        } else if(m("isWmunu")==1 || m("isWenu")==1){
            // not including h mass window selection.  can do with ntuples
            if(m("FatJetVdPhi")>2.9 && m("V_pt")>230){
                unsigned int nAddJets=0;
                *f["Jet_btagCMVA_max1_dR06"]=-99;
                *f["Jet_btagCMVA_max2_dR06"]=-99;
                for(int iJet=0;iJet<mInt("nJet");iJet++){
                    if(m("Jet_lepFilter",iJet)==0) continue;
                    TLorentzVector aJet;
                    aJet.SetPtEtaPhiM(m("Jet_pt",iJet),m("Jet_eta",iJet),m("Jet_phi",iJet),m("Jet_mass",iJet));
                    if(fatJetCand.DeltaR(aJet)>0.6){
                        if(m("Jet_Pt",iJet)>25 && abs(m("Jet_eta",iJet))<2.5){
                            nAddJets++;
                        }
                    } else {
                        float thisBTag=m(taggerName,iJet);
                        if(thisBTag>m("Jet_btagCMVA_max1_dR06")){
                            *f["Jet_btagCMVA_max2_dR06"]=m("Jet_btagCMVA_max1_dR06");
                            *f["Jet_btagCMVA_max1_dR06"]=thisBTag;
                        }else if(thisBTag>m("Jet_btagCMVA_max2_dR06")){
                            *f["Jet_btagCMVA_max2_dR06"]=thisBTag;
                        }
                    }
                }

                if(nAddJets>1){
                    *in["FatJetCand_index"]=-1;
                    *b["oneMergedJet"]=false;
                    return;
                } else {
                    *f["nJets25_dR06"]=nAddJets;
                }
            } else {
                *in["FatJetCand_index"]=-1;
                *b["oneMergedJet"]=false;
            }
        } else if(m("isZmm")==1 || m("isZee")==1){
            // not including h mass window selection.  can do with ntuples
            if( !(m("V_mass")>75 && m("V_mass")<105 && m("V_pt")>100)){  //since all precomputed, invert selection to cut
                *in["FatJetCand_index"]=-1;
                *b["oneMergedJet"]=false;
            }

        }
    }
}


void VHbbAnalysis::ComputeBoostedVariables(){
    if(mInt("FatJetCand_index")>-1){
        *f["nFatJets200"]=0;
        for(int iFatJet=0; iFatJet<mInt("nFatJet"); iFatJet++){
            if(m("FatJet_pt",iFatJet)>200 && fabs(m("FatJet_eta",iFatJet))<2.4) *f["nFatJets200"]=m("nFatJets200")+1;
        }
        *f["nJets30_0lep"]=0;
        *f["nJets30_2lep"]=0;
        *f["nJets25_dR06"]=0;
        for(int iJet=0; iJet<mInt("nJet"); iJet++){
            if(m("Jet_lepFilter",iJet)==0) continue;
            if(m("Jet_Pt",iJet)>30 && fabs(m("Jet_eta",iJet))<2.4){
                if( (m("Jet_puId",iJet)>6) || m("Jet_Pt",iJet)>50) *f["nJets30_0lep"]=m("nJets30_0lep")+1;
                if(iJet!=mInt("hJetInd2")&&iJet!=mInt("hJetInd1")&&(m("Jet_puId",iJet)>6 || m("Jet_Pt",iJet)>50)) *f["nJets30_2lep"]=m("nJets30_2lep")+1;
            }
        }
        //Compute number of b jets outside the fatjet
        TLorentzVector FatJet;
        FatJet.SetPtEtaPhiM(m("FatJet_pt",mInt("FatJetCand_index")),
            m("FatJet_eta",mInt("FatJetCand_index")),
            m("FatJet_phi",mInt("FatJetCand_index")),
            m("FatJet_msoftdrop",mInt("FatJetCand_index")) );

        int nBJetsOutsideFatJet = 0;
        for(int iJet=0;iJet<mInt("nJet");iJet++){
            if(m("Jet_lepFilter",iJet)
                && m("Jet_Pt",iJet)>25 
                && abs(m("Jet_eta",iJet))<2.5
                && m(taggerName,iJet)>m("tagWPL")){

                TLorentzVector aJet;
                aJet.SetPtEtaPhiM(m("Jet_pt",iJet),m("Jet_eta",iJet),m("Jet_phi",iJet),m("Jet_mass",iJet));
                if(FatJet.DeltaR(aJet)>0.8){//Is 0.6 right? fatjet 'radius' is 0.4, jet 'radius' is 0.2
                    nBJetsOutsideFatJet++;
                }
            }
        }
        *in["nBJetsOutsideFatJet"] = nBJetsOutsideFatJet;




        *f["FatJetCand_pt"]=m("FatJet_pt",mInt("FatJetCand_index"));
        *f["FatJetCand_eta"]=m("FatJet_eta",mInt("FatJetCand_index"));
        *f["FatJetCand_phi"]=m("FatJet_phi",mInt("FatJetCand_index"));
        *f["FatJetCand_tau21"]=m("FatJet_tau2",mInt("FatJetCand_index"))/m("FatJet_tau1",mInt("FatJetCand_index"));
        *f["FatJetCand_tau32"]=m("FatJet_tau3",mInt("FatJetCand_index"))/m("FatJet_tau2",mInt("FatJetCand_index"));
        *f["FatJetCand_tau1"]=m("FatJet_tau1",mInt("FatJetCand_index"));
        *f["FatJetCand_tau2"]=m("FatJet_tau2",mInt("FatJetCand_index"));
        *f["FatJetCand_tau3"]=m("FatJet_tau3",mInt("FatJetCand_index"));
        *f["FatJetCand_doubleB"]=m("FatJet_btagHbb",mInt("FatJetCand_index"));
        *f["FatJetCand_Msoftdrop_corrected"]=m("FatJet_msoftdrop",mInt("FatJetCand_index"));
        *f["FatJetCand_deepTagMD_bbvsLight"]=m("FatJet_deepTagMD_bbvsLight",mInt("FatJetCand_index"));
        float TvsQCD=m("FatJet_deepTagMD_TvsQCD",mInt("FatJetCand_index"));
        float HbbvsQCD=m("FatJet_deepTagMD_HbbvsQCD",mInt("FatJetCand_index"));
        *f["FatJetCand_deepTagMD_bbvsTop"]==1/(1+(TvsQCD/HbbvsQCD)*(1-HbbvsQCD)/(1-TvsQCD));
        //Compute boosted-specific BDT variables
         
        //define lepMetDPhi for other channels (temporary until channel-specific boosted BDTs are trained
        if(mInt("isWenu")==1 || mInt("isWmunu")==1){
            *f["lepMetDPhiBoostedBDT"] = m("lepMetDPhi");
        }
        else{
            *f["lepMetDPhiBoostedBDT"] = 0.2;
        }
    }
}

void VHbbAnalysis::UpdateJetPts() {
    for (int i = 0; i < mInt("nJet"); i++) {
        f["Jet_bReg"][i] = m("Jet_PtReg",i);
    }
}

std::pair<int,int> VHbbAnalysis::HighestPtGoodElectronsOppCharge(float min_pt, float max_rel_iso, float idcut, bool isOneLepton) {
    int first = -1;
    for (int i = 0; i<mInt("nElectron"); i++) {
        bool passEleIDCut = 0;
        if(isOneLepton){
            if (m("dataYear")==2018){
                passEleIDCut = m("Electron_mvaFall17V2Iso_WP80",i) >0 ;
            } else if(m("dataYear")==2017){
                passEleIDCut = m("Electron_mvaFall17V2Iso_WP80",i) > 0 ;
            } else if(m("dataYear")==2016){
                passEleIDCut = m("Electron_mvaFall17V2Iso_WP80",i) > 0;
            } else {
                std::cout<<"In HighestPtGoodElectronsOppCharge... what year is it? "<<m("dataYear")<<std::endl;
            }
        } else {
            if (m("dataYear")==2018){
                passEleIDCut = m("Electron_mvaFall17V2Iso_WP90",i) >0 ;
            } else if(m("dataYear")==2017){
                passEleIDCut = m("Electron_mvaFall17V2Iso_WP90",i) > 0;
            } else if(m("dataYear")==2016){
                passEleIDCut = m("Electron_mvaFall17V2Iso_WP90",i) >0;
            } else {
                std::cout<<"In HighestPtGoodElectronsOppCharge... what year is it? "<<m("dataYear")<<std::endl;
            }
       }
        if (fabs(m("Electron_eta",i)) < m("eletacut")
            && m("Electron_pt",i) > min_pt
            && m("Electron_pfRelIso03_all",i)< max_rel_iso
            //&& in["selLeptons_eleMVAIdSppring16GenPurp"][i] >= idcut
            && passEleIDCut)
        {
            // leptons are pt sorted
            if (first == -1) {
                first = i;
            } else if (mInt("Electron_charge",i)*mInt("Electron_charge",first) < 0) {
                return std::make_pair(first, i);
            }

        }
    }
    return std::make_pair(first, -1);
}

std::pair<int,int> VHbbAnalysis::HighestPtGoodMuonsOppCharge(float min_pt, float max_rel_iso, bool isOneLepton) {
    int first = -1;
    for (int i = 0; i<mInt("nMuon"); i++) {
        if (fabs(m("Muon_eta",i)) < m("muetacut")
            && m("Muon_pt",i) > min_pt
            && m("Muon_pfRelIso04_all",i) < max_rel_iso
            //&& in["selLeptons_looseIdPOG"][i] >= *f["muidcut"]
            && ( !isOneLepton || m("Muon_tightId",i) > 0) )
        {
            // leptons are pt sorted
            if (first == -1) {
                first = i;
            } else if (mInt("Muon_charge",i)*mInt("Muon_charge",first) < 0) {
                return std::make_pair(first, i);
            }
        }
    }
    return std::make_pair(first, -1);
}


// FIXME can this function be removed ??
bool VHbbAnalysis::ElectronSelection(int nLep){
    bool selectEvent=true;
    if(debug>1000){
        std::cout<<"nLep "<<nLep<<std::endl;
        std::cout<<"Running Wenu selections"<<std::endl;
        std::cout<<"*in[\"nselLeptons\"] "<<*in["nselLeptons"]<<std::endl;
        std::cout<<"d[\"selLeptons_pt\"][0] "<<f["selLeptons_pt"][0]<<std::endl;
        std::cout<<"in[\"selLeptons_pdgId\"] "<<in["selLeptons_pdgId"][0]<<std::endl;
        std::cout<<"d[\"selLeptons_relIso03\"] "<<f["selLeptons_relIso03"][0]<<std::endl;
        std::cout<<"*f[\"MET_Pt\"] "<<*f["MET_Pt"]<<std::endl;
        std::cout<<"*f[selLeptons_eleSieie_0] = "<<*f["selLeptons_eleSieie_0"]<<std::endl;
        std::cout<<"*f[hJets_pt_1] = "<<*f["hJets_pt_1"]<<std::endl;
    }
    // there is only one selected electron for Vtype == 3 which is the electron tag
    *in["elInd1"] = -1;
    *in["elInd2"] = -1;
    if(nLep==1){
        float elMaxPt = 0; // max pt of the electrons we select
        for (int i =0; i<mInt("nselLeptons"); i++) {
            if(fabs(mInt("selLeptons_pdgId",i))==11
                && m("selLeptons_pt",i)      > m("eptcut")
                && fabs(m("selLeptons_eta",i)) < m("eletacut")
                && m("selLeptons_relIso03",i)< m("erelisocut")
                //&& in["selLeptons_eleCutIdSpring15_25ns_v1"][i] >= *f["elidcut"]
                //&& in["selLeptons_tightId"][i] >= *f["elidcut"]
                //&& in["selLeptons_eleMVAIdSpring15Trig"][i] >= *f["elidcut"]
                && m("selLeptons_eleMVAIdSppring16GenPurp",i) >= m("elidcut")
                ){
                if (m("selLeptons_pt",i) > elMaxPt) {
                    elMaxPt = m("selLeptons_pt",i);
                    *in["elInd1"] = i;
                }
            }
        }
        if (mInt("elInd1") == -1) selectEvent = false;
    } else if(nLep==2){
        float elMaxPt = 0; // max pt of the electrons we select
        float elMax2Pt = 0; // 2nd max pt of the electrons we select
        for (int i =0; i<mInt("nselLeptons"); i++) {
            if(fabs(mInt("selLeptons_pdgId",i))==11
                && fabs(m("selLeptons_eta",i)) < m("eletacut")
                && m("selLeptons_relIso03",i)< m("erelisocut")
                && mInt("selLeptons_eleMVAIdSppring16GenPurp",i) >= m("elidcut")
                ){
                if (m("selLeptons_pt",i) > elMaxPt &&  m("selLeptons_pt",i) > m("e1ptcut")) {
                    elMaxPt = m("selLeptons_pt",i);
                    *in["elInd1"] = i;
                }
                else if (mInt("elInd1") > -1 && m("selLeptons_pt",i) > elMax2Pt && (mInt("selLeptons_charge",i)*mInt("selLeptons_charge",mInt("elInd1"))) < 0) {
                    elMax2Pt = m("selLeptons_pt",i);
                    *in["elInd2"] = i;
                }
            }
        }
        if (debug>1000) {
            std::cout<<"elInd1/2 = "<<mInt("elInd1")<<", "<<mInt("elInd2")<<std::endl;
        }
        if (mInt("elInd1") == -1 || mInt("elInd2") == -1) selectEvent = false;
    } else {
        std::cout<<"nLep is not 1 or 2 "<<nLep<<" - fail it."<<std::endl;
        selectEvent = false;
    }

    return selectEvent;
}


// FIXME can this function be removed ??
bool VHbbAnalysis::MuonSelection(int nLep){

    bool selectEvent=true;
    if(debug>1000){
        std::cout<<"nLep "<<nLep<<std::endl;
        std::cout<<"Running Wmunu selections"<<std::endl;
        std::cout<<"*in[\"nselLeptons\"] "<<*in["nselLeptons"]<<std::endl;
        std::cout<<"d[\"selLeptons_pt\"][0] "<<f["selLeptons_pt"][0]<<std::endl;
        std::cout<<"in[\"selLeptons_pdgId\"] "<<in["selLeptons_pdgId"][0]<<std::endl;
        std::cout<<"d[\"selLeptons_relIso04\"] "<<f["selLeptons_relIso04"][0]<<std::endl;
        std::cout<<"*d[\"MET_Pt\"] "<<*f["MET_Pt"]<<std::endl;
    }

    // there is only one selected electron for Vtype == 3 which is the electron tag
    // FIXME add configurable cuts

    if(nLep==1){
        *in["muInd1"] = -1;
        float muMaxPt = 0; // max pt of the muons we select
        for (int i=0; i<mInt("nselLeptons"); i++) {
            if(fabs(mInt("selLeptons_pdgId",i))==13
                && m("selLeptons_pt",i)      > m("muptcut")
                && fabs(m("selLeptons_eta",i)) < m("muetacut")
                && m("selLeptons_relIso04",i)< m("murelisocut")
                && mInt("selLeptons_tightId",i) >= m("muidcut")
                ){
                if (m("selLeptons_pt",i) > muMaxPt) {
                    muMaxPt = m("selLeptons_pt",i);
                    *in["muInd1"] = i;
                }
            }
        }
        if (mInt("muInd1") == -1) selectEvent = false;
    } else if(nLep==2) {
        *in["muInd1"] = -1;
        *in["muInd2"] = -1;
        float muMaxPt = 0; // max pt of the muons we select
        float muMax2Pt = 0; // 2nd max pt of the muons we select
        for (int i=0; i<mInt("nselLeptons"); i++) {
            if(fabs(mInt("selLeptons_pdgId",i))==13
                && fabs(m("selLeptons_eta",i)) < m("muetacut")
                && m("selLeptons_relIso04",i)< m("murelisocut")
                && mInt("selLeptons_looseIdPOG",i) >= m("muidcut")
                ){
                if (m("selLeptons_pt",i) > muMaxPt &&  m("selLeptons_pt",i) > m("mu1ptcut")) {
                    muMaxPt = m("selLeptons_pt",i);
                    *in["muInd1"] = i;
                } else if (mInt("muInd1")>-1 && m("selLeptons_pt",i) > muMax2Pt
                  && (mInt("selLeptons_charge",i)*mInt("selLeptons_charge",mInt("muInd1"))) < 0) {
                    muMax2Pt = m("selLeptons_pt",i);
                    *in["muInd2"] = i;
                }
            }
        }
        if (mInt("muInd1") == -1 || mInt("muInd2") == -1) selectEvent = false;
    } else {
        std::cout<<"nLep is not 1 or 2 "<<nLep<<" - fail it."<<std::endl;
        selectEvent = false;
    }

    return selectEvent;
}

int VHbbAnalysis::UpdatedVType() {

    m("Vtype_new");
    *f["Vtype_new"] = mInt("Vtype");

    // muon channels are good
    if (mInt("Vtype") == 0 || mInt("Vtype") == 2) {
        return *in["Vtype"];
    }

    // collect good electrons
    std::vector<int> good_electrons;
    for (int i = 0; i < mInt("nElectron"); i++) {
        if (m("Electron_pt",i) > 15
            && m("Electron_mvaSpring16GP_WP90",i) > 0)
        {
            good_electrons.push_back(i);
        }
    }

    // check for Vtype 1
    if (good_electrons.size() >= 2) {
        if (m("Electron_pt",good_electrons[0]) > 20) {
            int ele0_charge = mInt("Electron_charge",good_electrons[0]);
            for (unsigned i=1; i<good_electrons.size(); i++) {
                if (ele0_charge * mInt("Electron_charge",good_electrons[i]) < 0) {
                    *f["Vtype_new"] = 1;
                    return 1;
                }
            }
        } else {
            *f["Vtype_new"] = -1;
            return -1;  // (leading electron pt between 15 and 20)
        }
    }

    // can still be Wenu channel
    if (mInt("Vtype") == 3) {
        return 3;
    }

    // is Vtype 4 but shouldn't be?
    if (mInt("Vtype") == 4 && good_electrons.size()) {
        *f["Vtype_new"] = 5;
        return 5;
    }

    // isn't Vtype 4 but should be?
    // nope, because the new electrons are looser. Hence no migration from Vtype 5 to 4

    return *in["Vtype"];
}


bool VHbbAnalysis::PassVTypeAndTrigger(int vtype) {

    if (vtype < 0 || vtype > 4) {
        return false;
    }

    if(debug>10000) std::cout<<"dataYear doICHEP sampleIndex "<<*f["dataYear"]<<" "<<*f["doICHEP"]<<" "<<*in["sampleIndex"]<<std::endl;
    //FIXME configure for different channels
    //TODO it would be nicer to reverse the logic ie. if (trg1 || trg2 || trg3) {return true;} else if ....
    if (m("dataYear") == 2015) {
        // for 2015 V21 ntuples
        if (mInt("HLT_BIT_HLT_Ele23_WPLoose_Gsf_v") != 1
            && mInt("HLT_BIT_HLT_IsoMu20_v") != 1
            && mInt("HLT_BIT_HLT_IsoTkMu20_v") != 1
           ) {
            return false;
        }
    } else if (m("dataYear") == 2016 && m("doICHEP") == 1) {
        // for 2016 V22 ntuples there is no MC HLT simulation we can use, have to just apply the data trigger efficiency
        if (mInt("sampleIndex") == 0) {
            // regular triggers for 2016 data, but none applied to MC
            if (mInt("HLT_BIT_HLT_Ele27_eta2p1_WPLoose_Gsf_v") != 1
                && mInt("HLT_BIT_HLT_IsoMu22_v") != 1
                && mInt("HLT_BIT_HLT_IsoTkMu22_v") != 1
               ) {
                return false;
            }
        }
    } else if (m("dataYear") == 2016) {
        // 0-lepton - no MC efficiencies available, so apply trigger in data only
        if (vtype == 4 && mInt("sampleIndex") == 0
            && mInt("HLT_PFMET110_PFMHT110_IDTight") != 1
            && mInt("HLT_PFMET120_PFMHT120_IDTight") != 1
            && mInt("HLT_PFMET170_NoiseCleaned") != 1
            && mInt("HLT_PFMET170_BeamHaloCleaned") != 1
            && mInt("HLT_PFMET170_HBHECleaned") != 1
           ) {
            return false;
        }

        // 1-lepton
        if ((vtype == 2 || vtype == 3)
            && mInt("HLT_Ele27_eta2p1_WPTight_Gsf")!= 1
            && mInt("HLT_IsoMu24")                 != 1
            && mInt("HLT_IsoTkMu24")               != 1
           ) {
            return false;
        }
        // 2-lepton
        if ((vtype == 0 || vtype == 1)
            && mInt("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL")          != 1
            && mInt("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ")       != 1
            && mInt("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL")        != 1
            && mInt("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ")     != 1
            && mInt("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ") != 1
           ) {
            return false;
        }
    } else if (m("dataYear") == 2017) {
        // 0-lepton
        if (vtype == 4
            && mInt("HLT_PFMET120_PFMHT120_IDTight") != 1
            && mInt("HLT_PFMET120_PFMHT120_IDTight_PFHT60") != 1
           ) {
            return false;
        }
        // 1-lepton
        if ((vtype == 2 || vtype == 3)
            && mInt("HLT_Ele32_WPTight_Gsf_L1DoubleEG") != 1
            && mInt("HLT_IsoMu27") != 1
           ) {
            return false;
        }
        // 2-lepton
        if ((vtype == 0 || vtype == 1)
            && mInt("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8") != 1
            && mInt("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8") != 1
            && mInt("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL") != 1
           ) {
            return false;
        }
       // if ((vtype == 0 || vtype == 1)
       //     && *b["HLT_IsoMu27"] != 1
       //     && *b["HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL"] != 1
       //    ) {
       //     return false;
        //}
    } else if (m("dataYear") == 2018) {
        // 0-lepton
        if (vtype == 4
            && mInt("HLT_PFMET120_PFMHT120_IDTight") != 1
           ) {
            return false;
        }
        // 1-lepton
        if ((vtype == 2 || vtype == 3)
            && mInt("HLT_Ele32_WPTight_Gsf") != 1
            && mInt("HLT_IsoMu24") != 1
           ) {
            return false;
        }
        // 2-lepton
        if ((vtype == 0 || vtype == 1)
            && mInt("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8") != 1
            && mInt("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL") != 1
           ) {
            return false;
        }
    } else {
        std::cout<<"What is this?  Run 1?  Run 3?  2018? "<<std::endl;
    }
    return true;
}

//Member of VHbbAnalysis that convert the tagger type specified in settings.txt into the corresponding string
void VHbbAnalysis::SetTaggerName( float taggerType ){
    std::string name;

    if( m("taggerType")==0 ){
        name.append("Jet_btagDeepB");
    } else if(  m("taggerType")==1 ){
        name.append("Jet_btagCMVA");
    } else if( m("taggerType")==2 ){
        name.append("Jet_btagCSVV2");
    } else {
        std::cout << "Invalid tagger type, setting the tagger name to its defauld value: CMVA" << std::endl;
        name.append("Jet_btagCMVA");
    }

    taggerName = name;
}

//function used to get the EWK correction factor
float VHbbAnalysis::GetVHEWKCorrFactor( float V_pt, TH1D* hist ){
    if(debug>10000){
        std::cout<<"GetVHEWKCorrFactor V_pt"<<V_pt<<std::endl;
        std::cout<<"name of hist "<<hist->GetName()<<std::endl;
    }
    int ibin = hist->GetXaxis()->FindBin(V_pt);
    if(debug>10000){
        std::cout<<"ibin "<<ibin<<std::endl;
        std::cout<<"nbins "<< hist->GetNbinsX()<<std::endl;
    }
    if (ibin < 1){
        return hist->GetBinContent(1);
    }
    if (ibin > hist->GetNbinsX()){
        return hist->GetBinContent(hist->GetNbinsX());
    }

    return hist->GetBinContent(ibin);
}

std::pair<int,int> VHbbAnalysis::HighestPtBJets(){
    std::pair<int,int> pair(-1,-1);

    for(int i=0; i<mInt("nJet"); i++){
        if( (mInt("Jet_puId",i) > 6 || m("Jet_Pt",i)>50)
            && m("Jet_bReg",i)>m("j1ptCut")
            && m("Jet_btagCSVV2",i)>m("j1Btag")&&fabs(m("Jet_eta",i))<=m("JetEtaCut")) {
            if( pair.first == -1 ) {
                pair.first = i;
            } else if(m("Jet_bReg",pair.first)<m("Jet_bReg",i)){
                pair.first = i;
            }
        }
    }

    for(int i=0; i<mInt("nJet"); i++){
        if(i==pair.first) continue;
        if( (mInt("Jet_puId",i) > 6 || m("Jet_Pt",i)>50)
            && m("Jet_bReg",i)>m("j2ptCut")
            && m("Jet_btagCSVV2",i)>m("j2Btag")&&fabs(m("Jet_eta",i))<m("JetEtaCut")) {
            if( pair.second == -1 ) {
                pair.second = i;
            } else if(m("Jet_bReg",pair.second)<m("Jet_bReg",i)){
                pair.second = i;
            }
        }
    }

    return pair;
}


std::pair<int,int> VHbbAnalysis::HighestCSVBJets(float j1ptCut, float j2ptCut){
    std::pair<int,int> pair(-1,-1);

    for(int i=0; i<mInt("nJet"); i++){
        if( (mInt("Jet_puId",i) > 6 || m("Jet_Pt",i)>50)
            && m("Jet_bReg",i)>j1ptCut
            &&fabs(m("Jet_eta",i))<=m("JetEtaCut")) {
            if( pair.first == -1 ) {
                pair.first = i;
            } else if(m("Jet_btagCSVV2",pair.first)<m("Jet_btagCSVV2",i)){
                pair.first = i;
            }
        }
    }

    for(int i=0; i<mInt("nJet"); i++){
        if(i==pair.first) continue;
        if( (mInt("Jet_puId",i) > 6 || m("Jet_Pt",i)>50)
            && m("Jet_bReg",i)>j2ptCut
            &&fabs(m("Jet_eta",i))<m("JetEtaCut")) {
            if( pair.second == -1 ) {
                pair.second = i;
            } else if(m("Jet_btagCSVV2",pair.second)<m("Jet_btagCSVV2",i)){
                pair.second = i;
            }
        }
    }

    // different pt threshold can set the highest CSV value into pair.second
    if (pair.first > -1 && pair.second > -1 && m("Jet_btagCSVV2",pair.first) < m("Jet_btagCSVV2",pair.second)) {
        pair = std::make_pair(pair.second, pair.first);
    }

    return pair;
}


std::pair<int,int> VHbbAnalysis::HighestCMVABJets(float j1ptCut, float j2ptCut){
    std::pair<int,int> pair(-1,-1);

    for(int i=0; i<mInt("nJet"); i++){
        if( (mInt("Jet_puId",i) > 6 || m("Jet_Pt",i)>50)
            && m("Jet_bReg",i)>j1ptCut
            &&fabs(m("Jet_eta",i))<=m("JetEtaCut")) {
            if( pair.first == -1 ) {
                pair.first = i;
            } else if(m("Jet_btagCMVA",pair.first)<m("Jet_btagCMVA",i)){
                pair.first = i;
            }
        }
    }

    for(int i=0; i<mInt("nJet"); i++){
        if(i==pair.first) continue;
        if( (mInt("Jet_puId",i) > 6 || m("Jet_Pt",i)>50)
            && m("Jet_bReg",i)>j2ptCut
            &&fabs(m("Jet_eta",i))<m("JetEtaCut")) {
            if( pair.second == -1 ) {
                pair.second = i;
            } else if(m("Jet_btagCMVA",pair.second)<m("Jet_btagCMVA",i)){
                pair.second = i;
            }
        }
    }

    // different pt threshold can set the highest CMVA value into pair.second
    if (pair.first > -1 && pair.second > -1 && m("Jet_btagCMVA",pair.first) < m("Jet_btagCMVA",pair.second)) {
        pair = std::make_pair(pair.second, pair.first);
    }

    return pair;
}

std::pair<int,int> VHbbAnalysis::HighestDeepCSVBJets(float j1ptCut, float j2ptCut){
    std::pair<int,int> pair(-1,-1);

    for(int i=0; i<mInt("nJet"); i++){
        if( (mInt("Jet_puId",i) > 6 || m("Jet_Pt",i)>50)
            && m("Jet_bReg",i)>j1ptCut
            &&fabs(m("Jet_eta",i))<=m("JetEtaCut")) {
            if( pair.first == -1 ) {
                pair.first = i;
            } else if(m("Jet_btagDeepB",pair.first)<m("Jet_btagDeepB",i)){
                pair.first = i;
            }
        }
    }

    for(int i=0; i<mInt("nJet"); i++){
        if(i==pair.first) continue;
        if( (mInt("Jet_puId",i) > 6 || m("Jet_Pt",i)>50)
            && m("Jet_bReg",i)>j2ptCut
            &&fabs(m("Jet_eta",i))<m("JetEtaCut")) {
            if( pair.second == -1 ) {
                pair.second = i;
            } else if(m("Jet_btagDeepB",pair.second)<m("Jet_btagDeepB",i)){
                pair.second = i;
            }
        }
    }

    // different pt threshold can set the highest CMVA value into pair.second
    if (pair.first > -1 && pair.second > -1 && m("Jet_btagDeepB",pair.first) < m("Jet_btagDeepB",pair.second)) {
        pair = std::make_pair(pair.second, pair.first);
    }

    return pair;
}


//New function that should cover all tagger possibilities
std::pair<int,int> VHbbAnalysis::HighestTaggerValueBJets(float j1ptCut, float j2ptCut, std::string taggerName){
    std::pair<int,int> pair(-1,-1);

    for(int i=0; i<mInt("nJet"); i++){
        if(m("Jet_lepFilter",i) && (mInt("Jet_puId",i) > 6 || m("Jet_Pt",i)>50)
            && m("Jet_bReg",i)>j1ptCut
            &&fabs(m("Jet_eta",i))<=m("JetEtaCut")) {
            if( pair.first == -1 ) {
                pair.first = i;
            } else if(m(taggerName,pair.first)<m(taggerName,i)){
                pair.first = i;
            }
        }
    }

    for(int i=0; i<mInt("nJet"); i++){
        if(i==pair.first) continue;
        if(m("Jet_lepFilter",i) && (mInt("Jet_puId",i) > 6 ||m("Jet_Pt",i)>50)
            && m("Jet_bReg",i)>j2ptCut
            &&fabs(m("Jet_eta",i))<m("JetEtaCut")) {
            if( pair.second == -1 ) {
                pair.second = i;
            } else if(m(taggerName,pair.second)<m(taggerName,i)){
                pair.second = i;
            }
        }
    }

    // different pt threshold can set the highest CMVA value into pair.second
    if (pair.first > -1 && pair.second > -1 && m(taggerName,pair.first) < m(taggerName,pair.second)) {
        pair = std::make_pair(pair.second, pair.first);
    }

    return pair;
}

std::pair<int,int> VHbbAnalysis::HighestPtJJBJets(){
    std::pair<int,int> pair(-1,-1);

    // dont implement the btag CSV cuts until we've already selected the highest pT(jj) jets
    double maxPtJJ = 0.;
    for (int i=0; i<mInt("nJet"); i++) {
        if( (mInt("Jet_puId",i) > 6 || m("Jet_Pt",i)>50)
            && m("Jet_bReg",i)>m("j1ptCut")
            && fabs(m("Jet_eta",i))<m("JetEtaCut")) {
            TLorentzVector jet1;
            jet1.SetPtEtaPhiM(m("Jet_bReg",i),m("Jet_eta",i),m("Jet_phi",i),m("Jet_mass",i) * (m("Jet_bReg",i) / m("Jet_pt",i) ) );
            for (int j=0; j<mInt("nJet"); j++) {
                if (i == j) continue;
                if( (mInt("Jet_puId",j) > 6 || m("Jet_Pt",j)>50)
                    && m("Jet_bReg",j)>m("j2ptCut")
                    && fabs(m("Jet_eta",j))<m("JetEtaCut")) {
                    TLorentzVector jet2;
                    jet2.SetPtEtaPhiM(m("Jet_bReg",j),m("Jet_eta",j),m("Jet_phi",j),m("Jet_mass",j) * (m("Jet_bReg",j) / m("Jet_pt",j) ) );
                    TLorentzVector jj = jet1 + jet2;
                    double ptJJ = jj.Pt();
                    if (ptJJ >= maxPtJJ) {
                        if (j == pair.first && pair.second == i) {
                            // we've already picked this pair in the other order, order by CSV
                            // quick sanity check, make sure it's the same maxPtJJ
                            if (ptJJ != maxPtJJ) {
                                std::cout<<"Picked both orderings of highest pT(jj) jets, but the orderings don't have the same pT(jj)!!"<<std::endl;
                                std::cout<<"ptJJ = "<<ptJJ<<std::endl;
                                std::cout<<"maxPtJJ = "<<maxPtJJ<<std::endl;
                            }
                            if (m("Jet_btagCSVV2",j) > m("Jet_btagCSVV2",i)) {
                            //if (f["Jet_bReg"][j] > f["Jet_pt"][i]) {
                                pair.first = j;
                                pair.second = i;
                            }
                            else {
                                pair.first = i;
                                pair.second = j;
                            }
                        }
                        else {
                            pair.first = i;
                            pair.second = j;
                            maxPtJJ = ptJJ;
                        }
                    }
                }
            }
        }
    }
    // important to cut on CSV here to kill TTbar
    if(pair.first != -1){
        if (m("Jet_btagCSVV2",pair.first) < m("j1Btag")) {
            pair.first = -1;
        }
    }
    if(pair.second != -1){
        if (m("Jet_btagCSVV2",pair.second) < m("j2Btag")) {
            pair.second = -1;
        }
    }
    return pair;
}

double VHbbAnalysis::GetRecoTopMass(TLorentzVector Obj, bool isJet, int useMET, bool regPT, bool smearedPT) {
    if (regPT==1 && smearedPT==0){
        std::cout<<"Don't try to use regressed pT if you're not using smearing!"<<std::endl;    
    }
    // Try to reconstruct the top in ttbar with leptonic W decay
    // if isJet is true, construct top as given jet + closest lepton
    // if isJet is false, construct top as given lepton + closest lepton
    //
    // if useMET is 1, construct the top using the jet + lepton + met and take the transverse mass
    // if useMET is 2, construct the top using the jet + lepton + met, calculate met pZ by assuming mW

    double minDR = 999;
    int ObjClosestIndex = -1; // index of closest lepton if isJet, closet jet otherwise

    TLorentzVector Obj2; // closest lepton if isJet, closest jet otherwise
    TLorentzVector Top;
    //if (isJet) {
    //    // look in aleptons too FIXME
    //    // find the closest lepton to the given jet
    //    for (int i=0; i<*in["nselLeptons"]; i++) {
    //        TLorentzVector l;
    //        l.SetPtEtaPhiM(f["selLeptons_pt"][i], f["selLeptons_eta"][i], f["selLeptons_phi"][i], f["selLeptons_mass"][i] );
    //        double d1 = l.DeltaR(Obj);
    //        if (d1 <= minDR) {
    //            minDR = d1;
    //            ObjClosestIndex = i;
    //         }
    //    }
    //    if (ObjClosestIndex!=-1) {
    //        Obj2.SetPtEtaPhiM(f["selLeptons_pt"][ObjClosestIndex], f["selLeptons_eta"][ObjClosestIndex], f["selLeptons_phi"][ObjClosestIndex], f["selLeptons_mass"][ObjClosestIndex]);
    //    }
    //    else return -999;
    //}
    //else
    // find closest jet to the given lepton
    float thisPT=0;
    for (int i=0; i<mInt("nJet"); i++) {
        if(regPT){
            thisPT=m("Jet_bReg",i);
        } else if(smearedPT) {
            thisPT=m("Jet_Pt",i);
        } else {
            thisPT=m("Jet_pt",i);
        }
        if (thisPT< 30 || m(taggerName,i) < 0.5) continue; // only consider jets with some minimal preselection
        TLorentzVector j;
        if(regPT){
            j.SetPtEtaPhiM(thisPT, m("Jet_eta",i), m("Jet_phi",i), m("Jet_mass",i) * (m("Jet_bReg",i) / m("Jet_Pt",i) )  );
        } else {
            j.SetPtEtaPhiM(thisPT, m("Jet_eta",i), m("Jet_phi",i), m("Jet_mass",i)  );
        }
        double d1 = j.DeltaR(Obj);
        if (d1 <= minDR) {
            minDR = d1;
            ObjClosestIndex = i;
         }
    }
    if (ObjClosestIndex!=-1) {
        if(regPT){
            thisPT=m("Jet_bReg",ObjClosestIndex);
            Obj2.SetPtEtaPhiM(thisPT, m("Jet_eta",ObjClosestIndex), m("Jet_phi",ObjClosestIndex), m("Jet_mass",ObjClosestIndex) * (thisPT / m("Jet_Pt",ObjClosestIndex) ) );
        } else if (smearedPT){
            thisPT=m("Jet_Pt",ObjClosestIndex);
            Obj2.SetPtEtaPhiM(thisPT, m("Jet_eta",ObjClosestIndex), m("Jet_phi",ObjClosestIndex), m("Jet_mass",ObjClosestIndex));
        } else {
            thisPT=m("Jet_pt",ObjClosestIndex);
            Obj2.SetPtEtaPhiM(thisPT, m("Jet_eta",ObjClosestIndex), m("Jet_phi",ObjClosestIndex), m("Jet_mass",ObjClosestIndex));
        }
    } else return -999;


    if (useMET==0) {
        Top = Obj + Obj2;
    }

    if (useMET==1) {
        // try top = lep + jet + met
        // two-particle Mt = sqrt(Et**2 - pt**2)
        TLorentzVector MET;
        if (smearedPT){
            MET.SetPtEtaPhiM(m("MET_Pt"),0.,m("MET_Phi"),0.);
        } else {
            MET.SetPtEtaPhiM(m("MET_pt"),0.,m("MET_phi"),0.);
        }
        TLorentzVector Obj_transverse, Obj2_transverse; // can only consider transverse (x-y plane) 4-vector components if using MET
        Obj_transverse.SetPxPyPzE(Obj.Px(),Obj.Py(),0,TMath::Sqrt(TMath::Power(Obj.M(),2) + TMath::Power(Obj.Pt(),2)));
        Obj2_transverse.SetPxPyPzE(Obj2.Px(),Obj2.Py(),0,TMath::Sqrt(TMath::Power(Obj2.M(),2) + TMath::Power(Obj2.Pt(),2)));
        Top = Obj_transverse + Obj2_transverse + MET;
    }else if (useMET==2) {
        TLorentzVector MET;
        if (smearedPT){
            MET.SetPtEtaPhiM(m("MET_Pt"),0.,m("MET_Phi"),0.);
        } else {
            MET.SetPtEtaPhiM(m("MET_pt"),0.,m("MET_phi"),0.);
        }
        TLorentzVector lep, jet;
        if (isJet) {
            lep = Obj2;
            jet = Obj;
        }
        else {
            lep = Obj;
            jet = Obj2;
        }

        TLorentzVector neutrino = getNu4Momentum(lep, MET);
        Top = lep + jet + neutrino;
        TLorentzVector W = lep + neutrino;
    }

    return Top.M();
}

// count the number of WPs passed
int VHbbAnalysis::BtagWPForJet(int jetIndex){
    int nWPsPassed=0;
    if(jetIndex>-1){
        nWPsPassed+=int(m(taggerName,jetIndex)>m("tagWPL")); 
        if(nWPsPassed>0){
            nWPsPassed+=int(m(taggerName,jetIndex)>m("tagWPM")); 
            if(nWPsPassed>1){
                nWPsPassed+=int(m(taggerName,jetIndex)>m("tagWPT"));
            }
        } 
    }
    return nWPsPassed;
}

float VHbbAnalysis::ReWeightMC(int nPU){

double data2[50]={
3.72080752299e-07,
1.47718953836e-05,
5.79841527807e-05,
0.000132465155152,
0.000194029101597,
0.000266928406893,
0.000537390953436,
0.00249222471316,
0.00713605338046,
0.0147893534623,
0.0232230204065,
0.0317807128657,
0.0424387961713,
0.0541178963448,
0.0653538845508,
0.0749520482888,
0.0815302475861,
0.0845985629121,
0.0843143028301,
0.0807132646267,
0.0742418229119,
0.0659244351143,
0.0565195845436,
0.046354924125,
0.0359433811537,
0.026208951089,
0.0180189631761,
0.0117407079361,
0.00726631530962,
0.0042685813811,
0.00238135509614,
0.00126624280575,
0.000643957680542,
0.000313312284132,
0.000145563772599,
6.45017277731e-05,
2.73432971552e-05,
1.1242475393e-05,
4.66091032876e-06,
2.12192702606e-06,
1.19436640045e-06,
8.7124597584e-07,
7.62851560055e-07,
7.27278616135e-07,
7.15086085434e-07,
7.09204930279e-07,
7.0301548247e-07,
6.93785503818e-07,
6.80860353766e-07,
6.63703072683e-07
};

double mc2[50]={
0.000829312892165,
0.00124276115093,
0.00339329172857,
0.0040822471492,
0.00383036583662,
0.00659159291536,
0.00816022697836,
0.00943640805781,
0.0137777375057,
0.0170593913645,
0.0213193036616,
0.0247343182564,
0.0280848778784,
0.0323308482766,
0.0370394326746,
0.0456917732954,
0.055876288563,
0.0576956197619,
0.0625325292349,
0.059160374105,
0.0656650811434,
0.0678329020739,
0.0625142157078,
0.0548068434,
0.0503893308342,
0.0402098186314,
0.0374446995556,
0.0299661569297,
0.027202475816,
0.0219328403473,
0.0179586578161,
0.0142926732078,
0.00839941669255,
0.00522366398945,
0.00224457983859,
0.000779274967499,
0.000197066590772,
7.16031790944e-05,
0.0,
0.0,
0.0,
0.0,
0.0,
0.0,
0.0,
0.0,
0.0,
0.0,
0.0,
0.0
    };

  if(nPU<=0) return 1;
  if(nPU>50) return 1;
  if (mc2[nPU+1]>0) {
      return data2[nPU-1]/mc2[nPU-1];
  }
  else {
      return 1;
  }
}
float VHbbAnalysis::puWeight_2016(int i){
double puw[75]={0.366077,
                0.893925,
                1.197716,
                0.962700,
                1.120976,
                1.164859,
                0.795599,
                0.495824,
                0.742182,
                0.878856,
                0.964232,
                1.072499,
                1.125336,
                1.176027,
                1.202083,
                1.207644,
                1.200176,
                1.182682,
                1.143999,
                1.096633,
                1.065602,
                1.051166,
                1.051600,
                1.050630,
                1.049862,
                1.058173,
                1.072155,
                1.083030,
                1.095694,
                1.107871,
                1.094621,
                1.082620,
                1.041247,
                0.985752,
                0.910807,
                0.820923,
                0.716787,
                0.610013,
                0.503118,
                0.404841,
                0.309195,
                0.227920,
                0.163690,
                0.113180,
                0.077300,
                0.050922,
                0.031894,
                0.020094,
                0.012263,
                0.007426,
                0.004380,
                0.002608,
                0.001566,
                0.000971,
                0.000729,
                0.000673,
                0.000730,
                0.000949,
                0.001355,
                0.001894,
                0.003082,
                0.004097,
                0.004874,
                0.005256,
                0.005785,
                0.005515,
                0.005000,
                0.004410,
                0.004012,
                0.003548,
                0.003108,
                0.002702,
                0.002337,
                0.002025,
                0.001723,
};
if (i < 0) return 1.;
if (i > 74) return puw[74];

return puw[i];
}
float VHbbAnalysis::puWeight_2016Up(int i){
double puw[75]={0.356705,
                0.703864,
                1.133091,
                0.845953,
                1.019643,
                1.049070,
                0.725717,
                0.347752,
                0.500571,
                0.603038,
                0.632187,
                0.732831,
                0.827762,
                0.912492,
                0.959905,
                0.988468,
                1.024407,
                1.052811,
                1.051032,
                1.027287,
                1.005587,
                0.997986,
                1.014926,
                1.037522,
                1.058117,
                1.085500,
                1.120673,
                1.155086,
                1.192505,
                1.231453,
                1.246278,
                1.267854,
                1.259038,
                1.233748,
                1.181511,
                1.104795,
                1.002414,
                0.889422,
                0.769089,
                0.653879,
                0.532514,
                0.422607,
                0.329814,
                0.249947,
                0.188563,
                0.138161,
                0.096842,
                0.068648,
                0.047345,
                0.032491,
                0.021707,
                0.014530,
                0.009585,
                0.006168,
                0.004278,
                0.003071,
                0.002242,
                0.001923,
                0.002023,
                0.002397,
                0.003652,
                0.004779,
                0.005725,
                0.006274,
                0.007045,
                0.006866,
                0.006372,
                0.005755,
                0.005366,
                0.004866,
                0.004374,
                0.003905,
                0.003469,
                0.003090,
                0.002704,
};
if (i < 0) return 1.;
if (i > 74) return puw[74];

return puw[i];
}

float VHbbAnalysis::puWeight_2016Down(int i){
double puw[75]={0.379279,
                1.141192,
                1.259883,
                1.098912,
                1.250071,
                1.280857,
                0.920153,
                0.767715,
                1.092593,
                1.337441,
                1.486275,
                1.528307,
                1.497815,
                1.500955,
                1.497305,
                1.443719,
                1.367840,
                1.298948,
                1.227337,
                1.165694,
                1.125431,
                1.090590,
                1.064057,
                1.039956,
                1.019238,
                1.006031,
                0.996969,
                0.984879,
                0.972839,
                0.956550,
                0.914665,
                0.872225,
                0.807097,
                0.734251,
                0.651014,
                0.561384,
                0.466382,
                0.374546,
                0.288546,
                0.214509,
                0.149742,
                0.099889,
                0.064343,
                0.039586,
                0.023889,
                0.013821,
                0.007566,
                0.004156,
                0.002215,
                0.001187,
                0.000643,
                0.000384,
                0.000273,
                0.000243,
                0.000290,
                0.000396,
                0.000544,
                0.000783,
                0.001151,
                0.001604,
                0.002568,
                0.003340,
                0.003880,
                0.004079,
                0.004373,
                0.004058,
                0.003579,
                0.003068,
                0.002712,
                0.002327,
                0.001978,
                0.001667,
                0.001397,
                0.001172,
                0.000965,
};
if (i < 0) return 1.;
if (i > 74) return puw[74];

return puw[i];
}

// puWeight reweighting functions taken from https://github.com/GLP90/Xbb/blob/merge_silvio/interface/VHbbNameSpace.h
float VHbbAnalysis::puWeight_ichep(int i){

double puw[38]={0.00026954692859,
                0.00834201570744,
                0.0146551644559,
                0.0267593718187,
                0.045546866213,
                0.0351739785649,
                0.0389242914756,
                0.131305445658,
                0.337172065156,
                0.645651266938,
                0.908493559723,
                1.09739459449,
                1.26533979742,
                1.41586705341,
                1.53524109717,
                1.49032150282,
                1.39123249386,
                1.45737145535,
                1.38588035333,
                1.4464676134,
                1.23282560995,
                1.08387439504,
                1.02663484504,
                0.97464906611,
                0.829692105572,
                0.758597838706,
                0.554540190093,
                0.442016886186,
                0.291492210196,
                0.203297812168,
                0.131804681706,
                0.0834846669777,
                0.0680713812995,
                0.0497105409204,
                0.0496692836227,
                0.0581182475108,
                0.089209326104,
                0.0941178579122,
               };
//i = i - 1;
if (i < 0) return 1.;
if (i > 37) return puw[37];

return puw[i];

}

float VHbbAnalysis::puWeight_2016to2017(int i){
double puw[75]={0.351171,
                1.020525,
                0.956480,
                1.003314,
                0.914853,
                1.012825,
                0.643938,
                0.231796,
                0.178580,
                0.227180,
                0.217812,
                0.223676,
                0.218336,
                0.223980,
                0.241786,
                0.270820,
                0.318400,
                0.378066,
                0.432400,
                0.478011,
                0.520014,
                0.558275,
                0.601169,
                0.649210,
                0.708200,
                0.788037,
                0.888684,
                1.001937,
                1.129693,
                1.269717,
                1.394521,
                1.538883,
                1.660589,
                1.774132,
                1.862991,
                1.927538,
                1.959004,
                1.976609,
                1.981554,
                2.006673,
                2.022258,
                2.087826,
                2.247218,
                2.493959,
                2.908338,
                3.437898,
                4.003752,
                4.796440,
                5.631103,
                6.577854,
                7.446925,
                8.388868,
                9.230176,
                9.786471,
                10.952979,
                12.161871,
                12.581994,
                13.092740,
                13.610222,
                13.219221,
                14.640417,
                13.127559,
                10.518075,
                7.657630,
                5.728162,
                3.748675,
                2.364398,
                1.473255,
                0.963583,
                0.623496,
                0.406713,
                0.267636,
                0.177675,
                0.119636,
                0.079882
};
if (i < 0) return 1.;
if (i > 74) return puw[74];

return puw[i];
}

float VHbbAnalysis::puWeight_ichep_up(int i){

double puw[38]={0.000168728884,
                0.006518139648,
                0.012767671099,
                0.021849922433,
                0.041129948359,
                0.031318966078,
                0.031506395408,
                0.069139953987,
                0.201682923142,
                0.431573069513,
                0.676688315382,
                0.885061147726,
                1.02772732515,
                1.1496216224,
                1.26358384716,
                1.24734104964,
                1.19950685946,
                1.30515668683,
                1.28849179873,
                1.39604961327,
                1.24030942651,
                1.13881318563,
                1.12938083469,
                1.13608683795,
                1.04512825001,
                1.04965380671,
                0.848278572582,
                0.748408550686,
                0.548182746816,
                0.426699751655,
                0.308798449343,
                0.217654143858,
                0.197782762841,
                0.16222993513 ,
                0.183774325453,
                0.245159931575,
                0.426559360962,
                0.491832630157,
               };
//i = i - 1;
if (i < 0) return 1.;
if (i > 37) return puw[37];

return puw[i];

}

float VHbbAnalysis::puWeight_ichep_down(int i){

double puw[38]={0.000394948025924,
                0.010589291412,
                0.0168294422346,
                0.0324871095383,
                0.0512240404177,
                0.0392771251619,
                0.0585204632322,
                0.252037472998,
                0.543974217454,
                0.93479398165,
                1.17260170322,
                1.36170067318,
                1.5793616475,
                1.74509976454,
                1.86131835195,
                1.75449370357,
                1.57204882742,
                1.57954789714,
                1.44148253688,
                1.43725104611,
                1.16760012792,
                0.975000220575,
                0.8640144709,
                0.750072693636,
                0.574219048746,
                0.469714734739,
                0.306667353503,
                0.217186088114,
                0.126764763282,
                0.0784330110012,
                0.0451854756945,
                0.0252751549695,
                0.0180049218223,
                0.0113901195301,
                0.00987729209672,
                0.0103547199185,
                0.0158597867448,
                0.0210296659761,
               };

//i = i - 1;
if (i < 0) return 1.;
if (i > 37) return puw[37];

return puw[i];

}

// from https://twiki.cern.ch/twiki/bin/view/CMS/VHiggsBBCodeUtils#V_X_QCD_and_EWK_corrections
float VHbbAnalysis::ptWeightQCD(int nGenVbosons, float LHE_HT, int GenVbosons_pdgId){
    float SF = 1.;
    if (LHE_HT>100 && nGenVbosons==1){
        if (GenVbosons_pdgId == 23){ // Z
            SF =   ((LHE_HT>=100 && LHE_HT<200)*1.588 * ( 280.35 / (409.860000) ) + (LHE_HT>=200 && LHE_HT<400)*1.438 * ( 77.67 / ( 110.880000 )) + (LHE_HT>=400 && LHE_HT<600)*1.494 * (10.73 / (13.189 )) + (LHE_HT>=600)*1.139 * ( 4.116 / (4.524300) ));
        }
        if (abs(GenVbosons_pdgId) == 24){
            SF =   ((LHE_HT>=100 && LHE_HT<200)*1.588 * ( 1345 / (1.23 *  1.29e3) ) + (LHE_HT>=200 && LHE_HT<400)*1.438 * ( 359.7 / ( 1.23 *  3.86e2)) + (LHE_HT>=400 && LHE_HT<600)*1.494 * (48.91 / (1.23 * 47.9 )) + (LHE_HT>=600)*1.139 * ( 18.77 / (1.23 * 19.9) ));
        }
    }
    return SF>0?SF:0;
}

// weights correction for EWK NLO correction
// from https://twiki.cern.ch/twiki/bin/view/CMS/VHiggsBBCodeUtils#V_X_QCD_and_EWK_corrections
float VHbbAnalysis::ptWeightEWK(int nGenVbosons,float GenVbosons_pt,int Vtype,int GenVbosons_pdgId){
    float SF = 1.;
    if (nGenVbosons ==1){
        if (Vtype == 0 || Vtype == 1 || Vtype == 4 || Vtype == 5){
            if (GenVbosons_pdgId == 23){
                //for Z options
                if (GenVbosons_pt > 100. && GenVbosons_pt < 3000) SF = -0.1808051+6.04146*(TMath::Power((GenVbosons_pt+759.098),-0.242556));
            }
        } else if (Vtype == 2 || Vtype == 3){
            //for W options
            if (GenVbosons_pdgId == 24 || GenVbosons_pdgId == -24){
                if (GenVbosons_pt > 100. && GenVbosons_pt < 3000) SF = -0.830041+7.93714*(TMath::Power((GenVbosons_pt+877.978),-0.213831));
            }
        }
    }
    return SF>0?SF:0;
}

// Copied from http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/CMSSW/TopQuarkAnalysis/SingleTop/src/TopProducer.cc?revision=1.9&view=markup
TLorentzVector VHbbAnalysis::getNu4Momentum(const TLorentzVector& TLepton, const TLorentzVector& TMET){
    TLorentzVector Lepton;
    Lepton.SetPxPyPzE(TLepton.Px(), TLepton.Py(), TLepton.Pz(), TLepton.E());
    TLorentzVector MET;
    MET.SetPxPyPzE(TMET.Px(), TMET.Py(), 0., TMET.E());

    double  mW = 80.38;

    //std::vector<math::XYZTLorentzVector> result;
    std::vector<TLorentzVector> result;

    //  double Wmt = sqrt(pow(Lepton.et()+MET.pt(),2) - pow(Lepton.px()+MET.px(),2) - pow(Lepton.py()+MET.py(),2) );

    double MisET2 = (MET.Px()*MET.Px() + MET.Py()*MET.Py());
    double mu = (mW*mW)/2 + MET.Px()*Lepton.Px() + MET.Py()*Lepton.Py();
    double a  = (mu*Lepton.Pz())/(Lepton.Energy()*Lepton.Energy() - Lepton.Pz()*Lepton.Pz());
    double a2 = TMath::Power(a,2);
    double b  = (TMath::Power(Lepton.Energy(),2.)*(MisET2) - TMath::Power(mu,2.))/(TMath::Power(Lepton.Energy(),2) - TMath::Power(Lepton.Pz(),2));
    double pz1(0),pz2(0),pznu(0);
    // int nNuSol(0);

    //math::XYZTLorentzVector p4nu_rec;
    TLorentzVector p4nu_rec;
    TLorentzVector p4W_rec;
    TLorentzVector p4b_rec;
    TLorentzVector p4Top_rec;
    TLorentzVector p4lep_rec;

    p4lep_rec.SetPxPyPzE(Lepton.Px(),Lepton.Py(),Lepton.Pz(),Lepton.Energy());

    //math::XYZTLorentzVector p40_rec(0,0,0,0);

    if(a2-b > 0 ){
        //if(!usePositiveDeltaSolutions_)
        //  {
        //    result.push_back(p40_rec);
        //    return result;
        //  }
        double root = sqrt(a2-b);
        pz1 = a + root;
        pz2 = a - root;
        //nNuSol = 2;

        //if(usePzPlusSolutions_)pznu = pz1;
        //if(usePzMinusSolutions_)pznu = pz2;
        //if(usePzAbsValMinimumSolutions_){
          pznu = pz1;
          if(fabs(pz1)>fabs(pz2)) pznu = pz2;
        //}

        double Enu = sqrt(MisET2 + pznu*pznu);

        p4nu_rec.SetPxPyPzE(MET.Px(), MET.Py(), pznu, Enu);

        result.push_back(p4nu_rec);
    }else{
        //if(!useNegativeDeltaSolutions_){
        //  result.push_back(p40_rec);
        //  return result;
        //}
        //    double xprime = sqrt(mW;

        double ptlep = Lepton.Pt(),pxlep=Lepton.Px(),pylep=Lepton.Py(),metpx=MET.Px(),metpy=MET.Py();

        double EquationA = 1;
        double EquationB = -3*pylep*mW/(ptlep);
        double EquationC = mW*mW*(2*pylep*pylep)/(ptlep*ptlep)+mW*mW-4*pxlep*pxlep*pxlep*metpx/(ptlep*ptlep)-4*pxlep*pxlep*pylep*metpy/(ptlep*ptlep);
        double EquationD = 4*pxlep*pxlep*mW*metpy/(ptlep)-pylep*mW*mW*mW/ptlep;

        std::vector<long double> solutions = EquationSolve<long double>((long double)EquationA,(long double)EquationB,(long double)EquationC,(long double)EquationD);

        std::vector<long double> solutions2 = EquationSolve<long double>((long double)EquationA,-(long double)EquationB,(long double)EquationC,-(long double)EquationD);

        double deltaMin = 14000*14000;
        double zeroValue = -mW*mW/(4*pxlep);
        double minPx=0;
        double minPy=0;

        //std::cout<<"a "<<EquationA << " b " << EquationB  <<" c "<< EquationC <<" d "<< EquationD << std::endl;

        //if(usePxMinusSolutions_){
        for( int i =0; i< (int)solutions.size();++i){
            if(solutions[i]<0 ) continue;
            double p_x = (solutions[i]*solutions[i]-mW*mW)/(4*pxlep);
            double p_y = ( mW*mW*pylep + 2*pxlep*pylep*p_x -mW*ptlep*solutions[i])/(2*pxlep*pxlep);
            double Delta2 = (p_x-metpx)*(p_x-metpx)+(p_y-metpy)*(p_y-metpy);

            //std::cout<<"intermediate solution1 met x "<<metpx << " min px " << p_x  <<" met y "<<metpy <<" min py "<< p_y << std::endl;

            if(Delta2< deltaMin && Delta2 > 0){deltaMin = Delta2;
            minPx=p_x;
            minPy=p_y;}
        }
        //}
        //std::cout<<"a "<<EquationA << " b " << EquationB  <<" c "<< EquationC <<" d "<< EquationD << std::endl;

        //if(usePxMinusSolutions_){
//          for( int i =0; i< (int)solutions.size();++i){
//          if(solutions[i]<0 ) continue;
//          double p_x = (solutions[i]*solutions[i]-mW*mW)/(4*pxlep);
//          double p_y = ( mW*mW*pylep + 2*pxlep*pylep*p_x -mW*ptlep*solutions[i])/(2*pxlep*pxlep);
//          double Delta2 = (p_x-metpx)*(p_x-metpx)+(p_y-metpy)*(p_y-metpy);
//
//                //std::cout<<"intermediate solution1 met x "<<metpx << " min px " << p_x  <<" met y "<<metpy <<" min py "<< p_y << std::endl;
//
//               //std::cout<<"solution1 met x "<<metpx << " min px " << minPx  <<" met y "<<metpy <<" min py "<< minPy << std::endl;
//          }
        //}


        //if(usePxPlusSolutions_){
        for( int i =0; i< (int)solutions2.size();++i){
            if(solutions2[i]<0 ) continue;
            double p_x = (solutions2[i]*solutions2[i]-mW*mW)/(4*pxlep);
            double p_y = ( mW*mW*pylep + 2*pxlep*pylep*p_x +mW*ptlep*solutions2[i])/(2*pxlep*pxlep);
            double Delta2 = (p_x-metpx)*(p_x-metpx)+(p_y-metpy)*(p_y-metpy);
            //std::cout<<"intermediate solution2 met x "<<metpx << " min px " << minPx  <<" met y "<<metpy <<" min py "<< minPy << std::endl;
            if(Delta2< deltaMin && Delta2 > 0){deltaMin = Delta2;
                minPx=p_x;
                minPy=p_y;
            }
                  //std::cout<<"solution2 met x "<<metpx << " min px " << minPx  <<" met y "<<metpy <<" min py "<< minPy << std::endl;
        }
        //}

        double pyZeroValue= ( mW*mW*pxlep + 2*pxlep*pylep*zeroValue);
        double delta2ZeroValue= (zeroValue-metpx)*(zeroValue-metpx) + (pyZeroValue-metpy)*(pyZeroValue-metpy);

        if(deltaMin==14000*14000) return TLorentzVector(0,0,0,0);
        //if(deltaMin==14000*14000) return result.front();
        //    else std::cout << " test " << std::endl;

        if(delta2ZeroValue < deltaMin){
          deltaMin = delta2ZeroValue;
          minPx=zeroValue;
          minPy=pyZeroValue;}


        double mu_Minimum = (mW*mW)/2 + minPx*pxlep + minPy*pylep;
        double a_Minimum  = (mu_Minimum*Lepton.Pz())/(Lepton.Energy()*Lepton.Energy() - Lepton.Pz()*Lepton.Pz());
        pznu = a_Minimum;

        //if(!useMetForNegativeSolutions_){
        double Enu = sqrt(minPx*minPx+minPy*minPy + pznu*pznu);
        p4nu_rec.SetPxPyPzE(minPx, minPy, pznu , Enu);
        //}
        //else{
        //  pznu = a;
        //  double Enu = sqrt(metpx*metpx+metpy*metpy + pznu*pznu);
        //  p4nu_rec.SetPxPyPzE(metpx, metpy, pznu , Enu);
        //}
        result.push_back(p4nu_rec);
    }
    return result.front();
}

// W-jet NLO to LO re-weighting function from Z(ll)H(bb)
double VHbbAnalysis::LOtoNLOWeightBjetSplitEtabb(double etabb, int njets){
    double SF = 1.;
    if(etabb < 5){
        if(njets < 1){
            SF =   0.935422 + 0.0403162*etabb -0.0089026*etabb*etabb +0.0064324*etabb*etabb*etabb -0.000212443*etabb*etabb*etabb*etabb;
        }else if(njets == 1){
            SF =   0.962415 +0.0329463*etabb -0.0414479*etabb*etabb +0.0240993*etabb*etabb*etabb -0.00278271*etabb*etabb*etabb*etabb;
        }else if(njets >=2){
            SF =   (0.721265 -0.105643*etabb -0.0206835*etabb*etabb +0.00558626*etabb*etabb*etabb)*TMath::Exp(0.450244*etabb);
        }
    }
    return SF;
}

float VHbbAnalysis::getVPtCorrFactor(float V_pt, int sn, int sysVar) {
    // get reco pT(W) linear correction obtained from fit to 2016 data.
    // Separate corrections derived for TT, W+LF, and combination W+HF + single top
    // sysVar: -1 for down, 0 nominal, 1 up
    float yint = 1.0;
    float slope = 0.0;
    if (cursample->sampleName.find("TT_")!=std::string::npos) {
        // ttbar
        if (m("dataYear") == 2016) {
            yint = 1.064;
            if (sysVar==0) {
                // nominal
                slope = -0.000380;
            }
            else if (sysVar==1) {
                slope = -0.000469;
            }
            else if (sysVar==-1) {
                slope = -0.000291;
            }
        }
        else if (m("dataYear") == 2017) {
            yint = 1.103;
            if (sysVar==0) {
                // nominal
                slope = -6.05429e-04;
            }
            else if (sysVar==1) {
                slope = -6.05429e-04 - 8.01886e-05;
            }
            else if (sysVar==-1) {
                slope = -6.05429e-04 + 8.01886e-05;
            }
        }
    }
    else if ((mInt("isWmunu")==1||mInt("isWenu")==1)&&cursample->sampleNum%100==0&&(cursample->sampleName.find("WJets-HT")!=std::string::npos || cursample->sampleName.find("WJets_madgraph")!=std::string::npos || cursample->sampleName.find("WBJets_")!=std::string::npos || cursample->sampleName.find("WJets_BGenFilter")!=std::string::npos)){
        // W+LF
        if (m("dataYear") == 2016) {
            yint = 1.097;
            if (sysVar==0) {
                // nominal
                slope = -0.000575;
            }
            else if (sysVar==1) {
                slope = -0.000621;
            }
            else if (sysVar==-1) {
                slope = -0.000529;
            }
        }
        else if (m("dataYear") == 2017) {
            yint = 1.115;
            if (sysVar==0) {
                // nominal
                slope = -6.36789e-04;
            }
            else if (sysVar==1) {
                slope = -6.36789e-04 - 3.87586e-05;
            }
            else if (sysVar==-1) {
                slope = -6.36789e-04 + 3.87586e-05;
            }
        }
    }
    else if ((mInt("isWmunu")==1||mInt("isWenu")==1)&&(cursample->sampleNum%100==1||cursample->sampleNum%100==2)&&(cursample->sampleName.find("WJets-HT")!=std::string::npos || cursample->sampleName.find("WJets_madgraph")!=std::string::npos || cursample->sampleName.find("WBJets_")!=std::string::npos || cursample->sampleName.find("WJets_BGenFilter")!=std::string::npos || (cursample->sampleNum>=16&&cursample->sampleNum<=21))){
        // W+HF + single top
        if (m("dataYear") == 2016) {
            yint = 1.259;
            if (sysVar==0) {
                // nominal
                slope = -0.00167;
            }
            else if (sysVar==1) {
                slope = -0.00180;
            }
            else if (sysVar==-1) {
                slope = -0.00154;
            }
        }
        if (m("dataYear") == 2017) {
            yint = 1.337;
            if (sysVar==0) {
                // nominal
                slope = -1.56131e-03;
            }
            else if (sysVar==1) {
                slope = -1.56131e-03 - 1.45980e-04;
            }
            else if (sysVar==-1) {
                slope = -1.56131e-03 + 1.45980e-04;
            }
        }
    }
    return (yint + V_pt*slope);
}


void VHbbAnalysis::smearJets(float JERScale) {
    for (int i=0; i<mInt("nJet"); i++) {
        float corr_nominal = mInt("Jet_corr_JER",i);
        if (corr_nominal == -99) { continue; }
        //std::cout<<"original smearing correction: "<<corr_nominal<<std::endl;
        f["Jet_corr_JER"][i] = 1 + JERScale * (m("Jet_corr_JER",i)-1);
        f["Jet_corr_JERUp"][i] = 1 + JERScale * (m("Jet_corr_JERUp",i)-1);
        f["Jet_corr_JERDown"][i] = 1 + JERScale * (m("Jet_corr_JERDown",i)-1);
        //std::cout<<"smearing jet "<<JERScale<<" times the nominal JER by factor "<<(f["Jet_corr_JER"][i] / corr_nominal)<<std::endl;
        //std::cout<<"Jet_pt = "<<f["Jet_pt"][i]<<std::endl;
        f["Jet_Pt"][i] = m("Jet_Pt",i) * ( m("Jet_corr_JER",i) / corr_nominal);
        // std::cout<<"Re-smeared jet_pt = "<<f["Jet_pt"][i]<<std::endl;

        // re-evaluate jet energy regression on re-smeared jets
        //std::cout<<"Jet_pt_reg was: "<<f["Jet_pt_reg"][i]<<std::endl;
        f["Jet_pt_reg_Heppy"][i] = m("Jet_bReg",i);
        f["Jet_bReg"][i] = evaluateRegression(i);
        //std::cout<<"now it has been re-evaluated to: "<<f["Jet_pt_reg"][i]<<std::endl;
    }
}

float VHbbAnalysis::evaluateRegression(int i) {
    if (m("Jet_bReg",i) == -99) { return -99; }
    *f["hJets_pt_0"] = m("Jet_Pt",i);
    *f["hJets_eta_0"] = m("Jet_eta",i);
    TLorentzVector tmp;
    tmp.SetPtEtaPhiM(m("Jet_Pt",i),m("Jet_eta",i),m("Jet_phi",i),m("Jet_mass",i));
    *f["hJets_mt_0"] = tmp.Mt();
    //std::cout<<"4-vector Mt() is "<<tmp.Mt()<<std::endl;
    //float mt = TMath::Sqrt( TMath::Power(tmp.Et(),2) - TMath::Power(tmp.Pt(),2) );
    //std::cout<<"by-hand Mt is "<<mt<<std::endl;
    //*f["hJets_mt_0"] = mt;
    *f["hJets_leadTrackPt_0"] = m("Jet_leadTrackPt",i);

    if (m("Jet_leptonPtRel",i) > 0) {
        *f["hJets_leptonPtRel_0"] = m("Jet_leptonPtRel",i);
    } else {
        *f["hJets_leptonPtRel_0"] = 0.;
    }

    if (m("Jet_leptonPt",i) > 0) {
        *f["hJets_leptonPt_0"] = m("Jet_leptonPt",i);
    } else {
        *f["hJets_leptonPt_0"] =  0.;
    }
    if (m("Jet_leptonDeltaR",i) > 0) {
        *f["hJets_leptonDeltaR_0"] = m("Jet_leptonDeltaR",i);
    } else {
        *f["hJets_leptonDeltaR_0"] = 0.;
    }

    *f["hJets_neHEF_0"] = m("Jet_neHEF",i);
    *f["hJets_neEmEF_0"] = m("Jet_neEmEF",i);
    *f["hJets_vtxMass_0"] = m("Jet_vtxMass",i);
    *f["hJets_vtxPt_0"] = m("Jet_vtxPt",i);

    if (m("Jet_vtx3DVal",i) > 0) {
        *f["hJets_vtx3dL_0"] = m("Jet_vtx3DVal",i);
    }else {
        *f["hJets_vtx3dL_0"] = 0.;
    }
    *f["hJets_vtxNtracks_0"] = m("Jet_vtxNtracks",i);
    //*f["hJets_vtx3deL_0"] = f["Jet_vtx3DSig"][i];
    if (m("Jet_vtx3DSig",i) > 0) {
        *f["hJets_vtx3deL_0"] = m("Jet_vtx3DVal",i) / m("Jet_vtx3DSig",i);
    }else {
        *f["hJets_vtx3deL_0"] = 0;
    }
    return EvaluateRegression(bdtInfos["bdt_bjetreg"]);

}

void VHbbAnalysis::SetupFactorizedJECs(std::string variation) {
    std::string corrsToCalc[] = { "AbsoluteStatUp", "AbsoluteStatDown", "AbsoluteScaleUp", "AbsoluteScaleDown", "AbsoluteFlavMapUp", "AbsoluteFlavMapDown", "AbsoluteMPFBiasUp", "AbsoluteMPFBiasDown", "FragmentationUp", "FragmentationDown", "SinglePionECALUp", "SinglePionECALDown", "SinglePionHCALUp", "SinglePionHCALDown", "FlavorQCDUp", "FlavorQCDDown", "TimePtEtaUp", "TimePtEtaDown", "RelativeJEREC1Up", "RelativeJEREC1Down", "RelativeJEREC2Up", "RelativeJEREC2Down", "RelativeJERHFUp", "RelativeJERHFDown", "RelativePtBBUp", "RelativePtBBDown", "RelativePtEC1Up", "RelativePtEC1Down", "RelativePtEC2Up", "RelativePtEC2Down", "RelativePtHFUp", "RelativePtHFDown", "RelativeBalUp", "RelativeBalDown", "RelativeFSRUp", "RelativeFSRDown", "RelativeStatFSRUp", "RelativeStatFSRDown", "RelativeStatECUp", "RelativeStatECDown", "RelativeStatHFUp", "RelativeStatHFDown", "PileUpDataMCUp", "PileUpDataMCDown", "PileUpPtRefUp", "PileUpPtRefDown", "PileUpPtBBUp", "PileUpPtBBDown", "PileUpPtEC1Up", "PileUpPtEC1Down", "PileUpPtEC2Up", "PileUpPtEC2Down", "PileUpPtHFUp", "PileUpPtHFDown", "PileUpMuZeroUp", "PileUpMuZeroDown", "PileUpEnvelopeUp", "PileUpEnvelopeDown", "SubTotalPileUpUp", "SubTotalPileUpDown", "SubTotalRelativeUp", "SubTotalRelativeDown", "SubTotalPtUp", "SubTotalPtDown", "SubTotalScaleUp", "SubTotalScaleDown", "SubTotalAbsoluteUp", "SubTotalAbsoluteDown", "SubTotalMCUp", "SubTotalMCDown", "TotalUp", "TotalDown", "FlavorZJetUp", "FlavorZJetDown", "FlavorPhotonJetUp", "FlavorPhotonJetDown", "FlavorPureGluonUp", "FlavorPureGluonDown", "FlavorPureQuarkUp", "FlavorPureQuarkDown", "FlavorPureCharmUp", "FlavorPureCharmDown", "FlavorPureBottomUp", "FlavorPureBottomDown", "TimeRunBCDUp", "TimeRunBCDDown", "TimeRunEFUp", "TimeRunEFDown", "TimeRunGUp", "TimeRunGDown", "TimeRunHUp", "TimeRunHDown", "CorrelationGroupMPFInSituUp", "CorrelationGroupMPFInSituDown", "CorrelationGroupIntercalibrationUp", "CorrelationGroupIntercalibrationDown", "CorrelationGroupbJESUp", "CorrelationGroupbJESDown", "CorrelationGroupFlavorUp", "CorrelationGroupFlavorDown", "CorrelationGroupUncorrelatedUp", "CorrelationGroupUncorrelatedDown" };

    bool foundVar = false;
    for (int i=0; i < 102; i++) {
        if (variation == corrsToCalc[i].c_str()) foundVar = true;
    }

    if (!foundVar) return;

    for (int j=0; j < mInt("nJet"); j++) {
        float Jet_corr_var = m("Jet_corr_" + variation,j);
        float corr_nominal = m("Jet_corr",j);
        float scaleShift = Jet_corr_var / corr_nominal;
        f["Jet_corr_"+variation+"_ratio"][j] = scaleShift;
        float Jet_pt_nom = m("Jet_Pt",j);
        float Jet_pt_reg_nom = evaluateRegression(j); // probably safer to re-evaluate the nominal in case there is any residual bias in re-implementating of reg.
        f["Jet_Pt"][j] = m("Jet_Pt",j) * scaleShift;
        float Jet_pt_reg_var = evaluateRegression(j); // re-evaluate regression on top of variation
        f["Jet_Pt"][j] = Jet_pt_nom;
        f["Jet_pt_reg_corr"+variation+"_ratio"][j] = Jet_pt_reg_var / Jet_pt_reg_nom;
    }

}

double VHbbAnalysis::computeEventSFForDoubleLeptonTrig(std::string dataeff_lowptleg, std::string dataeff_highptleg, std::string mceff_lowptleg, std::string mceff_highptleg) {
    double dataeff_lowptleg_1 = m(dataeff_lowptleg,mInt("lepInd1"));
	double dataeff_lowptleg_2 = m(dataeff_lowptleg,mInt("lepInd2"));
	double dataeff_highptleg_1 = m(dataeff_highptleg,mInt("lepInd1"));
	double dataeff_highptleg_2 = m(dataeff_highptleg,mInt("lepInd2"));
	double mceff_lowptleg_1 = m(mceff_lowptleg,mInt("lepInd1"));
	double mceff_lowptleg_2 = m(mceff_lowptleg,mInt("lepInd2"));
	double mceff_highptleg_1 = m(mceff_highptleg,mInt("lepInd1"));
	double mceff_highptleg_2 = m(mceff_highptleg,mInt("lepInd2"));

	double effData_ = (std::pow(dataeff_lowptleg_1,2)*dataeff_highptleg_2 + std::pow(dataeff_lowptleg_2,2)*dataeff_highptleg_1)/(dataeff_lowptleg_1+dataeff_lowptleg_2);
	double effMC_ = (std::pow(mceff_lowptleg_1,2)*mceff_highptleg_2 + std::pow(mceff_lowptleg_2,2)*mceff_highptleg_1)/(mceff_lowptleg_1+mceff_lowptleg_2);
	double eff_ = effData_/effMC_;
    return eff_;
}
