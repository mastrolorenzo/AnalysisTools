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
    return;
}

// Check if events pass preselection.
// Returns true for events passing the preselection and false otherwise.
bool VHbbAnalysis::Preselection() {
    bool doCutFlowInPresel = int(*f["doCutFlow"]) < 0;

    if (*f["onlyEvenEvents"] && *f["onlyOddEvents"]) {
        std::cout << "Cannot set both onlyEvenEvents and onlyOddEvents to true!!" << std::endl;
        return false;
    }

    if (*f["onlyEvenEvents"]) {
        if (*in["event"] % 2 == 1) return false;
    } else if (*f["onlyOddEvents"]) {
        if (*in["event"] % 2 == 0) return false;
    }

    // stitch WJets inclusive sample to HT-binned samples
    if (cursample->sampleNum == 22 && *f["LHE_HT"] > 100) return false;
    // stitch ZJets inclusive sample to HT-binned samples
    if (cursample->sampleNum == 23 && *f["LHE_HT"] > 100) return false;

    // FIXME missing nGenStatus2bHad (must add gen particle loop--eventually push to post-processor)
    // use W+jets b-enriched samples but make sure all samples are orthogonal
    if (   cursample->sampleNum == 22 || cursample->sampleNum == 41 || cursample->sampleNum == 42
        || cursample->sampleNum == 43 || cursample->sampleNum == 44 || cursample->sampleNum == 45
        || cursample->sampleNum == 46 || cursample->sampleNum == 47
       ) {
        if (*f["LHE_Vpt"] > 100) {
            if (*in["LHE_Nb"] != 0) return false;
            //if (*in["LHE_Nb"] != 0 || *in["nGenStatus2bHad"] != 0) return false;
        }
    } else if (cursample->sampleNum == 48) {
        if (*f["LHE_Vpt"] < 100 || *f["LHE_Vpt"] > 200 || *in["LHE_Nb"] == 0) return false;
    } else if (cursample->sampleNum == 481) {
        if (*f["LHE_Vpt"] < 200 || *in["LHE_Nb"] == 0) return false;
    } else if (cursample->sampleNum == 49) {
        if (*f["LHE_Vpt"] < 100 || *f["LHE_Vpt"] > 200 ) return false;
        //if (*f["LHE_Vpt"] < 100 || *f["LHE_Vpt"] > 200 || *in["nGenStatus2bHad"] == 0) return false;
    } else if (cursample->sampleNum == 491) {
        if (*f["LHE_Vpt"] < 200) return false;
        //if (*f["LHE_Vpt"] < 200 || *in["nGenStatus2bHad"] == 0) return false;
    }

    if (cursample->sampleNum == 0) {
        //if (*f["json"] != 1 && !doCutFlowInPresel) return false;
        if (!doCutFlowInPresel) return false;
    }

    // Trigger Selections
    if(*f["reVType"]){
        *in["Vtype"]= UpdatedVType();
    }

    if (!PassVTypeAndTrigger(*in["Vtype"]) && !doCutFlowInPresel) {
        *in["controlSample"] = -1;
    }

    if (*in["Vtype"]<2) {
        if (*f["V_pt"] < *f["vptPreselCut"] && !doCutFlowInPresel) { // preserve 50 < Vpt < 100 for the 2-lepton channels
            return false;
        }
    } else if (*f["V_pt"] < *f["vptcut"] && !doCutFlowInPresel) {    // else cut away Vpt < 100
        return false;
    }

    // Heppy jet corrections for JER/JEC are full correction, it's easier to just use the
    // correction on top of the nominal

    // FIXME configure b-tagged via config?
    //    if (int(*f["doCMVA"]) != 0) {
    //        // use CMVAv2 discriminator instead of CSV
    //        f["Jet_btagCSVV2"][i] = f["Jet_btagCMVA"][i];
    //    }


    if (*f["smearJets"]) {
        for (int i = 0; i < *in["nJet"]; i++) {
            float JERScale = *f["JERScale"]; // apply JER smearing x times the nominal smearing amount
            if (JERScale != 1.0) {
                smearJets(JERScale);
            }

            // not sure about these names -- won't do just now FIXME
            //f["Jet_corr_JERUp_ratio"][i] = f["Jet_corr_JERUp"][i] / f["Jet_corr_JER"][i] ;
            //f["Jet_corr_JERDown_ratio"][i] = f["Jet_corr_JERDown"][i] / f["Jet_corr_JER"][i] ;
            //f["Jet_corr_JECUp_ratio"][i] = f["Jet_corr_JECUp"][i] / f["Jet_corr"][i] ;
            //f["Jet_corr_JECDown_ratio"][i] = f["Jet_corr_JECDown"][i] / f["Jet_corr"][i] ;
            //f["Jet_pt_reg_corrJECUp_ratio"][i] = f["Jet_pt_reg_corrJECUp"][i] / f["Jet_bReg"][i] ;
            //f["Jet_pt_reg_corrJECDown_ratio"][i] = f["Jet_pt_reg_corrJECDown"][i] / f["Jet_bReg"][i] ;
            //f["Jet_pt_reg_corrJERUp_ratio"][i] = f["Jet_pt_reg_corrJERUp"][i] / f["Jet_bReg"][i] ;
            //f["Jet_pt_reg_corrJERDown_ratio"][i] = f["Jet_pt_reg_corrJERDown"][i] / f["Jet_bReg"][i] ;
        }
    }

    if (*f["reReg"]) {
        for (int i = 0; i < *in["nJet"]; i++) {
            f["Jet_bReg"][i] = evaluateRegression(i);
        }
    }

    // add smearing on top of regression-->replace Jet_bReg
    for (int i = 0; i < *in["nJet"]; i++) {
        f["Jet_bReg"][i] = f["Jet_bReg"][i] * f["Jet_Pt"][i]/ f["Jet_pt"][i];
    }

    // FIXME why removed?
    //SetupFactorizedJECs(cursyst->name);

    // FIXME add Factorized JEC systs for MET
    //*f["MET_pt_JECUp_ratio"] = *f["met_shifted_JetEnUp_pt"] / *f["MET_pt"];
    //*f["MET_pt_JECDown_ratio"] = *f["met_shifted_JetEnDown_pt"] / *f["MET_pt"];

    // Preselect for two jets and one lepton which pass some minimum pt threshold
    int nPreselJets = 0;
    for (int i = 0; i < *in["nJet"]; i++) {
        if (f["Jet_bReg"][i] > *f["JetPtPresel"] && in["Jet_puId"][i] > 0 && b["Jet_lepFilter"][i]) nPreselJets++;
    }

    // int nPreselLep = 0; // FIXME needs to be configurable or possibly removed
    // for (int i = 0; i < *in["nselLeptons"]; i++) {
    //     if (f["selLeptons_pt"][i] > *f["LepPtPresel"]
    //         && fabs(f["selLeptons_eta"][i]) < 2.5
    //         && (fabs(in["selLeptons_pdgId"][i]) == 11 || fabs(in["selLeptons_pdgId"][i]) == 13)
    //        ) {
    //         nPreselLep++;
    //     }
    // }

    // remove nPreselLep cut for now since we want to pre-select for all channels
    // FIXME add boosted analyses and BDTs
    if (int(*f["doBoost"]) == 0) {
        if (nPreselJets < 2) return false;
    //    if (nPreselJets < 2 && !doCutFlowInPresel) return false;
    //} else {
    //    if (nPreselJets < 2 && *in["nFatjetAK08ungroomed"] < 1 && !doCutFlowInPresel) return false;
    }

    return true;
}

bool VHbbAnalysis::Analyze() {
    *in["sampleIndex"] = cursample->sampleNum;
    bool sel = true;
    bool doCutFlow = int(*f["doCutFlow"]) != 0;
    bool doOnlySignalRegion = bool(*f["doOnlySignalRegion"]);
    *in["controlSample"] = 0;
    *in["cutFlow"] = 0;

    // move fast requirements to the front
    if (debug > 1000) std::cout << "Imposing json and trigger requirements" << std::endl;
    if (*in["sampleIndex"] == 0) {
        *in["controlSample"] = -1;
    }

    // asking for doCutFlow instead of doOnlySignalRegion, as events with bad json will not go into CR
    if (!doCutFlow && *in["controlSample"] < 0) {
        return false;
    }

    if (!PassVTypeAndTrigger(*in["Vtype"])) {
        *in["controlSample"] = -1;
    }

    // asking for doCutFlow instead of doOnlySignalRegion, as events with bad trigger will not go into CR
    if (!doCutFlow && *in["controlSample"] < 0) {
        return false;
    }

    if (sel && *in["controlSample"] > -1) {
        *in["cutFlow"] += 1; // json and trigger selection
    }


//      _|                        _|
//      _|    _|_|    _|_|_|    _|_|_|_|    _|_|    _|_|_|      _|_|_|
//      _|  _|_|_|_|  _|    _|    _|      _|    _|  _|    _|  _|_|
//      _|  _|        _|    _|    _|      _|    _|  _|    _|      _|_|
//      _|    _|_|_|  _|_|_|        _|_|    _|_|    _|    _|  _|_|_|
//                    _|
//                    _|

    // leptons have to be selected before bjets, since we need to know which channel we're in
    *in["lepInd1"] = *in["muInd1"] = *in["elInd1"] = -1;
    *in["lepInd2"] = *in["muInd2"] = *in["elInd2"] = -1;
    *in["isZnn"] = 0;
    *in["isWmunu"] = 0;
    *in["isWenu"] = 0;
    *in["isZmm"] = 0;
    *in["isZee"] = 0;

    if (debug > 1000) std::cout << "Selecting leptons" << std::endl;


/*
    std::pair<int,int> good_muons_2lep = HighestPtGoodMuonsOppCharge(    *f["muptcut_2lepchan"], *f["murelisocut_2lepchan"]);
    std::pair<int,int> good_elecs_2lep = HighestPtGoodElectronsOppCharge(*f["elptcut_2lepchan"], *f["elrelisocut_2lepchan"]);
    std::pair<int,int> good_muons_1lep = HighestPtGoodMuonsOppCharge(    *f["muptcut_1lepchan"], *f["murelisocut_1lepchan"]);
    std::pair<int,int> good_elecs_1lep = HighestPtGoodElectronsOppCharge(*f["elptcut_1lepchan"], *f["elrelisocut_1lepchan"]);

    if (good_muons_2lep.second > -1 && (cursample->lepFlav == -1 || cursample->lepFlav == 0)) {
        *in["isZmm"] = 1;
        *in["lepInd1"] = *in["muInd1"] = good_muons_2lep.first;
        *in["lepInd2"] = *in["muInd2"] = good_muons_2lep.second;
    } else if (good_elecs_2lep.second > -1 && (cursample->lepFlav == -1 || cursample->lepFlav == 1)) {
        *in["isZee"] = 1;
        *in["lepInd1"] = *in["elInd1"] = good_elecs_2lep.first;
        *in["lepInd2"] = *in["elInd2"] = good_elecs_2lep.second;
    } else if (good_muons_1lep.first > -1 && (cursample->lepFlav == -1 || cursample->lepFlav == 0)) {
        *in["isWmunu"] = 1;
        *in["eventClass"] = 0;      // TODO is this still needed????
        *in["lepInd1"] = *in["muInd1"] = good_muons_1lep.first;
    } else if (good_elecs_1lep.first > -1 && (cursample->lepFlav == -1 || cursample->lepFlav == 1)) {
        *in["isWenu"] = 1;
        *in["eventClass"] = 1000;   // TODO is this still needed????
        *in["lepInd1"] = *in["elInd1"] = good_elecs_1lep.first;
    } else if (*f["MET_pt"] > *f["metcut_0lepchan"]
            && *in["Flag_goodVertices"]
            && *in["Flag_globalTightHalo2016Filter"]
            && *in["Flag_HBHENoiseFilter"]
            && *in["Flag_HBHENoiseIsoFilter"]
            && *in["Flag_EcalDeadCellTriggerPrimitiveFilter"]) {
        *in["isZnn"] = 1;
    } else {
        return false;  // TODO no channel selected => can return here  CORRECT??
    }
*/
    if (*in["Vtype"] == 0 && (cursample->lepFlav == -1 || cursample->lepFlav == 0)) {
        std::pair<int,int> good_muons_2lep = HighestPtGoodMuonsOppCharge(    m("muptcut_2lepchan"), m("murelisocut_2lepchan"));
        if (good_muons_2lep.second > -1) {
            *in["isZmm"] = 1;
            *in["lepInd1"] = *in["muInd1"] = good_muons_2lep.first;
            *in["lepInd2"] = *in["muInd2"] = good_muons_2lep.second;
        } else {
            *in["controlSample"] = -1;
        }
    } else if (*in["Vtype"] == 1 && (cursample->lepFlav == -1 || cursample->lepFlav == 1)) {
        std::pair<int,int> good_elecs_2lep = HighestPtGoodElectronsOppCharge(m("elptcut_2lepchan"), m("elrelisocut_2lepchan"), m("elidcut_2lepchan"));
        if (good_elecs_2lep.second > -1) {
            *in["isZee"] = 1;
            *in["lepInd1"] = *in["elInd1"] = good_elecs_2lep.first;
            *in["lepInd2"] = *in["elInd2"] = good_elecs_2lep.second;
        } else {
            *in["controlSample"] = -1;
        }
    } else if (*in["Vtype"] == 2 && (cursample->lepFlav == -1 || cursample->lepFlav == 0)) {
        std::pair<int,int> good_muons_1lep = HighestPtGoodMuonsOppCharge(    *f["muptcut_1lepchan"], *f["murelisocut_1lepchan"]);
        if (good_muons_1lep.first > -1) {
            *in["isWmunu"] = 1;
            *in["eventClass"] = 0;      // TODO is this still needed????
            *in["lepInd1"] = *in["muInd1"] = good_muons_1lep.first;
        } else {
            *in["controlSample"] = -1;
        }
    } else if (*in["Vtype"] == 3 && (cursample->lepFlav == -1 || cursample->lepFlav == 1)) {
        std::pair<int,int> good_elecs_1lep = HighestPtGoodElectronsOppCharge(*f["elptcut_1lepchan"], *f["elrelisocut_1lepchan"], *f["elidcut_1lepchan"]);
        if (good_elecs_1lep.first > -1) {
            *in["isWenu"] = 1;
            *in["eventClass"] = 1000;   // TODO is this still needed????
            *in["lepInd1"] = *in["elInd1"] = good_elecs_1lep.first;
        } else {
            *in["controlSample"] = -1;
        }
    } else if (*in["Vtype"] == 4) {
        if (*f["MET_pt"] > *f["metcut_0lepchan"]
            && *b["Flag_goodVertices"]
            && *b["Flag_globalTightHalo2016Filter"]
            && *b["Flag_HBHENoiseFilter"]
            && *b["Flag_HBHENoiseIsoFilter"]
            && *b["Flag_EcalDeadCellTriggerPrimitiveFilter"]) {
            *in["isZnn"] = 1;
        } else {
            *in["controlSample"] = -1;
        }
    } else {
        *in["controlSample"] = -1;
    }


    if (debug > 1000) {
        std::cout << "Vtype isZnn isWmunu isWenu isZmm isZee ";
        std::cout << *in["Vtype"] << " "
                  << *in["isZnn"] << " "
                  << *in["isWmunu"] << " "
                  << *in["isWenu"] << " "
                  << *in["isZmm"] << " "
                  << *in["isZee"]
                  << std::endl;
    }

    if (sel && *in["controlSample"] > -1) {
        *in["cutFlow"] += 1; // lepton selection
    }

    // asking for doCutFlow instead of doOnlySignalRegion, as events with bad leptons will not go into CR
    if (!doCutFlow && *in["controlSample"] < 0) {
        return false;
    }


//        _|              _|
//              _|_|    _|_|_|_|    _|_|_|
//        _|  _|_|_|_|    _|      _|_|
//        _|  _|          _|          _|_|
//        _|    _|_|_|      _|_|  _|_|_|
//        _|
//      _|


    // leptons had to be selected before bjets, since we need to know which channel we're in
    if (debug > 1000) std::cout << "Looping over jets" << std::endl;
    for (int i = 0; i < *in["nJet"]; i++) {
        if (int(*f["doReg"]) == 0) {
            // don't apply the regression. The easiest way to do this is to just reset the value of
            // the regressed jet pt to the nominal jet pt, since we build everything later from the
            // jet pt's
            f["Jet_bReg"][i] = f["Jet_Pt"][i];
        } else {
            // add smearing on top of regression-->replace Jet_bReg
            f["Jet_bReg"][i] = f["Jet_bReg"][i] * f["Jet_Pt"][i]/ f["Jet_pt"][i];
        }

        if (int(*f["doCMVA"]) != 0) {
           // use CMVAv2 discriminator instead of CSV
           f["Jet_btagCSVV2"][i] = f["Jet_btagCMVA"][i];
        }

        // Do this in pre-selection
        //// apply JER smearing x times the nominal smearing amount
        //float JERScale = *f["JERScale"];
        //if (JERScale != 1.0) {
        //    smearJets(JERScale);
        //} else if (int(*f["reReg"]) != 0) {
        //    f["Jet_bReg"][i] = evaluateRegression(i);
        //}
    }

    // keep track of each jet selection method separately
    float j1ptCut, j2ptCut, j1ptCSV;
    if (*in["isZnn"]) {
        j1ptCut = *f["j1ptCut_0lepchan"];
        j2ptCut = *f["j2ptCut_0lepchan"];
        j1ptCSV = *f["j1ptCSV_0lepchan"];
    } else if (*in["isWmunu"] || *in["isWenu"]) {
        j1ptCut = *f["j1ptCut_1lepchan"];
        j2ptCut = *f["j2ptCut_1lepchan"];
        j1ptCSV = *f["j1ptCSV_1lepchan"];
    } else if (*in["isZmm"] || *in["isZee"]) {
        j1ptCut = *f["j1ptCut_2lepchan"];
        j2ptCut = *f["j2ptCut_2lepchan"];
        j1ptCSV = *f["j1ptCSV_2lepchan"];
    }

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
    std::pair<int,int> bjets = HighestCSVBJets(j1ptCut, j2ptCut);

    // put CMVA cuts out of selection functions
    if (bjets.first != -1 && bjets.second != -1) {
        if (f["Jet_btagCMVA"][bjets.first] < j1ptCSV) {
            *in["controlSample"] = -1;
        } else if (f["Jet_btagCMVA"][bjets.second] < *f["j2ptCSV"]) {  // 2nd jet CSV is the same in all SR's
            *in["controlSample"] = -1;
        }
    }
    if (doOnlySignalRegion && *in["controlSample"] < 0) {
        return false;
    }

    if (bjets.first == -1) {
        if (int(*f["doBoost"]) == 0) {
            return false;  // TODO this must be the same for signal and control regions?
        }
        *in["hJetInd1"] = 0;
    } else {
        *in["hJetInd1"] = bjets.first;
    }

    if (bjets.second == -1) {
        if (int(*f["doBoost"]) == 0) {
            return false;  // TODO this must be the same for signal and control regions?
            *in["hJetInd2"] = 1;
        } else if (*in["nJet"] < 2) {
            // assume nJet (AK4 jets) > 0 even for boosted events.
            // Is this reasonable? If not code needs some overhaul to handle boosted analysis.
            *in["hJetInd2"] = 0;
        }
    } else {
        *in["hJetInd2"] = bjets.second;
    }

    if (int(*f["doBoost"]) != 0) {
        // need to do some filtering on b-tagging for events that we won't use in the boosted analysis
        // (i.e. there is no fat jet) or file size for W+jets is ridiculous for boosted analysis
        if (*ui["nFatJet"] < 1) {
            if (f["Jet_btagCMVA"][*in["hJetInd1"]] < *f["j1ptCSV"] || f["Jet_btagCMVA"][*in["hJetInd2"]] < *f["j2ptCSV"]) {
                return false;  // TODO this must be the same for signal and control regions?
            }
        }
    }

    if (sel && *in["controlSample"] > -1) {
        *in["cutFlow"] += 1; // selected jets
    }


//      _|        _|
//      _|_|_|          _|_|_|    _|_|_|    _|_|_|      _|  _|_|    _|_|      _|_|_|    _|_|
//      _|    _|  _|  _|    _|  _|    _|  _|_|          _|_|      _|_|_|_|  _|        _|    _|
//      _|    _|  _|  _|    _|  _|    _|      _|_|      _|        _|        _|        _|    _|
//      _|    _|  _|    _|_|_|    _|_|_|  _|_|_|        _|          _|_|_|    _|_|_|    _|_|
//                          _|        _|
//                      _|_|      _|_|


    if (debug > 1000) {
        std::cout << "nJet = " << *in["nJet"] << std::endl;
        std::cout << "hJetInd1 = " << *in["hJetInd1"] << std::endl;
        std::cout << "hJetInd2 = " << *in["hJetInd2"] << std::endl;
        std::cout << "found two bjets with pt and CSV "
                  << f["Jet_bReg"][*in["hJetInd1"]] << " "
                  << f["Jet_btagCSVV2"][*in["hJetInd1"]] << " "
                  << f["Jet_bReg"][*in["hJetInd2"]] << " "
                  << f["Jet_btagCSVV2"][*in["hJetInd2"]] << " "
                  << std::endl;
    }

    // Reconstruct Higgs
    TLorentzVector HJ1, HJ2, Hbb;
    HJ1.SetPtEtaPhiM(f["Jet_bReg"][*in["hJetInd1"]], f["Jet_eta"][*in["hJetInd1"]], f["Jet_phi"][*in["hJetInd1"]], f["Jet_mass"][*in["hJetInd1"]] * (f["Jet_bReg"][*in["hJetInd1"]] / f["Jet_pt"][*in["hJetInd1"]]));
    HJ2.SetPtEtaPhiM(f["Jet_bReg"][*in["hJetInd2"]], f["Jet_eta"][*in["hJetInd2"]], f["Jet_phi"][*in["hJetInd2"]], f["Jet_mass"][*in["hJetInd2"]] * (f["Jet_bReg"][*in["hJetInd2"]] / f["Jet_pt"][*in["hJetInd2"]]));
    Hbb = HJ1 + HJ2;

    TLorentzVector HJ1_noreg, HJ2_noreg, Hbb_noreg;
    HJ1_noreg.SetPtEtaPhiM(f["Jet_pt"][*in["hJetInd1"]], f["Jet_eta"][*in["hJetInd1"]], f["Jet_phi"][*in["hJetInd1"]], f["Jet_mass"][*in["hJetInd1"]]);
    HJ2_noreg.SetPtEtaPhiM(f["Jet_pt"][*in["hJetInd2"]], f["Jet_eta"][*in["hJetInd2"]], f["Jet_phi"][*in["hJetInd2"]], f["Jet_mass"][*in["hJetInd2"]]);
    Hbb_noreg = HJ1_noreg + HJ2_noreg;

    *f["H_pt"] = Hbb.Pt();
    m("isWmunu");
    m("isWenu");
    m("H_pt");
    m("hptcut_1lepchan");
    m("controlSample");
    if (*in["isZnn"] && *f["H_pt"] < *f["hptcut_0lepchan"]) {
        *in["controlSample"] = -1;
    } else if ((m("isWmunu") || m("isWenu")) && m("H_pt") < m("hptcut_1lepchan")) {
        *in["controlSample"] = -1;
    } else if (sel && *in["controlSample"] > -1) {
        *in["cutFlow"] += 1; // pT(jj) cut
    }

    if(debug>1000) std::cout<<"doOnlySignalRegion controlSample "<<doOnlySignalRegion<<" "<<*in["controlSample"]<<std::endl;
    if (doOnlySignalRegion && *in["controlSample"] < 0) {
        return false;
    }


//      _|      _|
//      _|      _|      _|  _|_|    _|_|      _|_|_|    _|_|
//      _|      _|      _|_|      _|_|_|_|  _|        _|    _|
//        _|  _|        _|        _|        _|        _|    _|
//          _|          _|          _|_|_|    _|_|_|    _|_|


    // channel specific selection and vector boson reconstruction
    TLorentzVector V, GenBJ1, GenBJ2, GenBJJ;
    TLorentzVector W_withNuFromMWCon;

    if (*in["isZee"] == 1 || *in["isZmm"] == 1) {
        TLorentzVector lep1, lep2;
        if (*in["isZee"] == 1){
            lep1.SetPtEtaPhiM(f["Electron_pt"][*in["lepInd1"]], f["Electron_eta"][*in["lepInd1"]], f["Electron_phi"][*in["lepInd1"]], f["Electron_mass"][*in["lepInd1"]]);
            lep2.SetPtEtaPhiM(f["Electron_pt"][*in["lepInd2"]], f["Electron_eta"][*in["lepInd2"]], f["Electron_phi"][*in["lepInd2"]], f["Electron_mass"][*in["lepInd2"]]);
        } else if (*in["isZmm"] == 1){
            lep1.SetPtEtaPhiM(f["Muon_pt"][*in["lepInd1"]], f["Muon_eta"][*in["lepInd1"]], f["Muon_phi"][*in["lepInd1"]], f["Muon_mass"][*in["lepInd1"]]);
            lep2.SetPtEtaPhiM(f["Muon_pt"][*in["lepInd2"]], f["Muon_eta"][*in["lepInd2"]], f["Muon_phi"][*in["lepInd2"]], f["Muon_mass"][*in["lepInd2"]]);
        } else {
            std::cout<<"This shouldn't happen.  isZee  isZmm "<<*in["isZee"]<<" "<<*in["isZmm"]<<std::endl;
        }
        V = lep1 + lep2;
    } else if (*in["isWenu"] == 1 || *in["isWmunu"] == 1) {
        if (*in["isWenu"] == 1){
            *f["selLeptons_pt_0"] = f["Electron_pt"][*in["lepInd1"]];
            *f["selLeptons_eta_0"] = f["Electron_eta"][*in["lepInd1"]];
            *f["selLeptons_phi_0"] = f["Electron_phi"][*in["lepInd1"]];
            *f["selLeptons_mass_0"] = f["Electron_mass"][*in["lepInd1"]];
            *in["selLeptons_pdgId_0"] = in["Electron_charge"][*in["lepInd1"]]*11;
        } else if (*in["isWmunu"] == 1){
            *f["selLeptons_pt_0"] = f["Muon_pt"][*in["lepInd1"]];
            *f["selLeptons_eta_0"] = f["Muon_eta"][*in["lepInd1"]];
            *f["selLeptons_phi_0"] = f["Muon_phi"][*in["lepInd1"]];
            *f["selLeptons_mass_0"] = f["Muon_mass"][*in["lepInd1"]];
            *in["selLeptons_pdgId_0"] = in["Muon_charge"][*in["lepInd1"]]*13;
        } else {
            std::cout<<"This shouldn't happen.  isWenu  isWmunu "<<*in["isWenu"]<<" "<<*in["isWmunu"]<<std::endl;
        }

        //FIXME why do we do this?  Cut-flow?  CP
        if (*in["lepInd1"] == -1) {
            // not Wenu or Wmunu, use preselected lepton
            *in["lepInd1"] = 0;
        }

        if (debug > 1000) {
            std::cout << "cutting on dphi(lep, met)" << std::endl;
            std::cout << *in["lepInd1"] << " "
            //          << f["selLeptons_phi"][*in["lepInd1"]] << " "
            //          << f["selLeptons_eta"][*in["lepInd1"]] << " "
            //          << f["selLeptons_mass"][*in["lepInd1"]]
                      << std::endl;
        }

        *f["lepMetDPhi"] = fabs(EvalDeltaPhi(*f["selLeptons_phi_0"], *f["MET_phi"]));

        if (*in["isWmunu"] && *f["lepMetDPhi"] > *f["muMetDPhiCut"]) {
            *in["controlSample"] = -1;
        } else if (*in["isWenu"] && *f["lepMetDPhi"] > *f["elMetDPhiCut"]) {
            *in["controlSample"] = -1;
        }

        if (doOnlySignalRegion && *in["controlSample"] < 0) {
            return false;
        }

        if (sel && *in["controlSample"] > -1) {
            *in["cutFlow"] += 1;
        }

        TLorentzVector Lep, MET;
        // Reconstruct W
        if (debug > 1000) {
            std::cout << "met " << *f["MET_pt"] << " " << *f["MET_phi"] << std::endl;
        }

        MET.SetPtEtaPhiM(*f["MET_pt"], 0., *f["MET_phi"], 0.); // Eta/M don't affect calculation of W.pt and W.phi
        Lep.SetPtEtaPhiM(*f["selLeptons_pt_0"], *f["selLeptons_eta_0"], *f["selLeptons_phi_0"], *f["selLeptons_mass_0"]);
        double cosPhi12 = (Lep.Px()*MET.Px() + Lep.Py()*MET.Py()) / (Lep.Pt() * MET.Pt()); // cos of the angle between the lepton and the missing energy
        *f["V_mt"] = TMath::Sqrt(2*Lep.Pt()*MET.Pt() * (1 - cosPhi12));
        *f["Lep_HJ1_dPhi"] = Lep.DeltaPhi(HJ1);
        *f["Lep_HJ2_dPhi"] = Lep.DeltaPhi(HJ2);
        V = MET + Lep;

        TLorentzVector neutrino = getNu4Momentum(Lep, MET);
        W_withNuFromMWCon = neutrino + Lep;

        //FIXME we should reduce this code to keep only what we use

        *f["Top1_mass_fromLepton_regPT"] = GetRecoTopMass(Lep, false, 0, true); // construct top mass from closest jet to lepton
        *f["Top1_mass_fromLepton_regPT_wMET"] = GetRecoTopMass(Lep, false, 1, true); // construct top mass from closest jet to lepton
        *f["Top1_mass_fromLepton_regPT_w4MET"] = GetRecoTopMass(Lep, false, 2, true); // construct top mass from closest jet to lepton
    } else if (*in["isZnn"] == 1) {
        TLorentzVector MET;
        MET.SetPtEtaPhiM(*f["MET_pt"], 0., *f["MET_phi"], 0.); // Eta/M don't affect calculation of V.pt and V.phi
        V = MET;
    } else {
        std::cerr << "not any selection... can I go now?" << std::endl;  // haha, this is a good one.
    }


//      _|      _|      _|      _|    _|      _|        _|                                                  _|      _|
//      _|      _|      _|      _|    _|      _|  _|        _|_|_|      _|_|    _|_|_|  _|_|      _|_|_|  _|_|_|_|        _|_|_|    _|_|_|
//      _|      _|  _|_|_|_|_|  _|_|_|_|      _|_|      _|  _|    _|  _|_|_|_|  _|    _|    _|  _|    _|    _|      _|  _|        _|_|
//        _|  _|        _|      _|    _|      _|  _|    _|  _|    _|  _|        _|    _|    _|  _|    _|    _|      _|  _|            _|_|
//          _|          _|      _|    _|      _|    _|  _|  _|    _|    _|_|_|  _|    _|    _|    _|_|_|      _|_|  _|    _|_|_|  _|_|_|


    // Calculate V, H and, V+H kinematics

    *f["V_pt"] = V.Pt();

    if (*f["V_pt"] < *f["vptcut"]) {
        *in["controlSample"] = -1;
    } else if (sel && *in["controlSample"] > -1) {
        *in["cutFlow"] += 1; // pT(W) cut
    }

    if (*in["isZee"] == 1 || *in["isZmm"] == 1) {
        *f["V_mass"] = V.M();
        if (*f["V_mass"] < *f["zmasslow"] || *f["V_mass"] > *f["zmasshigh"]) {
            *in["controlSample"] = -1;
        }
    }

    if (doOnlySignalRegion && *in["controlSample"] < 0) {
        return false;
    }

    // di-jet kinematics
    *f["HJ1_HJ2_dPhi"] = HJ1.DeltaPhi(HJ2);
    *f["HJ1_HJ2_dEta"] = fabs(HJ1.Eta() - HJ2.Eta());
    *f["HJ1_HJ2_dR"] = HJ1.DeltaR(HJ2);
    *f["JJEtaBal"] = (fabs(f["Jet_eta"][*in["hJetInd1"]] + f["Jet_eta"][*in["hJetInd2"]])) / (fabs(f["Jet_eta"][*in["hJetInd1"]] - f["Jet_eta"][*in["hJetInd2"]]));
    *f["H_mass_step2"] = *f["H_mass"];
    *f["H_mass_noreg"] = Hbb_noreg.M();
    *f["H_mass"] = Hbb.M(); // mass window cut? regression applied in FinishEvent
    if (cursyst->name != "nominal") {
        *f[Form("H_mass_%s", cursyst->name.c_str())] = Hbb.M();
        *f[Form("H_pt_%s", cursyst->name.c_str())] = Hbb.Pt();
        *in[Form("nAddJets252p9_puid_%s", cursyst->name.c_str())] = *in["nAddJets252p9_puid"];
    }

    if (cursyst->name != "nominal") {
        for (int i = 0; i < *in["nJet"]; i++) {
            f[Form("Jet_btagCMVA_%s", cursyst->name.c_str())][i] = f["Jet_btagCMVA"][i];
        }
    }

    // Now we can calculate whatever we want (transverse) with V and H four-vectors
    *f["jjVPtRatio"] = *f["H_pt"] / *f["V_pt"];
    *f["HVdPhi"] = fabs(Hbb.DeltaPhi(V));
    *f["HVdEta"] = fabs(Hbb.Eta() - V.Eta());
    if (*in["isWmunu"] == 1 || *in["isWenu"] == 1) {
        *f["HVdEta_4MET"] = fabs(Hbb.Eta() -  W_withNuFromMWCon.Eta());
    }


//                                                                    _|
//        _|_|_|    _|_|    _|_|_|      _|_|    _|  _|_|    _|_|_|  _|_|_|_|    _|_|    _|  _|_|
//      _|    _|  _|_|_|_|  _|    _|  _|_|_|_|  _|_|      _|    _|    _|      _|    _|  _|_|
//      _|    _|  _|        _|    _|  _|        _|        _|    _|    _|      _|    _|  _|
//        _|_|_|    _|_|_|  _|    _|    _|_|_|  _|          _|_|_|      _|_|    _|_|    _|
//            _|
//        _|_|

    // FIXME  for the moment pretty much all the gen stuff is broken
    // it relies on branches that don't exist anymore
    //
    ////// Compare gen kinematics for b jets for signal vs. ttbar
    ////if (*in["sampleIndex"] != 0) {
    ////    if (*in["nGenBQuarkFromH"] > 1) {
    ////        // signal event
    ////        GenBJ1.SetPtEtaPhiM(f["GenBQuarkFromH_pt"][0], f["GenBQuarkFromH_eta"][0], f["GenBQuarkFromH_phi"][0], f["GenBQuarkFromH_mass"][0]);
    ////        GenBJ2.SetPtEtaPhiM(f["GenBQuarkFromH_pt"][1], f["GenBQuarkFromH_eta"][1], f["GenBQuarkFromH_phi"][1], f["GenBQuarkFromH_mass"][1]);
    ////    } else if (*in["nGenBQuarkFromTop"] > 0) {
    ////        GenBJ1.SetPtEtaPhiM(f["GenBQuarkFromTop_pt"][0], f["GenBQuarkFromTop_eta"][0], f["GenBQuarkFromTop_phi"][0], f["GenBQuarkFromTop_mass"][0]);
    ////        if (*in["nGenBQuarkFromTop"] > 1) {
    ////            GenBJ2.SetPtEtaPhiM(f["GenBQuarkFromTop_pt"][1], f["GenBQuarkFromTop_eta"][1], f["GenBQuarkFromTop_phi"][1], f["GenBQuarkFromTop_mass"][1]);
    ////        }
    ////    }

    ////    *f["GenBJ1_pt"] = GenBJ1.Pt();
    ////    *f["GenBJ1_eta"] = GenBJ1.Eta();
    ////    *f["GenBJ1_phi"] = GenBJ1.Phi();
    ////    *f["GenBJ1_mass"] = GenBJ1.M();
    ////    *f["GenBJ2_pt"] = GenBJ2.Pt();
    ////    *f["GenBJ2_eta"] = GenBJ2.Eta();
    ////    *f["GenBJ2_phi"] = GenBJ2.Phi();
    ////    *f["GenBJ2_mass"] = GenBJ2.M();

    ////    GenBJJ = GenBJ1 + GenBJ2;
    ////    *f["GenBJJ_pt"] = GenBJJ.Pt();
    ////    *f["GenBJJ_eta"] = GenBJJ.Eta();
    ////    *f["GenBJJ_phi"] = GenBJJ.Phi();
    ////    *f["GenBJJ_mass"] = GenBJJ.M();
    ////    *f["GenBJJ_dPhi"] = GenBJ2.DeltaPhi(GenBJ1);
    ////    *f["GenBJJ_dR"] = GenBJ2.DeltaR(GenBJ1);
    ////    *f["GenBJJ_dEta"] = fabs(GenBJ1.Eta() - GenBJ2.Eta());

    ////    TLorentzVector GenLep1, GenLep2; // closest gen lep to either jet1 or jet2. Sometimes these could be the same lepton.
    ////    double minDR1 = 999;
    ////    double minDR2 = 999;
    ////    double minSecDR2 = 999; //second closest gen lepton to second jet
    ////    int GenLepIndex1 = -1; // index of the lepton closest to jet 1
    ////    int GenLepIndex2 = -1; // index of the lepton closest to jet 2
    ////    int GenLepSecIndex2 = -1; // index of the lepton second closest to jet 2
    ////    for (int i = 0; i < *in["nGenLep"]; i++) {
    ////       TLorentzVector gl;
    ////       gl.SetPtEtaPhiM(f["GenLep_pt"][i], f["GenLep_eta"][i], f["GenLep_phi"][i], f["GenLep_mass"][i] );
    ////       double DR1 = gl.DeltaR(GenBJ1);
    ////       double DR2 = gl.DeltaR(GenBJ2);

    ////       if (DR1 <= minDR1) {
    ////           minDR1 = DR1;
    ////           GenLepIndex1 = i;
    ////       }

    ////       if (DR2 <= minDR2) {
    ////           minSecDR2 = minDR2;
    ////           GenLepSecIndex2 = GenLepIndex2;
    ////           minDR2 = DR2;
    ////           GenLepIndex2 = i;
    ////       } else if (DR2 < minSecDR2) {
    ////           minSecDR2 = DR2;
    ////           GenLepSecIndex2 = i;
    ////       }
    ////    }

    ////    if (GenLepIndex1 == GenLepIndex2) {
    ////        // don't allow us to use the same lepton for each jet
    ////        GenLepIndex2 = GenLepSecIndex2;
    ////    }

    ////    *in["GenLepIndex1"] = GenLepIndex1;
    ////    *in["GenLepIndex2"] = GenLepIndex2;

    ////    if (GenLepIndex1 != -1) {
    ////        GenLep1.SetPtEtaPhiM(f["GenLep_pt"][GenLepIndex1], f["GenLep_eta"][GenLepIndex1], f["GenLep_phi"][GenLepIndex1], f["GenLep_mass"][GenLepIndex1] );
    ////        *f["GenLep_GenBJ1_dR"] = GenLep1.DeltaR(GenBJ1);
    ////        *f["GenLep_GenBJ1_dEta"] = fabs(GenLep1.Eta() - GenBJ1.Eta());
    ////        *f["GenLep_GenBJ1_dPhi"] = GenLep1.DeltaPhi(GenBJ1);

    ////        // try to reconstruct the top mass, although we've lost the neutrino so it will be shifted left
    ////        TLorentzVector GenTop1 = GenLep1 + GenBJ1;
    ////        *f["GenTop1_mass"] = GenTop1.M();
    ////    } else {
    ////        *f["GenLep_GenBJ1_dR"] = -999;
    ////        *f["GenLep_GenBJ1_dEta"] = -999;
    ////        *f["GenLep_GenBJ1_dPhi"] = -999;
    ////        *f["GenTop1_mass"] = -999;
    ////    }

    ////    if (GenLepIndex2 != -1) {
    ////        GenLep2.SetPtEtaPhiM(f["GenLep_pt"][GenLepIndex2], f["GenLep_eta"][GenLepIndex2], f["GenLep_phi"][GenLepIndex2], f["GenLep_mass"][GenLepIndex2] );
    ////        *f["GenLep_GenBJ2_dR"] = GenLep2.DeltaR(GenBJ2);
    ////        *f["GenLep_GenBJ2_dEta"] = (GenLep2.Eta(), GenBJ2.Eta());
    ////        *f["GenLep_GenBJ2_dPhi"] = GenLep2.DeltaPhi(GenBJ2);

    ////        // try to reconstruct the top mass, although we've lost the neutrino so it will be shifted left
    ////        TLorentzVector GenTop2 = GenLep2 + GenBJ2;
    ////        *f["GenTop2_mass"] = GenTop2.M();
    ////    } else {
    ////        *f["GenLep_GenBJ2_dR"] = -999;
    ////        *f["GenLep_GenBJ2_dEta"] = -999;
    ////        *f["GenLep_GenBJ2_dPhi"] = -999;
    ////        *f["GenTop2_mass"] = -999;
    ////    }

    ////    // construct Gen W
    ////    TLorentzVector GenW1, GenW2;
    ////    if (*in["nGenVbosons"] > 0 && fabs(in["GenVbosons_pdgId"][0]) == 24) {
    ////        GenW1.SetPtEtaPhiM(f["GenVbosons_pt"][0], f["GenVbosons_eta"][0], f["GenVbosons_phi"][0], f["GenVbosons_mass"][0]);
    ////        *f["GenW_GenBJJ_dPhi"] = GenW1.DeltaPhi(GenBJJ);
    ////        *f["GenW_GenBJJ_dEta"] = fabs(GenW1.Eta() - GenBJJ.Eta());
    ////        // grab both W's in ttbar events
    ////        if (*in["nGenVbosons"] > 1 && fabs(in["GenVbosons_pdgId"][1]) == 24) {
    ////            GenW2.SetPtEtaPhiM(f["GenVbosons_pt"][1], f["GenVbosons_eta"][1], f["GenVbosons_phi"][1], f["GenVbosons_mass"][1]);
    ////        }
    ////    } else {
    ////        *f["GenW_GenBJJ_dPhi"] = -999;
    ////        *f["GenW_GenBJJ_dEta"] = -999;
    ////    }

    ////    std::vector<TLorentzVector> genWQuarks; // gen quarks from hadronic gen W decay
    ////    for (int i = 0; i < *in["nGenWZQuark"]; i++) {
    ////        TLorentzVector v;
    ////        v.SetPtEtaPhiM(f["GenWZQuark_pt"][i], f["GenWZQuark_eta"][i], f["GenWZQuark_phi"][i], f["GenWZQuark_mass"][i]);
    ////        genWQuarks.push_back(v);
    ////    }

    ////    //int nSelectedJetsMatched = 0; // count the number (0, 1, 2) of selected jets matched to the real bottom quarks
    ////    // Match Jets with Gen B Jets from Higgs/Tops
    ////    for (int i = 0; i < *in["nJet"]; i++) {
    ////        in["Jet_genJetMatchId"][i] = 0; // 0 if no gen match, 1 for pt-leading b-jet, 2 for pt sub-leading b-jet, 3 if matched to jet from hadronic W decay
    ////        TLorentzVector Jet;
    ////        Jet.SetPtEtaPhiM(f["Jet_bReg"][i], f["Jet_eta"][i], f["Jet_phi"][i], f["Jet_mass"][i] * (f["Jet_bReg"][i] / f["Jet_pt"][i]));

    ////        //double dR1 = Jet.DeltaR(GenHJ1);
    ////        //double dR2 = Jet.DeltaR(GenHJ2);
    ////        double dR1 = 999;
    ////        if (GenBJ1.Pt() > 0) dR1 = Jet.DeltaR(GenBJ1);
    ////        double dR2 = 999;
    ////        if (GenBJ2.Pt() > 0) dR2 = Jet.DeltaR(GenBJ2);

    ////        // try to match the jet to one of the jets from hadronic W decay
    ////        double dR3 = 999;
    ////        for (int j = 0; j < (int) genWQuarks.size(); j++) {
    ////            double Jet_genWQuarkDR = Jet.DeltaR(genWQuarks[j]);
    ////            if (Jet_genWQuarkDR < dR3) {
    ////                dR3 = Jet_genWQuarkDR;
    ////            }
    ////        }

    ////        f["Jet_genWQuarkDR"][i] = dR3;
    ////        if (dR3 < std::min(dR1, dR2) && dR3 < 0.5) {
    ////            in["Jet_genJetMatchId"][i] = 3;
    ////        } else if (dR1 <= dR2 && dR1 < 0.5) {
    ////            in["Jet_genJetMatchId"][i] = 1;
    ////        } else if (dR2 < 0.5) {
    ////            in["Jet_genJetMatchId"][i] = 2;
    ////        }

    ////        f["Jet_genHJetMinDR"][i] = std::min(dR1, dR2);

    ////        if (i == *in["hJetInd1"]) {
    ////            *f["hJet1_matchedMinDR"] = f["Jet_genHJetMinDR"][i];
    ////        } else if (i == *in["hJetInd2"]) {
    ////            *f["hJet2_matchedMinDR"] = f["Jet_genHJetMinDR"][i];
    ////        }
    ////    }
    ////}

    if (debug > 1000) std::cout << "counting additional jets and leptons" << std::endl;


//                      _|        _|  _|    _|      _|                                _|        _|            _|
//        _|_|_|    _|_|_|    _|_|_|      _|_|_|_|        _|_|    _|_|_|      _|_|_|  _|                      _|
//      _|    _|  _|    _|  _|    _|  _|    _|      _|  _|    _|  _|    _|  _|    _|  _|        _|            _|
//      _|    _|  _|    _|  _|    _|  _|    _|      _|  _|    _|  _|    _|  _|    _|  _|        _|            _|
//        _|_|_|    _|_|_|    _|_|_|  _|      _|_|  _|    _|_|    _|    _|    _|_|_|  _|        _|    _|      _|
//                                                                                              _|  _|
//                                                                                            _|


    // count the number of additional leptons and jets, then cut on this number
    int nAddLep = 0;

    // 15 to 30 by 1 GeV, 1.5 to 3 w/ 0.1 in eta
    //std::vector<int> ptCuts = {15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35};
    //std::vector<double> etaCuts = {1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.1,3.2,3.3,3.4,3.5};
    std::vector<int> ptCuts = {25};
    std::vector<double> etaCuts = {2.9};

    for (int i = 0; i < (int) ptCuts.size(); i++) {
        for (int j = 0; j < (int) etaCuts.size(); j++) {
            float maxPt = 0; // max. pT of additional jets
            int nAddJet_tmp = 0;
            std::string eta_cut = Form("%.1f", etaCuts[j]); // convert, say, 2.5 to '2p5'
            std::replace(eta_cut.begin(), eta_cut.end(), '.', 'p');
            for (int k = 0; k < m("nJet"); k++) {
                if (k == m("hJetInd1") || k == m("hJetInd2")) continue;
                if (m("Jet_bReg",k) > ptCuts[i] && fabs(m("Jet_eta",k)) < etaCuts[j] && m("Jet_puId",k) > 0) {
                    nAddJet_tmp++;
                    if (f["Jet_bReg"][k] > maxPt) {
                        maxPt = f["Jet_bReg"][k];
                        *f[Form("AddJets%i%s_puid_leadJet_pt", ptCuts[i], eta_cut.c_str())] = maxPt;
                        *f[Form("AddJets%i%s_puid_leadJet_eta", ptCuts[i], eta_cut.c_str())] = f["Jet_eta"][k];
                        *f[Form("AddJets%i%s_puid_leadJet_phi", ptCuts[i], eta_cut.c_str())] = f["Jet_phi"][k];
                        *f[Form("AddJets%i%s_puid_leadJet_btagCSV", ptCuts[i], eta_cut.c_str())] = f["Jet_btagCSVV2"][k];
                    }
                }
            }
            *in[Form("nAddJets%i%s_puid", ptCuts[i], eta_cut.c_str())] = nAddJet_tmp;
        }
    }

    // count additional leptons (check both collections, which are exclusive)
    for (int i = 0; i < m("nMuon"); i++) {
        if (i == *in["muInd1"]) continue; // don't look at the lepton we've selected from the W
        if (i == *in["muInd2"]) continue; // don't look at the lepton we've selected from the Z
        if (f["Muon_pt"][i] > 15 && fabs(f["Muon_eta"][i]) < 2.5 && f["Muon_pfRelIso04_all"][i] < 0.1) {
            nAddLep++;
        }
    }

    for (int i = 0; i < m("nElectron"); i++) {
        if (i == *in["elInd1"]) continue; // don't look at the lepton we've selected from the W
        if (i == *in["elInd2"]) continue; // don't look at the lepton we've selected from the Z
        if (f["Electron_pt"][i] > 15 && fabs(f["Electron_eta"][i]) < 2.5 && f["Electron_pfRelIso03_all"][i] < 0.1) {
            nAddLep++;
        }
    }


    *in["nAddLeptons"] = nAddLep;
    if ((*in["isZnn"] || *in["isWmunu"] || *in["isWenu"]) && nAddLep >= m("nAddLeptonsCut")) {
        *in["controlSample"] = -1;
    } else if (sel && *in["controlSample"] > -1) {
        *in["cutFlow"] += 1; // additional lepton veto
    }

    if ((*in["isZnn"] || *in["isWmunu"] || *in["isWenu"]) && *in["nAddJets252p9_puid"] >= *f["nAddJetsCut"]) {
        *in["controlSample"] = -1;
    } else if (sel && *in["controlSample"] > -1) {
        *in["cutFlow"] += 1; // additional jet veto
    }

    if (*in["isZnn"] && *f["HVdPhi"] < *f["HVDPhiCut_0lepchan"]) {
        *in["controlSample"] = -1;
    } else if ((*in["isWmunu"] || *in["isWenu"]) && *f["HVdPhi"] < *f["HVDPhiCut_1lepchan"]) {
        *in["controlSample"] = -1;
    } else if ((*in["isZmm"] || *in["isZee"]) && m("HVdPhi") < m("HVDPhiCut_2lepchan")) {
        *in["controlSample"] = -1;
    } else if (sel && *in["controlSample"] > -1) {
        *in["cutFlow"] += 1; // dPhi(jj,W) cut
    }

    if (doOnlySignalRegion && *in["controlSample"] < 0) {
        return false;
    }


    // Iterate over the jet collection and count the number of central jets and
    // the number of jets azimuthally within 0.5 radians of the reconstructed
    // Z boson (MET). The latter is used to reject QCD events.
    int nJetsCloseToMet = 0, nJetsCentral = 0;
    for (int i = 0; i < m("nJet"); i++) {
        if (f["Jet_bReg"][i] > 30 && in["Jet_puId"][i] >= 4) {
            if (fabs(f["Jet_eta"][i]) < 2.4) {
                nJetsCentral += 1;
            }
            float absDeltaPhiJetMet = fabs(EvalDeltaPhi(f["Jet_phi"][i], *f["V_phi"]));
            if (absDeltaPhiJetMet < 0.5) {
                nJetsCloseToMet += 1;
            }
        }
    }

    float absDeltaPhiHiggsJet1Met = fabs(EvalDeltaPhi(f["Jet_phi"][*in["hJetInd1"]], *f["V_phi"]));
    float absDeltaPhiHiggsJet2Met = fabs(EvalDeltaPhi(f["Jet_phi"][*in["hJetInd2"]], *f["V_phi"]));
    float minAbsDeltaPhiHiggsJetsMet = std::min(absDeltaPhiHiggsJet1Met, absDeltaPhiHiggsJet2Met);
    float absDeltaPhiMetTrackMet = fabs(EvalDeltaPhi(*f["V_phi"], m("TkMET_phi")));
    *f["minMetjDPhi"] = minAbsDeltaPhiHiggsJetsMet;
    *f["MetTkMetDPhi"] = absDeltaPhiMetTrackMet;

    if (*in["isZnn"] && fabs(m("minMetjDPhi")) < m("minMetjDPhiCut")) {
        *in["controlSample"] = -1;
    } else if (*in["isZnn"] && fabs(m("MetTkMetDPhi")) > m("minMetjDPhiCut")) {
        *in["controlSample"] = -1;
    }

    // if selected for the signal region, return here
    // (if the Hmass cut is not applied, this steals events from the control regions)
    // if (*in["controlSample"] == 0) return true;



//                                      _|                          _|                                                    _|
//        _|_|_|    _|_|    _|_|_|    _|_|_|_|  _|  _|_|    _|_|    _|        _|_|_|    _|_|_|  _|_|_|  _|_|    _|_|_|    _|    _|_|      _|_|_|
//      _|        _|    _|  _|    _|    _|      _|_|      _|    _|  _|      _|_|      _|    _|  _|    _|    _|  _|    _|  _|  _|_|_|_|  _|_|
//      _|        _|    _|  _|    _|    _|      _|        _|    _|  _|          _|_|  _|    _|  _|    _|    _|  _|    _|  _|  _|            _|_|
//        _|_|_|    _|_|    _|    _|      _|_|  _|          _|_|    _|      _|_|_|      _|_|_|  _|    _|    _|  _|_|_|    _|    _|_|_|  _|_|_|
//                                                                                                              _|
//                                                                                                              _|

    if(debug>1000) std::cout<<"about to apply control region selection"<<std::endl;

    // Control Samples
    float absDeltaPhiHiggsMet = *f["HVdPhi"];   // fabs(EvalDeltaPhi(*f["HCMVAV2_reg_phi"], *f["V_phi"]));
    float higgsJet1CMVA = f["Jet_btagCMVA"][*in["hJetInd1"]];
    float higgsJet2CMVA = f["Jet_btagCMVA"][*in["hJetInd2"]];
    float minCSVA = std::min(higgsJet1CMVA, higgsJet2CMVA);
    float maxCSVA = std::max(higgsJet1CMVA, higgsJet2CMVA);
    float CSVAL = -0.5884, CSVAM = 0.4432, CSVAT = 0.9432;
    float V_mass = *f["V_mass"], H_mass = *f["H_mass"];


    // 0-lepton
    bool base0LepCSSelection = (
        // Vector Boson Cuts
        // (*f["Vtype"] == 2 || *f["Vtype"] == 3 || *f["Vtype"] == 4)
        (*in["isWmunu"] || *in["isWenu"] || *in["isZnn"])
        && *in["cutFlow"] >= 2
        && m("V_pt") > 170
        // Higgs Boson Cuts
        && H_mass < 500
        //&& *f["HCMVAV2_pt"] > 120 CP REMOVED
        // Higgs Jet Cuts
        && f["Jet_bReg"][*in["hJetInd1"]] > 60
        && f["Jet_bReg"][*in["hJetInd2"]] > 35
        && higgsJet2CMVA > CSVAL
        && fabs(f["Jet_eta"][*in["hJetInd1"]]) < 2.4
        && fabs(f["Jet_eta"][*in["hJetInd2"]]) < 2.4
        //&& in["Jet_id"][*in["hJetInd1"]] >= 4
        //&& in["Jet_id"][*in["hJetInd2"]] >= 4
        && in["Jet_puId"][*in["hJetInd1"]] >= 4
        && in["Jet_puId"][*in["hJetInd2"]] >= 4
        // Leading Jet Cuts
        //&& in["Jet_id"][0] >= 4
        && in["Jet_puId"][0] >= 4
        // MET Filters
        && m("Flag_goodVertices")
        && m("Flag_globalTightHalo2016Filter")
        && m("Flag_HBHENoiseFilter")
        && m("Flag_HBHENoiseIsoFilter")
        && m("Flag_EcalDeadCellTriggerPrimitiveFilter")
    );


    if (base0LepCSSelection) {
        // if (*f["Vtype"] == 2 || *f["Vtype"] == 3) {
        if (*in["isWmunu"] || *in["isWenu"]) {
            // FIXME nselLeptons needs to be replaced
            //if (minAbsDeltaPhiHiggsJetsMet < 1.57 && *in["nselLeptons"] >= 1 && higgsJet1CMVA > CSVAM && nJetsCentral >= 4 && absDeltaPhiHiggsMet > 2) {
            if (minAbsDeltaPhiHiggsJetsMet < 1.57 && higgsJet1CMVA > CSVAM && nJetsCentral >= 4 && absDeltaPhiHiggsMet > 2) {
                *in["controlSample"] = 1; // TTbar Control Sample Index
            }
        // } else if (*f["Vtype"] == 4) {
        } else if (*in["isZnn"]) {
            bool vetoHiggsMassWindow = m("H_mass") < 60 || m("H_mass") > 160;
            // FIXME nselLeptons needs to be replaced
            //if (nJetsCloseToMet == 0 && absDeltaPhiMetTrackMet < 0.5 && *in["nselLeptons"] == 0 && absDeltaPhiHiggsMet > 2) {
            if (nJetsCloseToMet == 0 && absDeltaPhiMetTrackMet < 0.5 && absDeltaPhiHiggsMet > 2) {
                if (higgsJet1CMVA < CSVAM && nJetsCentral <= 3) {
                    *in["controlSample"] = 2; // Z+Light Control Sample Index
                } else if (higgsJet1CMVA > CSVAT && nJetsCentral == 2 && vetoHiggsMassWindow) {
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
        *in["cutFlow"] >= 2
        && (*in["isWmunu"] || *in["isWenu"])
        && f["Jet_bReg"][*in["hJetInd1"]] > *f["j1ptCut_1lepchan"]
        && f["Jet_bReg"][*in["hJetInd2"]] > *f["j1ptCut_1lepchan"]
        && *f["MET_pt"] > *f["metcut_1lepchan"]
        && *in["nAddLeptons"] == 0
        && *f["V_pt"] > *f["vptcut"]
        && H_mass < 250
        && *f["H_pt"] > *f["hptcut_1lepchan"]
    );

    float maxCSV = std::max(f["Jet_btagCSVV2"][*in["hJetInd1"]], f["Jet_btagCSVV2"][*in["hJetInd2"]]);

    // htJet30 is not currenty in nanoAOD, need to calculate it
    *f["htJet30"] = 0.;
    for (int i=0; i<m("nJet");i++) {
        if (m("Jet_pt",i)>30. && m("Jet_lepFilter",i) && m("Jet_puId",i)) {
            *f["htJet30"] = *f["htJet30"] + m("Jet_pt",i);
        }
    }
    for (int i=0; i<m("nMuon"); i++) {
        if (m("Muon_pt",i)>5 && m("Muon_pfRelIso04_all",i)<0.4) {
            *f["htJet30"] = *f["htJet30"] + m("Muon_pt",i);
        }
    }
    for (int i=0; i<m("nElectron"); i++) {
        if (m("Electron_pt",i)>5 && m("Electron_pfRelIso03_all",i)<0.4) {
            *f["htJet30"] = *f["htJet30"] + m("Electron_pt",i);
        }
    }

    if (base1LepCSSelection) {
        if (maxCSV > CSVAT){ //ttbar or W+HF
            if (m("nAddJets252p9_puid") > 1.5) { //ttbar
                *in["controlSample"] = 11;
            } else if (m("nAddJets252p9_puid") < 0.5 && m("MET_pt")/sqrt(m("htJet30")) > 2.) { //W+HF // remove mass window so we can use the same ntuple for VV, just be careful that we always avoid overlap with SR
                *in["controlSample"] = 13;
            }
        } else if (maxCSV > CSVAL && maxCSV < CSVAM && m("MET_pt")/sqrt(m("htJet30")) > 2.) { //W+LF
            *in["controlSample"] = 12;
        }

        if (*in["sampleIndex"] == 0) {
            std::cout << "data CS event " << *in["controlSample"]
                      << " maxCSV " << maxCSV
                      << " nAddJets252p9_puid " << *in["nAddJets252p9_puid"]
                      << " MET_pt " << m("MET_pt")
                      << " MET_sumEt " << m("MET_sumEt")
                      << " H_mass " << m("H_mass")
                      << std::endl;
        }
    }

    // 2-lepton
    bool base2LepCSSelection = (    // implementing table 15 of AN2015_168_v12, page 75
        *in["cutFlow"] >= 2         // require Vtype and trigger
        && (*in["isZmm"] || *in["isZee"])
        && *f["V_pt"] > 50
    );

    if (base2LepCSSelection) {
        if (////////////////////////// ttbar
            maxCSVA > CSVAT
            && minCSVA > CSVAL
            && (V_mass > 10 && (V_mass < 75 || V_mass > 120))
        ) {
            *in["controlSample"] = 21;
        } else if (/////////////////// Z + LF
            maxCSVA < CSVAL
            && minCSVA < CSVAL
            && absDeltaPhiHiggsMet > 2.5
            && (V_mass > 75 && V_mass < 105)
            && (H_mass > 90 && H_mass <= 150)
        ) {
            *in["controlSample"] = 22;
        } else if (/////////////////// Z + HF
            maxCSVA > CSVAT
            && minCSVA > CSVAL
            && *f["MET_pt"] < 60
            && absDeltaPhiHiggsMet > 2.5
            && (V_mass > 85 && V_mass <= 97)
            && (H_mass <= 90 || H_mass > 150)
        ) {
            *in["controlSample"] = 23;
        }
    }
    // end of 2-lepton


    if (doCutFlow && *in["cutFlow"] >= *in["doCutFlow"]) {
        // keep all preselected events for cutflow
        return true;
    } else {
        return *in["controlSample"] > -1;
    }
}


void VHbbAnalysis::FinishEvent() {

    // b tag weights
    if (bTagCalibReader) {

        if (debug>101) std::cout<<"evaluating btag scale factors"<<std::endl;

        // TODO uncertainties
        float bTagWeight = 1.;
        for (int i=0; i<m("nJet"); i++)
        {
            if (m("Jet_puId", i) > 0
                && fabs(m("Jet_eta", i))<=m("JetEtaCut")
                && m("Jet_pt", i)>m("JetPtPresel"))
            {
                int hadron_flav = m("Jet_hadronFlavour", i);
                auto flav = (hadron_flav==5) ? BTagEntry::FLAV_B :
                            (hadron_flav==4) ? BTagEntry::FLAV_C :
                            BTagEntry::FLAV_UDSG;
                bTagWeight *= bTagCalibReader->eval_auto_bounds(
                    "central",
                    flav,
                    fabs(m("Jet_eta", i)),
                    m("Jet_pt", i),
                    m("Jet_btagDeepB", i)
                );
            }
        }
        *f["bTagWeight"] = bTagWeight;
    }

    //if (bool(*f["doCutFlow"])) {
    //    ofile->cd();
    //    outputTree->Fill();
    //    return;
    //}
    // General use variables

    if (cursample->sampleNum == 49 || cursample->sampleNum == 491) {
        *f["nProcEvents"] = cursample->CountWeighted->GetBinContent(1);
    }
    else {
        *f["nProcEvents"] = cursample->processedEvents;
    }
    //if (*in["sampleIndex"] != 0) {
    //    for (int i=0; i<*in["nLHE_weights_scale"];i++) {
    //        TH1F *CountWeightedLHEWeightScale = cursample->CountWeightedLHEWeightScale;
    //        f["LHE_weights_scale_normwgt"][i] = *f["nProcEvents"] / CountWeightedLHEWeightScale->GetBinContent(CountWeightedLHEWeightScale->FindBin(i));
    //    }
    //    for (int i=0; i<*in["nLHE_weights_pdf"];i++) {
    //        TH1F *CountWeightedLHEWeightPdf = cursample->CountWeightedLHEWeightPdf;
    //        if (CountWeightedLHEWeightPdf->GetBinContent(CountWeightedLHEWeightPdf->FindBin(i)) != 0) {
    //            f["LHE_weights_pdf_normwgt"][i] = *f["nProcEvents"] / CountWeightedLHEWeightPdf->GetBinContent(CountWeightedLHEWeightPdf->FindBin(i));
    //        }
    //        else {
    //            f["LHE_weights_pdf_normwgt"][i] = 1.0;
    //        }
    //        //cout<<f["LHE_weights_pdf_normwgt"][i]<<" = "<<*f["nProcEvents"]<<" / "<<CountWeightedLHEWeightPdf->GetBinContent(CountWeightedLHEWeightPdf->FindBin(i))<<endl;
    //    }
    //}
    //cout<<f["LHE_weights_pdf_normwgt"][5]<<endl;

    if(*in["sampleIndex"]!=0){
        if (*f["genWeight"] > 0) {
            *f["weight"] = cursample->intWeight;
        } else {
            *f["weight"] = cursample->intWeight * -1;
        }
        if (cursample->sampleNum == 49 || cursample->sampleNum == 491) {
            // special prescription for WJets_BGenFilter sample
            *f["weight"] = *f["weight"]*fabs(*f["genWeight"]);
        }
    }
    else {
        *f["weight"] = 1.0;
    }

    *f["weight_ptQCD"] = 1.0;
    *f["weight_ptEWK"] = 1.0;
    // FIXME nominal must be last
    if(cursyst->name=="nominal"){
        ofile->cd();
        if (debug>10000) std::cout<<"filling output tree"<<std::endl;
        outputTree->Fill();
    }
    return;

    // FIXME ... um... ignoring everything?

    if(*in["sampleIndex"]!=0){
        if (*f["doICHEP"] != 1) {
            *f["weight_PU"] = *f["puWeight"];
            *f["weight_PUUp"] = *f["puWeightUp"] / *f["puWeight"];
            *f["weight_PUDown"] = *f["puWeightDown"] / *f["puWeight"];
        }
        else {
            //*f["weight_PU"] = *f["puWeight"];
            //*f["weight_PU"]=ReWeightMC(int(*f["nTrueInt"])+0.5); // it is now a float (continuous distribution?) so we round to the nearest int
            *f["weight_PU"]=puWeight_ichep(int(*f["nTrueInt"])); // it is now a float (continuous distribution?) so we round to the nearest int
            *f["weight_PUUp"]=(puWeight_ichep_up(int(*f["nTrueInt"]))) / *f["weight_PU"]; // it is now a float (continuous distribution?) so we round to the nearest int
            *f["weight_PUDown"]=(puWeight_ichep_down(int(*f["nTrueInt"]))) / *f["weight_PU"]; // it is now a float (continuous distribution?) so we round to the nearest int
        }
        if (*in["nGenTop"]==0 && *in["nGenVbosons"]>0) {
            // only apply to Z/W+jet samples
            *f["weight_ptQCD"]=ptWeightQCD(*in["nGenVbosons"], *f["LHE_HT"], in["GenVbosons_pdgId"][0]);
            *f["weight_ptEWK"]=ptWeightEWK(*in["nGenVbosons"], f["GenVbosons_pt"][0], *f["Vtype"], in["GenVbosons_pdgId"][0]);
        }
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

    if (*f["do2015"] == 1) {
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

    else if (*f["doICHEP"] == 1) {
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
    if (cursample->sampleNum==22 || cursample->sampleNum==44 || cursample->sampleNum==45 || cursample->sampleNum==46 || cursample->sampleNum==47
        || cursample->sampleNum==48 || cursample->sampleNum==49 || cursample->sampleNum==41 || cursample->sampleNum==42 || cursample->sampleNum==43) {
        if (*f["LHE_HT"]>100 && *f["LHE_HT"]<200) {
            if (*f["LHE_Vpt"] > WBjets_ptVMin && *f["LHE_Vpt"] < WBjets_ptVMax && *in["LHE_Nb"] > 0) {
                if (cursample->sampleNum == 48) {
                    WJetStitchWeight = (1 - weightWBjetsHT100);
                }
                else {
                    WJetStitchWeight = weightWBjetsHT100;
                }
            }
            else if (*f["LHE_Vpt"] > WjetsBgen_ptVMin && *f["LHE_Vpt"] < WjetsBgen_ptVMax && *in["LHE_Nb"] == 0 && *in["nGenStatus2bHad"] > 0) {
                if (cursample->sampleNum == 49) {
                    WJetStitchWeight = (1 - weightWjetsBgenHT100);
                }
                else {
                    WJetStitchWeight = weightWjetsBgenHT100;
                }
            }
        }
        if (*f["LHE_HT"]>200 && *f["LHE_HT"]<400) {
            if (*f["LHE_Vpt"] > WBjets_ptVMin && *f["LHE_Vpt"] < WBjets_ptVMax && *in["LHE_Nb"] > 0) {
                if (cursample->sampleNum == 48) {
                    WJetStitchWeight = (1 - weightWBjetsHT200);
                }
                else {
                    WJetStitchWeight = weightWBjetsHT200;
                }
            }
            else if (*f["LHE_Vpt"] > WjetsBgen_ptVMin && *f["LHE_Vpt"] < WjetsBgen_ptVMax && *in["LHE_Nb"] == 0 && *in["nGenStatus2bHad"] > 0) {
                if (cursample->sampleNum == 49) {
                    WJetStitchWeight = (1 - weightWjetsBgenHT200);
                }
                else {
                    WJetStitchWeight = weightWjetsBgenHT200;
                }
            }
        }
        if (*f["LHE_HT"]>400 && *f["LHE_HT"]<600) {
            if (*f["LHE_Vpt"] > WBjets_ptVMin && *f["LHE_Vpt"] < WBjets_ptVMax && *in["LHE_Nb"] > 0) {
                if (cursample->sampleNum == 48) {
                    WJetStitchWeight = (1 - weightWBjetsHT400);
                }
                else {
                    WJetStitchWeight = weightWBjetsHT400;
                }
            }
            else if (*f["LHE_Vpt"] > WjetsBgen_ptVMin && *f["LHE_Vpt"] < WjetsBgen_ptVMax && *in["LHE_Nb"] == 0 && *in["nGenStatus2bHad"] > 0) {
                if (cursample->sampleNum == 49) {
                    WJetStitchWeight = (1 - weightWjetsBgenHT400);
                }
                else {
                    WJetStitchWeight = weightWjetsBgenHT400;
                }
            }
        }
        if (*f["LHE_HT"]>600 && *f["LHE_HT"]<800) {
            if (*f["LHE_Vpt"]  > WBjets_ptVMin && *f["LHE_Vpt"] < WBjets_ptVMax && *in["LHE_Nb"] > 0) {
                if (cursample->sampleNum == 48) {
                    WJetStitchWeight = (1 - weightWBjetsHT600);
                }
                else {
                    WJetStitchWeight = weightWBjetsHT600;
                }
            }
        else if (*f["LHE_Vpt"] > WjetsBgen_ptVMin && *f["LHE_Vpt"] < WjetsBgen_ptVMax && *in["LHE_Nb"] == 0 && *in["nGenStatus2bHad"] > 0) {
                if (cursample->sampleNum == 49) {
                    WJetStitchWeight = (1 - weightWjetsBgenHT600);
                }
                else {
                    WJetStitchWeight = weightWjetsBgenHT600;
                }
            }
        }
        if (*f["LHE_HT"]>800 && *f["LHE_HT"]<1200) {
            if (*f["LHE_Vpt"] > WBjets_ptVMin && *f["LHE_Vpt"] < WBjets_ptVMax && *in["LHE_Nb"] > 0) {
                if (cursample->sampleNum == 48) {
                    WJetStitchWeight = (1 - weightWBjetsHT800);
                }
                else {
                    WJetStitchWeight = weightWBjetsHT800;
                }
            }
            else if (*f["LHE_Vpt"] > WjetsBgen_ptVMin && *f["LHE_Vpt"] < WjetsBgen_ptVMax && *in["LHE_Nb"] == 0 && *in["nGenStatus2bHad"] > 0) {
                if (cursample->sampleNum == 49) {
                    WJetStitchWeight = (1 - weightWjetsBgenHT800);
                }
                else {
                    WJetStitchWeight = weightWjetsBgenHT800;
                }
            }
        }
        if (*f["LHE_HT"]>1200 && *f["LHE_HT"]<2500) {
            if (*f["LHE_Vpt"] > WBjets_ptVMin && *f["LHE_Vpt"] < WBjets_ptVMax && *in["LHE_Nb"] > 0) {
                if (cursample->sampleNum == 48) {
                    WJetStitchWeight = (1 - weightWBjetsHT1200);
                }
                else {
                    WJetStitchWeight = weightWBjetsHT1200;
                }
            }
            else if (*f["LHE_Vpt"] > WjetsBgen_ptVMin && *f["LHE_Vpt"] < WjetsBgen_ptVMax && *in["LHE_Nb"] == 0 && *in["nGenStatus2bHad"] > 0) {
                if (cursample->sampleNum == 49) {
                    WJetStitchWeight = (1 - weightWjetsBgenHT1200);
                }
                else {
                    WJetStitchWeight = weightWjetsBgenHT1200;
                }
            }
        }
        if (*f["LHE_HT"]>2500) {
            if (*f["LHE_Vpt"] > WBjets_ptVMin && *f["LHE_Vpt"] < WBjets_ptVMax && *in["LHE_Nb"] > 0) {
                if (cursample->sampleNum == 48) {
                    WJetStitchWeight = (1 - weightWBjetsHT2500);
                }
                else {
                    WJetStitchWeight = weightWBjetsHT2500;
                }
            }
            else if (*f["LHE_Vpt"] > WjetsBgen_ptVMin && *f["LHE_Vpt"] < WjetsBgen_ptVMax && *in["LHE_Nb"] == 0 && *in["nGenStatus2bHad"] > 0) {
                if (cursample->sampleNum == 49) {
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

    // we need to just save the bTagWeight since we only want to apply it
    // for the nominal shape
    //if (cursyst->name == "nominal") {
    //    *f["weight"] = *f["weight"] * *f["bTagWeight"];
    //}

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
        *in["sampleIndex"] = *in["sampleIndex"]*100;
        *in["sampleIndex_sel"] = *in["sampleIndex_sel"]*100;
        *in["sampleIndex_GenJetSum"] = *in["sampleIndex_GenJetSum"]*100;
        *in["sampleIndex_GenBJetSum"] = *in["sampleIndex_GenBJetSum"]*100;
        *in["sampleIndex_GenJetSumNB"] = *in["sampleIndex_GenJetSumNB"]*100;

        if (fabs(in["Jet_mcFlavour"][*in["hJetInd1"]]) == 5)  *in["bMCFlavorSumSelected"]=*in["bMCFlavorSumSelected"]+1;
        if (fabs(in["Jet_mcFlavour"][*in["hJetInd2"]]) == 5)  *in["bMCFlavorSumSelected"]=*in["bMCFlavorSumSelected"]+1;

        for(int iJet=0;iJet<*in["nJet"];iJet++){
            if(fabs(in["Jet_mcFlavour"][iJet])==5) *in["bMCFlavorSum"]=*in["bMCFlavorSum"]+1;
        }


        for(int indGJ=0; indGJ<*in["nGenJet"]; indGJ++){
            *in["bGenJetSum"]=*in["bGenJetSum"]+1;
        }

        for(int indGJ=0; indGJ<*in["nGenJet"]; indGJ++){
            *in["bGenJetBSum"]=*in["bGenJetBSum"]+in["GenJet_numBHadrons"][indGJ];
        }

        for(int indGJ=0; indGJ<*in["nGenJet"]; indGJ++){
            if (in["GenJet_numBHadrons"][indGJ] > 0) {
                *in["bGenBJetSum"]=*in["bGenBJetSum"] + 1;
            }
        }

        //if(*in["bMCFlavorSum"]==1){
        //    *in["sampleIndex"] = *in["sampleIndex"] + 1;
        //}else if(*in["bMCFlavorSum"]>1){
        //    *in["sampleIndex"] = *in["sampleIndex"] + 2;
        //}

        if(*in["bMCFlavorSumSelected"]==1){
            *in["sampleIndex_sel"] = *in["sampleIndex_sel"] + 1;
        }else if(*in["bMCFlavorSumSelected"]>1){
            *in["sampleIndex_sel"] = *in["sampleIndex_sel"] + 2;
        }

        if(*in["bGenJetBSum"]==1){
            *in["sampleIndex_GenJetSum"] = *in["sampleIndex_GenJetSum"] + 1;
        }else if(*in["bGenJetBSum"]>1){
            *in["sampleIndex_GenJetSum"] = *in["sampleIndex_GenJetSum"] + 2;
        }

        if(*in["bGenBJetSum"]==1){
            *in["sampleIndex_GenBJetSum"] = *in["sampleIndex_GenBJetSum"] + 1;
        }else if(*in["bGenBJetSum"]>1){
            *in["sampleIndex_GenBJetSum"] = *in["sampleIndex_GenBJetSum"] + 2;
        }

        if(*in["bGenJetSum"]==1){
            *in["sampleIndex_GenJetSumNB"] = *in["sampleIndex_GenJetSumNB"] + 1;
        }else if(*in["bGenJetSum"]>1){
            *in["sampleIndex_GenJetSumNB"] = *in["sampleIndex_GenJetSumNB"] + 2;
        }

        // default now the number of b's counted
        //*in["sampleIndex"]=*in["sampleIndex_GenJetSum"];
        *in["sampleIndex"]=*in["sampleIndex_GenBJetSum"];
    }

    // Split by boost category
    if(*in["isWenu"]) {
        if(*f["V_pt"] >= 100 && *f["V_pt"] < 150) *in["eventClass"] += 2;
        else if(*f["V_pt"] >= 150) *in["eventClass"] += 1;
        else *in["eventClass"] += 3;
    }
    else if(*in["isWmunu"]) {
        if(*f["V_pt"] >= 100 && *f["V_pt"] < 130) *in["eventClass"] += 3;
        else if(*f["V_pt"] >= 130 && *f["V_pt"] < 180) *in["eventClass"] += 2;
        else if(*f["V_pt"] >= 180) *in["eventClass"] += 1;
        else *in["eventClass"] += 4;
    }

    if (*in["isWmunu"]) {
        *f["selLeptons_relIso_0"] = f["selLeptons_relIso04"][*in["lepInd1"]];
    }
    else if (*in["isWenu"]) {
        *f["selLeptons_relIso_0"] = f["selLeptons_relIso03"][*in["lepInd1"]];
    }
    else if (*in["isZmm"]) {
        *f["selLeptons_relIso_0"] = f["selLeptons_relIso04"][*in["lepInd1"]];
        *f["selLeptons_relIso_1"] = f["selLeptons_relIso04"][*in["lepInd2"]];
    }
    else if (*in["isZee"]) {
        *f["selLeptons_relIso_0"] = f["selLeptons_relIso03"][*in["lepInd1"]];
        *f["selLeptons_relIso_1"] = f["selLeptons_relIso03"][*in["lepInd2"]];
    }


    if(*in["sampleIndex"]!=0){
        if (*in["nGenLep"] == 0) {
            // gen lep is originally a tau?
            *f["selLeptons_genLepDR_0"] = -1;
        } else {
            TLorentzVector GenLep, El;
            GenLep.SetPtEtaPhiM(f["GenLep_pt"][0],f["GenLep_eta"][0],f["GenLep_phi"][0],f["GenLep_mass"][0]);
            El.SetPtEtaPhiM(f["selLeptons_pt"][*in["lepInd1"]], f["selLeptons_eta"][*in["lepInd1"]], f["selLeptons_phi"][*in["lepInd1"]], f["selLeptons_mass"][*in["lepInd1"]]);
            *f["selLeptons_genLepDR_0"] = El.DeltaR(GenLep);
        }
    }

    // Reconstruct Higgs and W and recalculate variables ourselves
    if(debug>1000) std::cout<<"Making composite candidates"<<std::endl;
    TLorentzVector MET,Lep,W,HJ1,HJ2,Hbb;
    MET.SetPtEtaPhiM(*f["MET_pt"], 0., *f["MET_phi"], 0.); // Eta/M don't affect calculation of W.pt and W.phi
    Lep.SetPtEtaPhiM(f["selLeptons_pt"][*in["lepInd1"]], f["selLeptons_eta"][*in["lepInd1"]], f["selLeptons_phi"][*in["lepInd1"]], f["selLeptons_mass"][*in["lepInd1"]]);
    W = MET + Lep;

    HJ1.SetPtEtaPhiM(f["Jet_bReg"][*in["hJetInd1"]], f["Jet_eta"][*in["hJetInd1"]], f["Jet_phi"][*in["hJetInd1"]], f["Jet_mass"][*in["hJetInd1"]] * (f["Jet_bReg"][*in["hJetInd1"]] / f["Jet_pt"][*in["hJetInd1"]] ) );
    HJ2.SetPtEtaPhiM(f["Jet_bReg"][*in["hJetInd2"]], f["Jet_eta"][*in["hJetInd2"]], f["Jet_phi"][*in["hJetInd2"]], f["Jet_mass"][*in["hJetInd2"]] * (f["Jet_bReg"][*in["hJetInd2"]] / f["Jet_pt"][*in["hJetInd2"]] ) );
    Hbb = HJ1 + HJ2;

    // We already calculate these in Analyze()
    //*d["H_mass"] = Hbb.M();
    //*d["H_pt"] = Hbb.Pt();
    //*d["V_pt"] = W.Pt();
    //*d["HVdPhi"] = Hbb.DeltaPhi(W);

    // Set variables used by the BDT
    //*f["H_mass_f"] = (float) *f["H_mass"];
    //*f["H_pt_f"] = (float) *f["H_pt"];
    //*f["V_pt_f"] = (float) *f["V_pt"];
    *f["hJets_btagCSV_0"] = (float) f["Jet_btagCSVV2"][*in["hJetInd1"]];
    *f["hJets_btagCSV_1"] = (float) f["Jet_btagCSVV2"][*in["hJetInd2"]];
    //*f["HVdPhi_f"] = (float) *f["HVdPhi"];
    *f["H_dEta"] = fabs(f["Jet_eta"][*in["hJetInd1"]] - f["Jet_eta"][*in["hJetInd2"]]);

    //*f["hJets_mt_0"] = HJ1.Mt();
    //*f["hJets_mt_1"] = HJ2.Mt();
    //*f["H_dR"] = (float) HJ1.DeltaR(HJ2);
    //*f["absDeltaPullAngle"] = 0.; //FIXME what is this in the new ntuples??
    *f["hJets_pt_0"] = (float) f["Jet_bReg"][*in["hJetInd1"]];
    *f["hJets_pt_1"] = (float) f["Jet_bReg"][*in["hJetInd2"]];
    //*f["MET_sumEt_f"] = (float) *f["MET_sumEt"]; // is this the right variable??
    *f["nAddJet_f"] = (float) *in["nAddJets252p9_puid"];
    *f["nAddLep_f"] = (float) *in["nAddLeptons"];
    *f["isWenu_f"] = (float) *in["isWenu"];
    *f["isWmunu_f"] = (float) *in["isWmunu"];
    *f["softActivityVH_njets5_f"] = (float) *in["softActivityVH_njets5"];
    //*f["nPVs_f"] = (float) *in["nPVs"];

    if (BDTisSet) {
        if(debug>5000) {std::cout<<"Evaluating BDT..."<<std::endl; }
        for (unsigned int i=0; i<bdtInfos.size(); i++) {

            BDTInfo tmpBDT = bdtInfos[i];
            if(debug>5000) {
                PrintBDTInfoValues(tmpBDT);
                std::cout<<"BDT evaluates to: "<<tmpBDT.reader->EvaluateMVA(tmpBDT.bdtmethod)<<std::endl;
            }
            std::string bdtname(tmpBDT.bdtname);
            if (cursyst->name != "nominal") {
                bdtname.append("_");
                bdtname.append(cursyst->name);
            }
            if (*f["doVPtReweighting"] > 0) {
                float V_pt_nom = *f["V_pt"];
                *f["V_pt"] = *f["V_pt_VPtCorrUp"];
                *f[bdtname+"_VPtCorrUp"] = tmpBDT.reader->EvaluateMVA(tmpBDT.bdtmethod);
                *f["V_pt"] = *f["V_pt_VPtCorrDown"];
                *f[bdtname+"_VPtCorrDown"] = tmpBDT.reader->EvaluateMVA(tmpBDT.bdtmethod);
                *f["V_pt"] = V_pt_nom;
            }
            *f[bdtname] = tmpBDT.reader->EvaluateMVA(tmpBDT.bdtmethod);
        }
    }

    if(jet1EnergyRegressionIsSet && jet2EnergyRegressionIsSet) {
        if(debug>10000) {
            std::cout<<"Evaluating the Jet Energy Regression..."<<std::endl;
            PrintBDTInfoValues(jet1EnergyRegression);
            PrintBDTInfoValues(jet2EnergyRegression);
        }
        double r1Pt = evaluateRegression(*in["hJetInd1"]);
        double r2Pt = evaluateRegression(*in["hJetInd2"]);

        *f["Jet1_regWeight"] = r1Pt/(*f["hJets_pt_0"]);
        *f["Jet2_regWeight"] = r2Pt/(*f["hJets_pt_1"]);

        TLorentzVector hJ1_reg = TLorentzVector();
        TLorentzVector hJ2_reg = TLorentzVector();

        *f["Jet1_pt_reg"] = r1Pt;
        *f["Jet2_pt_reg"] = r2Pt;
        hJ1_reg.SetPtEtaPhiM(r1Pt, f["Jet_eta"][*in["hJetInd1"]], f["Jet_phi"][*in["hJetInd1"]], f["Jet_mass"][*in["hJetInd1"]] * (r1Pt/f["Jet_pt"][*in["hJetInd1"]]) );
        hJ2_reg.SetPtEtaPhiM(r2Pt, f["Jet_eta"][*in["hJetInd2"]], f["Jet_phi"][*in["hJetInd2"]], f["Jet_mass"][*in["hJetInd2"]] * (r2Pt/f["Jet_pt"][*in["hJetInd2"]]) );

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
    *f["CS_SF"] = 1.0;
    if (*in["sampleIndex"]==2200 || *in["sampleIndex"]==4400 || *in["sampleIndex"]==4500 || *in["sampleIndex"]==4600 || *in["sampleIndex"]==4700 || *in["sampleIndex"]==4800 || *in["sampleIndex"]==4900 || *in["sampleIndex"] == 48100 || *in["sampleIndex"] == 49100 || *in["sampleIndex"]==4100 || *in["sampleIndex"]==4200 || *in["sampleIndex"]==4300) {
        *f["CS_SF"] = *f["SF_Wj0b"];
    }
    /*else if (*in["sampleIndex"]==2201 || *in["sampleIndex"]==4401 || *in["sampleIndex"]==4501 || *in["sampleIndex"]==4601 || *in["sampleIndex"]==4701 || *in["sampleIndex"]==4801 || *in["sampleIndex"]==4901 || *in["sampleIndex"]==2202 || *in["sampleIndex"]==4402 || *in["sampleIndex"]==4502 || *in["sampleIndex"]==4602 || *in["sampleIndex"]==4702 || *in["sampleIndex"]==4802 || *in["sampleIndex"]==4902) {
        *f["CS_SF"] = *f["SF_WHF"];
    }*/
    else if (*in["sampleIndex"]==2201 || *in["sampleIndex"]==4401 || *in["sampleIndex"]==4501 || *in["sampleIndex"]==4601 || *in["sampleIndex"]==4701 || *in["sampleIndex"]==4801 || *in["sampleIndex"]==4901 || *in["sampleIndex"] == 48101 || *in["sampleIndex"] == 49101 || *in["sampleIndex"]==4101 || *in["sampleIndex"]==4201 || *in["sampleIndex"]==4301) {
        *f["CS_SF"] = *f["SF_Wj1b"];
    }
    else if (*in["sampleIndex"]==2202 || *in["sampleIndex"]==4402 || *in["sampleIndex"]==4502 || *in["sampleIndex"]==4602 || *in["sampleIndex"]==4702 || *in["sampleIndex"]==4802 || *in["sampleIndex"]==4902 || *in["sampleIndex"] == 48102 || *in["sampleIndex"] == 49102 || *in["sampleIndex"]==4102 || *in["sampleIndex"]==4202 || *in["sampleIndex"]==4302) {
        *f["CS_SF"] = *f["SF_Wj2b"];
    }
    else if (*in["sampleIndex"]==50 || *in["sampleIndex"]==51 || *in["sampleIndex"]==52 || *in["sampleIndex"]==12 || *in["sampleIndex"]==120 || *in["sampleIndex"]==500) {
        *f["CS_SF"] = *f["SF_TT"];
    }

    if (*in["sampleIndex"]!=0) {
        *f["Lep_SF"] = 1.0;
        if (*in["isWmunu"] == 1) {
            if (*f["do2015"] == 1) {
                // used for 2015 analysis
                *f["Lep_SF"] = f["selLeptons_SF_IsoTight"][*in["lepInd1"]] * f["selLeptons_SF_IdCutTight"][*in["lepInd1"]] * f["selLeptons_SF_HLT_RunD4p3"][*in["lepInd1"]];
            }
            else if (*f["doICHEP"] == 1) {
                // for 2016 analysis
                // old V23 prescription
                //*f["Lep_SF"] = f["SF_MuIDLoose"][*in["lepInd1"]] * f["SF_MuIsoTight"][*in["lepInd1"]] * (0.0673 * f["SF_SMuTrig_Block1"][*in["lepInd1"]] + 0.9327 * f["SF_SMuTrig_Block2"][*in["lepInd1"]] );

                *f["Lep_SF"] = f["selLeptons_SF_IdCutTight"][*in["lepInd1"]] * f["selLeptons_SF_IsoTight"][*in["lepInd1"]] * f["selLeptons_SF_trk_eta"][*in["lepInd1"]] * (0.945*f["selLeptons_SF_HLT_RunD4p3"][*in["lepInd1"]] + 0.055*f["selLeptons_SF_HLT_RunD4p2"][*in["lepInd1"]])  ;

            }
            else {
                //*f["Lep_SF"] = ( (20.1/36.4) * f["SF_MuIDTightBCDEF"][*in["lepInd"]] + (16.3/36.4) * f["SF_MuIDTightGH"][*in["lepInd"]]) * ( (20.1/36.4) * f["SF_MuIsoTightBCDEF"][*in["lepInd"]] + (16.3/36.4) * f["SF_MuIsoTightGH"][*in["lepInd"]] ) ;
                *f["Lep_SF"] = ( (20.1/36.4) * f["SF_MuIDTightBCDEF"][*in["lepInd1"]] + (16.3/36.4) * f["SF_MuIDTightGH"][*in["lepInd1"]]) * ( (20.1/36.4) * f["SF_MuIsoTightBCDEF"][*in["lepInd1"]] + (16.3/36.4) * f["SF_MuIsoTightGH"][*in["lepInd1"]] ) *  ( (20.1/36.4) * f["SF_MuTriggerBCDEF"][*in["lepInd1"]] + (16.3/36.4) * f["SF_MuTriggerGH"][*in["lepInd1"]]) * ( (20.1/36.4) * f["SF_MuTrackerBCDEF"][*in["lepInd1"]] + (16.3/36.4) * f["SF_MuTrackerGH"][*in["lepInd1"]]);
                *f["Lep_SFUp"] = ( (20.1/36.4) * (f["SF_MuIDTightBCDEF"][*in["lepInd1"]] + f["SF_MuIDTightBCDEF_err"][*in["lepInd1"]]) + (16.3/36.4) * (f["SF_MuIDTightGH"][*in["lepInd1"]] + f["SF_MuIDTightGH_err"][*in["lepInd1"]]) ) * ( (20.1/36.4) * (f["SF_MuIsoTightBCDEF"][*in["lepInd1"]] + f["SF_MuIsoTightBCDEF_err"][*in["lepInd1"]] ) + (16.3/36.4) * (f["SF_MuIsoTightGH"][*in["lepInd1"]] + f["SF_MuIsoTightGH_err"][*in["lepInd1"]]) ) *  ( (20.1/36.4) * (f["SF_MuTriggerBCDEF"][*in["lepInd1"]] + f["SF_MuTriggerBCDEF_err"][*in["lepInd1"]] )+ (16.3/36.4) * (f["SF_MuTriggerGH"][*in["lepInd1"]]) + f["SF_MuTriggerGH_err"][*in["lepInd1"]]) * ( (20.1/36.4) * (f["SF_MuTrackerBCDEF"][*in["lepInd1"]] + f["SF_MuTrackerBCDEF_err"][*in["lepInd1"]])+ (16.3/36.4) * (f["SF_MuTrackerGH"][*in["lepInd1"]] + f["SF_MuTrackerGH_err"][*in["lepInd1"]] ));
                *f["Lep_SFUp"] = *f["Lep_SFUp"] / *f["Lep_SF"];
                *f["Lep_SFDown"] = ( (20.1/36.4) * (f["SF_MuIDTightBCDEF"][*in["lepInd1"]] - f["SF_MuIDTightBCDEF_err"][*in["lepInd1"]]) + (16.3/36.4) * (f["SF_MuIDTightGH"][*in["lepInd1"]] - f["SF_MuIDTightGH_err"][*in["lepInd1"]]) ) * ( (20.1/36.4) * (f["SF_MuIsoTightBCDEF"][*in["lepInd1"]] - f["SF_MuIsoTightBCDEF_err"][*in["lepInd1"]] ) + (16.3/36.4) * (f["SF_MuIsoTightGH"][*in["lepInd1"]] - f["SF_MuIsoTightGH_err"][*in["lepInd1"]]) ) *  ( (20.1/36.4) * (f["SF_MuTriggerBCDEF"][*in["lepInd1"]] - f["SF_MuTriggerBCDEF_err"][*in["lepInd1"]] )+ (16.3/36.4) * (f["SF_MuTriggerGH"][*in["lepInd1"]]) - f["SF_MuTriggerGH_err"][*in["lepInd1"]]) * ( (20.1/36.4) * (f["SF_MuTrackerBCDEF"][*in["lepInd1"]] - f["SF_MuTrackerBCDEF_err"][*in["lepInd1"]])+ (16.3/36.4) * (f["SF_MuTrackerGH"][*in["lepInd1"]] - f["SF_MuTrackerGH_err"][*in["lepInd1"]] ));
                *f["Lep_SFDown"] = *f["Lep_SFDown"] / *f["Lep_SF"];
            }
        }
        else if (*in["isWenu"] == 1) {
           if (*f["do2015"] == 1) {
               // used for 2015 analysis
               *f["Lep_SF"] = f["selLeptons_SF_IsoTight"][*in["lepInd1"]] * f["selLeptons_SF_IdMVATight"][*in["lepInd1"]] * f["SF_HLT_Ele23_WPLoose"][*in["lepInd1"]];
           }
           else if (*f["doICHEP"] ==1) {
               // for 2016 analysis
               // old for V23
               //*f["Lep_SF"] =   f["SF_ElIdMVATrigWP80"][*in["lepInd1"]] * f["SF_HLT_Ele23_WPLoose"][*in["lepInd1"]];
               *f["Lep_SF"] = f["SF_ElIdMVATrigWP80"][*in["lepInd1"]] * f["EffHLT_Ele27_WPLoose_Eta2p1"][*in["lepInd1"]] * f["SF_egammaEffi_tracker"][*in["lepInd1"]];
           }
           else {
               *f["Lep_SF"] = f["SF_SingleElTrigger"][*in["lepInd1"]] * f["SF_ElIdIso"][*in["lepInd1"]] *  f["SF_egammaEffi_tracker"][*in["lepInd1"]];
               *f["Lep_SFUp"] = (f["SF_SingleElTrigger"][*in["lepInd1"]] + f["SF_SingleElTrigger_err"][*in["lepInd1"]] )* (f["SF_ElIdIso"][*in["lepInd1"]] + f["SF_ElIdIso_err"][*in["lepInd1"]] ) *  (f["SF_egammaEffi_tracker"][*in["lepInd1"]] + f["SF_egammaEffi_tracker_err"][*in["lepInd1"]] );
                *f["Lep_SFUp"] = *f["Lep_SFUp"] / *f["Lep_SF"];
               *f["Lep_SFDown"] = (f["SF_SingleElTrigger"][*in["lepInd1"]] - f["SF_SingleElTrigger_err"][*in["lepInd1"]] )* (f["SF_ElIdIso"][*in["lepInd1"]] - f["SF_ElIdIso_err"][*in["lepInd1"]] ) *  (f["SF_egammaEffi_tracker"][*in["lepInd1"]] - f["SF_egammaEffi_tracker_err"][*in["lepInd1"]] );
                *f["Lep_SFDown"] = *f["Lep_SFDown"] / *f["Lep_SF"];

               // failing SF to correct for data/MC differences in lepton veto yields
               *f["LepFail_SF"] = 1.0;
               *f["LepFail_SFUp"] = 1.0;
               *f["LepFail_SFDown"] = 1.0;
               if (*in["nselLeptons"] > 1) {
                   for (int i=0; i < *in["nselLeptons"]; i++) {
                       if (i == *in["lepInd1"]) continue;
                       *f["LepFail_SF"] = f["SF_EleVeto"][i];
                       *f["LepFail_SFUp"] = f["SF_EleVeto"][i] + f["SF_EleVeto_err"][i];
                       *f["LepFail_SFUp"] = *f["LepFail_SFUp"] / *f["LepFail_SF"];
                       *f["LepFail_SFDown"] = f["SF_EleVeto"][i] - f["SF_EleVeto_err"][i];
                       *f["LepFail_SFDown"] = *f["LepFail_SFDown"] / *f["LepFail_SF"];
                       break; // apply for highest pT lepton besides selected one
                   }
               }

               //*f["Lep_SF"] = f["SF_ElIdMVATrigWP80"][*in["lepInd1"]] * f["SF_egammaEffi_tracker"][*in["lepInd1"]];
               //*f["Lep_SF"] = 1.0;
           }
        }
        // in V23 for 2016 they changed the names of all the btag weights x.x
        if (*f["do2015"] != 1) {
            if (int(*f["doCMVA"]) == 0) {
                *f["bTagWeight"] = *f["btagWeightCSV"];
                *f["bTagWeight_JESUp"] = *f["btagWeightCSV_up_jes"];
                *f["bTagWeight_JESDown"] = *f["btagWeightCSV_down_jes"];
                *f["bTagWeight_HFUp"] = *f["btagWeightCSV_up_hf"];
                *f["bTagWeight_HFDown"] = *f["btagWeightCSV_down_hf"];
                *f["bTagWeight_LFUp"] = *f["btagWeightCSV_up_lf"];
                *f["bTagWeight_LFDown"] = *f["btagWeightCSV_down_lf"];
                *f["bTagWeight_HFStats1Up"] = *f["btagWeightCSV_up_hfstats1"];
                *f["bTagWeight_HFStats1Down"] = *f["btagWeightCSV_down_hfstats1"];
                *f["bTagWeight_LFStats1Up"] = *f["btagWeightCSV_up_lfstats1"];
                *f["bTagWeight_LFStats1Down"] = *f["btagWeightCSV_down_lfstats1"];
                *f["bTagWeight_HFStats2Up"] = *f["btagWeightCSV_up_hfstats2"];
                *f["bTagWeight_HFStats2Down"] = *f["btagWeightCSV_down_hfstats2"];
                *f["bTagWeight_LFStats2Up"] = *f["btagWeightCSV_up_lfstats2"];
                *f["bTagWeight_LFStats2Down"] = *f["btagWeightCSV_down_lfstats2"];
                *f["bTagWeight_cErr1Up"] = *f["btagWeightCSV_up_cferr1"];
                *f["bTagWeight_cErr1Down"] = *f["btagWeightCSV_down_cferr1"];
                *f["bTagWeight_cErr2Up"] = *f["btagWeightCSV_up_cferr2"];
                *f["bTagWeight_cErr2Down"] = *f["btagWeightCSV_down_cferr2"];
            }
            else {
                *f["bTagWeight"] = *f["btagWeightCMVAV2"];
                *f["bTagWeight_JESUp"] = *f["btagWeightCMVAV2_up_jes"];
                *f["bTagWeight_JESDown"] = *f["btagWeightCMVAV2_down_jes"];
                *f["bTagWeight_HFUp"] = *f["btagWeightCMVAV2_up_hf"];
                *f["bTagWeight_HFDown"] = *f["btagWeightCMVAV2_down_hf"];
                *f["bTagWeight_LFUp"] = *f["btagWeightCMVAV2_up_lf"];
                *f["bTagWeight_LFDown"] = *f["btagWeightCMVAV2_down_lf"];
                *f["bTagWeight_HFStats1Up"] = *f["btagWeightCMVAV2_up_hfstats1"];
                *f["bTagWeight_HFStats1Down"] = *f["btagWeightCMVAV2_down_hfstats1"];
                *f["bTagWeight_LFStats1Up"] = *f["btagWeightCMVAV2_up_lfstats1"];
                *f["bTagWeight_LFStats1Down"] = *f["btagWeightCMVAV2_down_lfstats1"];
                *f["bTagWeight_HFStats2Up"] = *f["btagWeightCMVAV2_up_hfstats2"];
                *f["bTagWeight_HFStats2Down"] = *f["btagWeightCMVAV2_down_hfstats2"];
                *f["bTagWeight_LFStats2Up"] = *f["btagWeightCMVAV2_up_lfstats2"];
                *f["bTagWeight_LFStats2Down"] = *f["btagWeightCMVAV2_down_lfstats2"];
                *f["bTagWeight_cErr1Up"] = *f["btagWeightCMVAV2_up_cferr1"];
                *f["bTagWeight_cErr1Down"] = *f["btagWeightCMVAV2_down_cferr1"];
                *f["bTagWeight_cErr2Up"] = *f["btagWeightCMVAV2_up_cferr2"];
                *f["bTagWeight_cErr2Down"] = *f["btagWeightCMVAV2_down_cferr2"];
            }
        }

        if (*f["doCutFlow"]>0 && *in["cutFlow"]<5) {
            // lepton scale factor calculation will break for some events in the cutflow before lepton selection
            *f["weight"] = *f["weight"] * *f["weight_PU"] * *f["bTagWeight"] * *f["CS_SF"] * *f["weight_ptQCD"] * *f["weight_ptEWK"];
        }
        else if (*f["do2015"]==1 || *f["doICHEP"]==1){
            *f["weight"] = *f["weight"] * *f["weight_PU"] * *f["bTagWeight"] * *f["CS_SF"] * *f["weight_ptQCD"] * *f["weight_ptEWK"] * *f["Lep_SF"];
        }
        else {
            // bTagWeight is broken in V25 Heppy ntuples
            *f["weight"] = *f["weight"] * *f["weight_PU"] * *f["CS_SF"] * *f["weight_ptQCD"] * *f["weight_ptEWK"] * *f["Lep_SF"];
        }

        // Add NLO to LO W+jet re-weighting from Z(ll)H(bb)
        float deta_bb = fabs(f["Jet_eta"][*in["hJetInd1"]] - f["Jet_eta"][*in["hJetInd2"]]);
        int sampleIndex = *in["sampleIndex"];
        float WJetNLOWeight = 1.0;
        float weight_ptQCD = *f["weight_ptQCD"];
        if (sampleIndex<4100 || sampleIndex>4902) { WJetNLOWeight = 1.0; }
        else if (sampleIndex==4100 || sampleIndex==4200 || sampleIndex==4300 || sampleIndex==4400 || sampleIndex==4500 || sampleIndex==4600 || sampleIndex==4700 || sampleIndex==4800 || sampleIndex==4900 || sampleIndex==48100 || sampleIndex==49100) {
            WJetNLOWeight = LOtoNLOWeightBjetSplitEtabb(deta_bb, 0);
            WJetNLOWeight = (WJetNLOWeight/weight_ptQCD)*1.21;
        }
        else if (sampleIndex==4101 || sampleIndex==4201 || sampleIndex==4301 || sampleIndex==4401 || sampleIndex==4501 || sampleIndex==4601 || sampleIndex==4701 || sampleIndex==4801 || sampleIndex==4901 || sampleIndex==48101 || sampleIndex==49101) {
            WJetNLOWeight = LOtoNLOWeightBjetSplitEtabb(deta_bb, 1);
            WJetNLOWeight = (WJetNLOWeight/weight_ptQCD)*1.21;
        }
        else if (sampleIndex==4102 || sampleIndex==4202 || sampleIndex==4302 || sampleIndex==4402 || sampleIndex==4502 || sampleIndex==4602 || sampleIndex==4702 || sampleIndex==4802 || sampleIndex==4902 || sampleIndex==48102 || sampleIndex==49102) {
            WJetNLOWeight = LOtoNLOWeightBjetSplitEtabb(deta_bb, 2);
            WJetNLOWeight = (WJetNLOWeight/weight_ptQCD)*1.21;
        }
        *f["weight"] = *f["weight"] * WJetNLOWeight;
        *f["WJetNLOWeight"] = WJetNLOWeight;

    }

    // FIXME nominal must be last
    if(cursyst->name=="nominal"){
        ofile->cd();
        if (debug>10000) std::cout<<"filling output tree"<<std::endl;
        outputTree->Fill();
    }

    if(debug>100) std::cout<<"Ending FinishEvent()"<<std::endl;
    return;
}

void VHbbAnalysis::TermAnalysis(){
    if(debug>10) std::cout<<"START TermAnalysis()"<<std::endl;
    ofile->cd();
    outputTree->Write();
    ofile->Close();
    if(debug>10) std::cout<<"DONE TermAnalysis()"<<std::endl;
    return;
}

std::pair<int,int> VHbbAnalysis::HighestPtGoodElectronsOppCharge(float min_pt, float max_rel_iso, float idcut) {
    int first = -1;
    for (int i = 0; i<*ui["nElectron"]; i++) {
        if (fabs(f["Electron_eta"][i]) < *f["eletacut"]
            && f["Electron_pt"][i] > min_pt
            && f["Electron_pfRelIso03_all"][i]< max_rel_iso
            //&& in["selLeptons_eleMVAIdSppring16GenPurp"][i] >= idcut
            && b["Electron_mvaSpring16GP_WP90"][i] > 0)
        {
            // leptons are pt sorted
            if (first == -1) {
                first = i;
            } else if (in["Electron_charge"][i]*in["Electron_charge"][first] < 0) {
                return std::make_pair(first, i);
            }
        }
    }
    return std::make_pair(first, -1);
}

std::pair<int,int> VHbbAnalysis::HighestPtGoodMuonsOppCharge(float min_pt, float max_rel_iso) {
    int first = -1;
    for (int i = 0; i<*ui["nMuon"]; i++) {
        if (fabs(f["Muon_eta"][i]) < *f["muetacut"]
            && f["Muon_pt"][i] > min_pt
            && f["Muon_pfRelIso04_all"][i] < max_rel_iso
            //&& in["selLeptons_looseIdPOG"][i] >= *f["muidcut"]
            && b["Muon_mediumId"][i] > 0)
        {
            // leptons are pt sorted
            if (first == -1) {
                first = i;
            } else if (in["Muon_charge"][i]*in["Muon_charge"][first] < 0) {
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
        std::cout<<"*f[\"MET_pt\"] "<<*f["MET_pt"]<<std::endl;
        std::cout<<"*f[selLeptons_eleSieie_0] = "<<*f["selLeptons_eleSieie_0"]<<std::endl;
        std::cout<<"*f[hJets_pt_1] = "<<*f["hJets_pt_1"]<<std::endl;
    }
    // there is only one selected electron for Vtype == 3 which is the electron tag
    *in["elInd1"] = -1;
    *in["elInd2"] = -1;
    if(nLep==1){
        float elMaxPt = 0; // max pt of the electrons we select
        for (int i =0; i<*in["nselLeptons"]; i++) {
            if(fabs(in["selLeptons_pdgId"][i])==11
                && f["selLeptons_pt"][i]      > *f["eptcut"]
                && fabs(f["selLeptons_eta"][i]) < *f["eletacut"]
                && f["selLeptons_relIso03"][i]< *f["erelisocut"]
                //&& in["selLeptons_eleCutIdSpring15_25ns_v1"][i] >= *f["elidcut"]
                //&& in["selLeptons_tightId"][i] >= *f["elidcut"]
                //&& in["selLeptons_eleMVAIdSpring15Trig"][i] >= *f["elidcut"]
                && in["selLeptons_eleMVAIdSppring16GenPurp"][i] >= *f["elidcut"]
                ){
                if (f["selLeptons_pt"][i] > elMaxPt) {
                    elMaxPt = f["selLeptons_pt"][i];
                    *in["elInd1"] = i;
                }
            }
        }
        if (*in["elInd1"] == -1) selectEvent = false;
    } else if(nLep==2){
        float elMaxPt = 0; // max pt of the electrons we select
        float elMax2Pt = 0; // 2nd max pt of the electrons we select
        for (int i =0; i<*in["nselLeptons"]; i++) {
            if(fabs(in["selLeptons_pdgId"][i])==11
                && fabs(f["selLeptons_eta"][i]) < *f["eletacut"]
                && f["selLeptons_relIso03"][i]< *f["erelisocut"]
                && in["selLeptons_eleMVAIdSppring16GenPurp"][i] >= *f["elidcut"]
                ){
                if (f["selLeptons_pt"][i] > elMaxPt &&  f["selLeptons_pt"][i] > *f["e1ptcut"]) {
                    elMaxPt = f["selLeptons_pt"][i];
                    *in["elInd1"] = i;
                }
                else if (*in["elInd1"] > -1 && f["selLeptons_pt"][i] > elMax2Pt && (in["selLeptons_charge"][i]*in["selLeptons_charge"][*in["elInd1"]]) < 0) {
                    elMax2Pt = f["selLeptons_pt"][i];
                    *in["elInd2"] = i;
                }
            }
        }
        if (debug>1000) {
            std::cout<<"elInd1/2 = "<<*in["elInd1"]<<", "<<*in["elInd2"]<<std::endl;
        }
        if (*in["elInd1"] == -1 || *in["elInd2"] == -1) selectEvent = false;
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
        std::cout<<"*d[\"MET_pt\"] "<<*f["MET_pt"]<<std::endl;
    }

    // there is only one selected electron for Vtype == 3 which is the electron tag
    // FIXME add configurable cuts

    if(nLep==1){
        *in["muInd1"] = -1;
        float muMaxPt = 0; // max pt of the muons we select
        for (int i=0; i<*in["nselLeptons"]; i++) {
            if(fabs(in["selLeptons_pdgId"][i])==13
                && f["selLeptons_pt"][i]      > *f["muptcut"]
                && fabs(f["selLeptons_eta"][i]) < *f["muetacut"]
                && f["selLeptons_relIso04"][i]< *f["murelisocut"]
                && in["selLeptons_tightId"][i] >= *f["muidcut"]
                ){
                if (f["selLeptons_pt"][i] > muMaxPt) {
                    muMaxPt = f["selLeptons_pt"][i];
                    *in["muInd1"] = i;
                }
            }
        }
        if (*in["muInd1"] == -1) selectEvent = false;
    } else if(nLep==2) {
        *in["muInd1"] = -1;
        *in["muInd2"] = -1;
        float muMaxPt = 0; // max pt of the muons we select
        float muMax2Pt = 0; // 2nd max pt of the muons we select
        for (int i=0; i<*in["nselLeptons"]; i++) {
            if(fabs(in["selLeptons_pdgId"][i])==13
                && fabs(f["selLeptons_eta"][i]) < *f["muetacut"]
                && f["selLeptons_relIso04"][i]< *f["murelisocut"]
                && in["selLeptons_looseIdPOG"][i] >= *f["muidcut"]
                ){
                if (f["selLeptons_pt"][i] > muMaxPt &&  f["selLeptons_pt"][i] > *f["mu1ptcut"]) {
                    muMaxPt = f["selLeptons_pt"][i];
                    *in["muInd1"] = i;
                } else if (*in["muInd1"]>-1 && f["selLeptons_pt"][i] > muMax2Pt
                  && (in["selLeptons_charge"][i]*in["selLeptons_charge"][*in["muInd1"]]) < 0) {
                    muMax2Pt = f["selLeptons_pt"][i];
                    *in["muInd2"] = i;
                }
            }
        }
        if (*in["muInd1"] == -1 || *in["muInd2"] == -1) selectEvent = false;
    } else {
        std::cout<<"nLep is not 1 or 2 "<<nLep<<" - fail it."<<std::endl;
        selectEvent = false;
    }

    return selectEvent;
}

int VHbbAnalysis::UpdatedVType() {

    *f["Vtype_new"] = *in["Vtype"];

    // muon channels are good
    if (*in["Vtype"] == 0 || *in["Vtype"] == 2) {
        return *in["Vtype"];
    }

    // collect good electrons
    std::vector<int> good_electrons;
    for (int i = 0; i < *ui["nElectron"]; i++) {
        if (f["Electron_pt"][i] > 15
            && b["Electron_mvaSpring16GP_WP90"][i] > 0)
        {
            good_electrons.push_back(i);
        }
    }

    // check for Vtype 1
    if (good_electrons.size() >= 2) {
        if (f["Electron_pt"][good_electrons[0]] > 20) {
            int ele0_charge = in["Electron_charge"][good_electrons[0]];
            for (int i=1; i<good_electrons.size(); i++) {
                if (ele0_charge * in["Electron_charge"][good_electrons[i]] < 0) {
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
    if (*in["Vtype"] == 3) {
        return 3;
    }

    // is Vtype 4 but shouldn't be?
    if (*in["Vtype"] == 4 && good_electrons.size()) {
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
    if (*f["dataYear"] == 2015) {
        // for 2015 V21 ntuples
        if (*in["HLT_BIT_HLT_Ele23_WPLoose_Gsf_v"] != 1
            && *in["HLT_BIT_HLT_IsoMu20_v"] != 1
            && *in["HLT_BIT_HLT_IsoTkMu20_v"] != 1
           ) {
            return false;
        }
    } else if (*f["dataYear"] == 2016 && *f["doICHEP"] == 1) {
        // for 2016 V22 ntuples there is no MC HLT simulation we can use, have to just apply the data trigger efficiency
        if (*in["sampleIndex"] == 0) {
            // regular triggers for 2016 data, but none applied to MC
            if (*in["HLT_BIT_HLT_Ele27_eta2p1_WPLoose_Gsf_v"] != 1
                && *in["HLT_BIT_HLT_IsoMu22_v"] != 1
                && *in["HLT_BIT_HLT_IsoTkMu22_v"] != 1
               ) {
                return false;
            }
        }
    } else if (m("dataYear") == 2016) {
        // 0-lepton
        if (vtype == 4
            && *b["HLT_PFMET110_PFMHT110_IDTight"] != 1
            && *b["HLT_PFMET120_PFMHT120_IDTight"] != 1
            //&& *b["HLT_PFMET170_NoiseCleaned"] != 1
            //&& *b["HLT_PFMET170_HBHE_BeamHaloCleaned"] != 1
            //&& *b["HLT_PFMET170_HBHECleaned"] != 1
           ) {
            return false;
        }
        // 1-lepton
        if ((vtype == 2 || vtype == 3)
            //&& *b["HLT_Ele27_eta2p1_WPTight_Gsf"] != 1
            && *b["HLT_IsoMu24"] != 1
            //&& *b["HLT_IsoTkMu24"] != 1
           ) {
            return false;
        }
        // 2-lepton
        if ((vtype == 0 || vtype == 1)
            && *b["HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL"] != 1
            && *b["HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ"] != 1
            //&& *b["HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL"] != 1
            //&& *b["HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ"] != 1
            && *b["HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ"] != 1
           ) {
            return false;
        }
    } else if (*f["dataYear"] == 2017) {
        // 0-lepton
        if (vtype == 4
            && *b["HLT_PFMET110_PFMHT110_IDTight"] != 1
            && *b["HLT_PFMETTypeOne120_PFMHT120_IDTight"] != 1
           ) {
            return false;
        }
        // 1-lepton
        if ((vtype == 2 || vtype == 3)
            && *b["HLT_Ele32_WPTight_Gsf_L1DoubleEG"] != 1
            && *b["HLT_Ele32_WPTight_Gsf"] != 1
            && *b["HLT_IsoMu24"] != 1
           ) {
            return false;
        }
        // 2-lepton
        if ((vtype == 0 || vtype == 1)
            && *b["HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL"] != 1
            && *b["HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL"] != 1
           ) {
            return false;
        }
    } else {
        std::cout<<"What is this?  Run 1?  Run 3?  2018? "<<std::endl;
    }
    return true;
}

std::pair<int,int> VHbbAnalysis::HighestPtBJets(){
    std::pair<int,int> pair(-1,-1);

    for(int i=0; i<*in["nJet"]; i++){
        if(in["Jet_puId"][i] > 0
            && f["Jet_bReg"][i]>*f["j1ptCut"]
            && f["Jet_btagCSVV2"][i]>*f["j1ptCSV"]&&fabs(f["Jet_eta"][i])<=*f["JetEtaCut"]) {
            if( pair.first == -1 ) {
                pair.first = i;
            } else if(f["Jet_bReg"][pair.first]<f["Jet_pt"][i]){
                pair.first = i;
            }
        }
    }

    for(int i=0; i<*in["nJet"]; i++){
        if(i==pair.first) continue;
        if(in["Jet_puId"][i] > 0
            && f["Jet_bReg"][i]>*f["j2ptCut"]
            && f["Jet_btagCSVV2"][i]>*f["j2ptCSV"]&&fabs(f["Jet_eta"][i])<*f["JetEtaCut"]) {
            if( pair.second == -1 ) {
                pair.second = i;
            } else if(f["Jet_bReg"][pair.second]<f["Jet_pt"][i]){
                pair.second = i;
            }
        }
    }

    return pair;
}


std::pair<int,int> VHbbAnalysis::HighestCSVBJets(float j1ptCut, float j2ptCut){
    std::pair<int,int> pair(-1,-1);

    for(int i=0; i<*in["nJet"]; i++){
        if(in["Jet_puId"][i] > 0
            && f["Jet_bReg"][i]>j1ptCut
            &&fabs(f["Jet_eta"][i])<=*f["JetEtaCut"]) {
            if( pair.first == -1 ) {
                pair.first = i;
            } else if(f["Jet_btagCSVV2"][pair.first]<f["Jet_btagCSVV2"][i]){
                pair.first = i;
            }
        }
    }

    for(int i=0; i<*in["nJet"]; i++){
        if(i==pair.first) continue;
        if(in["Jet_puId"][i] > 0
            && f["Jet_bReg"][i]>j2ptCut
            &&fabs(f["Jet_eta"][i])<*f["JetEtaCut"]) {
            if( pair.second == -1 ) {
                pair.second = i;
            } else if(f["Jet_btagCSVV2"][pair.second]<f["Jet_btagCSVV2"][i]){
                pair.second = i;
            }
        }
    }

    // different pt threshold can set the highest CSV value into pair.second
    if (pair.first > -1 && pair.second > -1 && f["Jet_btagCSVV2"][pair.first] < f["Jet_btagCSVV2"][pair.second]) {
        pair = std::make_pair(pair.second, pair.first);
    }

    return pair;
}

std::pair<int,int> VHbbAnalysis::HighestPtJJBJets(){
    std::pair<int,int> pair(-1,-1);

    // dont implement the btag CSV cuts until we've already selected the highest pT(jj) jets
    double maxPtJJ = 0.;
    for (int i=0; i<*in["nJet"]; i++) {
        if(in["Jet_puId"][i] > 0
            && f["Jet_bReg"][i]>*f["j1ptCut"]
            && fabs(f["Jet_eta"][i])<*f["JetEtaCut"]) {
            TLorentzVector jet1;
            jet1.SetPtEtaPhiM(f["Jet_bReg"][i],f["Jet_eta"][i],f["Jet_phi"][i],f["Jet_mass"][i] * (f["Jet_bReg"][i] / f["Jet_pt"][i] ) );
            for (int j=0; j<*in["nJet"]; j++) {
                if (i == j) continue;
                if(in["Jet_puId"][j] > 0
                    && f["Jet_bReg"][j]>*f["j2ptCut"]
                    && fabs(f["Jet_eta"][j])<*f["JetEtaCut"]) {
                    TLorentzVector jet2;
                    jet2.SetPtEtaPhiM(f["Jet_bReg"][j],f["Jet_eta"][j],f["Jet_phi"][j],f["Jet_mass"][j] * (f["Jet_bReg"][j] / f["Jet_pt"][j] ) );
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
                            if (f["Jet_btagCSVV2"][j] > f["Jet_btagCSVV2"][i]) {
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
    if (f["Jet_btagCSVV2"][pair.first] < *f["j1ptCSV"]) pair.first = -1;
    if (f["Jet_btagCSVV2"][pair.second] < *f["j2ptCSV"]) pair.second = -1;
    return pair;
}

double VHbbAnalysis::GetRecoTopMass(TLorentzVector Obj, bool isJet, int useMET, bool regPT) {

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
    for (int i=0; i<*in["nJet"]; i++) {
        if(regPT){
            thisPT=f["Jet_bReg"][i];
        } else {
            thisPT=f["Jet_pt"][i];
        }
        if (thisPT< 30 || f["Jet_btagCSVV2"][i] < 0.5) continue; // only consider jets with some minimal preselection
        TLorentzVector j;
        j.SetPtEtaPhiM(thisPT, f["Jet_eta"][i], f["Jet_phi"][i], f["Jet_mass"][i] * (f["Jet_bReg"][i] / f["Jet_pt"][i] )  );
        double d1 = j.DeltaR(Obj);
        if (d1 <= minDR) {
            minDR = d1;
            ObjClosestIndex = i;
         }
    }
    if (ObjClosestIndex!=-1) {
        if(regPT){
            thisPT=f["Jet_bReg"][ObjClosestIndex];
            Obj2.SetPtEtaPhiM(thisPT, f["Jet_eta"][ObjClosestIndex], f["Jet_phi"][ObjClosestIndex], f["Jet_mass"][ObjClosestIndex] * (thisPT / f["Jet_pt"][ObjClosestIndex] ) );
        } else {
            thisPT=f["Jet_pt"][ObjClosestIndex];
            Obj2.SetPtEtaPhiM(thisPT, f["Jet_eta"][ObjClosestIndex], f["Jet_phi"][ObjClosestIndex], f["Jet_mass"][ObjClosestIndex]);
        }
    } else return -999;


    if (useMET==0) {
        Top = Obj + Obj2;
    }

    if (useMET==1) {
        // try top = lep + jet + met
        // two-particle Mt = sqrt(Et**2 - pt**2)
        TLorentzVector MET;
        MET.SetPtEtaPhiM(*f["MET_pt"],0.,*f["MET_phi"],0.);
        TLorentzVector Obj_transverse, Obj2_transverse; // can only consider transverse (x-y plane) 4-vector components if using MET
        Obj_transverse.SetPxPyPzE(Obj.Px(),Obj.Py(),0,TMath::Sqrt(TMath::Power(Obj.M(),2) + TMath::Power(Obj.Pt(),2)));
        Obj2_transverse.SetPxPyPzE(Obj2.Px(),Obj2.Py(),0,TMath::Sqrt(TMath::Power(Obj2.M(),2) + TMath::Power(Obj2.Pt(),2)));
        Top = Obj_transverse + Obj2_transverse + MET;
    }else if (useMET==2) {
        TLorentzVector MET;
        MET.SetPtEtaPhiM(*f["MET_pt"],0.,*f["MET_phi"],0.);
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
    int nNuSol(0);

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
        nNuSol = 2;

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
          for( int i =0; i< (int)solutions.size();++i){
          if(solutions[i]<0 ) continue;
          double p_x = (solutions[i]*solutions[i]-mW*mW)/(4*pxlep);
          double p_y = ( mW*mW*pylep + 2*pxlep*pylep*p_x -mW*ptlep*solutions[i])/(2*pxlep*pxlep);
          double Delta2 = (p_x-metpx)*(p_x-metpx)+(p_y-metpy)*(p_y-metpy);

                //std::cout<<"intermediate solution1 met x "<<metpx << " min px " << p_x  <<" met y "<<metpy <<" min py "<< p_y << std::endl;

               //std::cout<<"solution1 met x "<<metpx << " min px " << minPx  <<" met y "<<metpy <<" min py "<< minPy << std::endl;
          }
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

float VHbbAnalysis::getVPtCorrFactor(float V_pt) {
    return (1.15 - 0.00086*V_pt);
}
float VHbbAnalysis::getVPtCorrFactorUp(float V_pt) {
    return (1.24 - 0.0017*V_pt);
}
float VHbbAnalysis::getVPtCorrFactorDown(float V_pt) {
    return (1.16 - 0.000068*V_pt);
}


void VHbbAnalysis::smearJets(float JERScale) {
    for (int i=0; i<*in["nJet"]; i++) {
        float corr_nominal = f["Jet_corr_JER"][i];
        if (corr_nominal == -99) { continue; }
        //std::cout<<"original smearing correction: "<<corr_nominal<<std::endl;
        f["Jet_corr_JER"][i] = 1 + JERScale * (f["Jet_corr_JER"][i]-1);
        f["Jet_corr_JERUp"][i] = 1 + JERScale * (f["Jet_corr_JERUp"][i]-1);
        f["Jet_corr_JERDown"][i] = 1 + JERScale * (f["Jet_corr_JERDown"][i]-1);
        //std::cout<<"smearing jet "<<JERScale<<" times the nominal JER by factor "<<(f["Jet_corr_JER"][i] / corr_nominal)<<std::endl;
        //std::cout<<"Jet_pt = "<<f["Jet_pt"][i]<<std::endl;
        f["Jet_pt"][i] = f["Jet_pt"][i] * ( f["Jet_corr_JER"][i] / corr_nominal);
        // std::cout<<"Re-smeared jet_pt = "<<f["Jet_pt"][i]<<std::endl;

        // re-evaluate jet energy regression on re-smeared jets
        //std::cout<<"Jet_pt_reg was: "<<f["Jet_pt_reg"][i]<<std::endl;
        f["Jet_pt_reg_Heppy"][i] = f["Jet_bReg"][i];
        f["Jet_bReg"][i] = evaluateRegression(i);
        //std::cout<<"now it has been re-evaluated to: "<<f["Jet_pt_reg"][i]<<std::endl;
    }
}

float VHbbAnalysis::evaluateRegression(int i) {
    if (f["Jet_bReg"][i] == -99) { return -99; }
    *f["hJets_pt_0"] = float(f["Jet_pt"][i]);
    *f["hJets_eta_0"] = float(f["Jet_eta"][i]);
    TLorentzVector tmp;
    tmp.SetPtEtaPhiM(f["Jet_pt"][i],f["Jet_eta"][i],f["Jet_phi"][i],f["Jet_mass"][i]);
    *f["hJets_mt_0"] = tmp.Mt();
    //std::cout<<"4-vector Mt() is "<<tmp.Mt()<<std::endl;
    float mt = TMath::Sqrt( TMath::Power(tmp.Et(),2) - TMath::Power(tmp.Pt(),2) );
    //std::cout<<"by-hand Mt is "<<mt<<std::endl;
    //*f["hJets_mt_0"] = mt;
    *f["hJets_leadTrackPt_0"] = float(f["Jet_leadTrackPt"][i]);

    if (f["Jet_leptonPtRel"][i] > 0) {
        *f["hJets_leptonPtRel_0"] = float(f["Jet_leptonPtRel"][i]);
    } else {
        *f["hJets_leptonPtRel_0"] = 0.;
    }

    if (f["Jet_leptonPt"][i] > 0) {
        *f["hJets_leptonPt_0"] = float(f["Jet_leptonPt"][i]);
    } else {
        *f["hJets_leptonPt_0"] =  0.;
    }
    if (f["Jet_leptonDeltaR"][i] > 0) {
        *f["hJets_leptonDeltaR_0"] = float(f["Jet_leptonDeltaR"][i]);
    } else {
        *f["hJets_leptonDeltaR_0"] = 0.;
    }

    *f["hJets_neHEF_0"] = float(f["Jet_neHEF"][i]);
    *f["hJets_neEmEF_0"] = float(f["Jet_neEmEF"][i]);
    *f["hJets_vtxMass_0"] = float(f["Jet_vtxMass"][i]);
    *f["hJets_vtxPt_0"] = float(f["Jet_vtxPt"][i]);

    if (f["Jet_vtx3DVal"][i] > 0) {
        *f["hJets_vtx3dL_0"] = f["Jet_vtx3DVal"][i];
    }else {
        *f["hJets_vtx3dL_0"] = 0.;
    }
    *f["hJets_vtxNtracks_0"] = float(f["Jet_vtxNtracks"][i]);
    //*f["hJets_vtx3deL_0"] = f["Jet_vtx3DSig"][i];
    if (f["Jet_vtx3DSig"][i] > 0) {
        *f["hJets_vtx3deL_0"] = f["Jet_vtx3DVal"][i] / f["Jet_vtx3DSig"][i];
    }else {
        *f["hJets_vtx3deL_0"] = 0;
    }
    //*f["nPVs_f"] = (float) *in["nPVs"];
    //PrintBDTInfoValues(jet1EnergyRegression);
    return jet1EnergyRegression.reader->EvaluateRegression(jet1EnergyRegression.bdtmethod)[0];

}

void VHbbAnalysis::SetupFactorizedJECs(std::string variation) {
    std::string corrsToCalc[] = { "AbsoluteStatUp", "AbsoluteStatDown", "AbsoluteScaleUp", "AbsoluteScaleDown", "AbsoluteFlavMapUp", "AbsoluteFlavMapDown", "AbsoluteMPFBiasUp", "AbsoluteMPFBiasDown", "FragmentationUp", "FragmentationDown", "SinglePionECALUp", "SinglePionECALDown", "SinglePionHCALUp", "SinglePionHCALDown", "FlavorQCDUp", "FlavorQCDDown", "TimePtEtaUp", "TimePtEtaDown", "RelativeJEREC1Up", "RelativeJEREC1Down", "RelativeJEREC2Up", "RelativeJEREC2Down", "RelativeJERHFUp", "RelativeJERHFDown", "RelativePtBBUp", "RelativePtBBDown", "RelativePtEC1Up", "RelativePtEC1Down", "RelativePtEC2Up", "RelativePtEC2Down", "RelativePtHFUp", "RelativePtHFDown", "RelativeBalUp", "RelativeBalDown", "RelativeFSRUp", "RelativeFSRDown", "RelativeStatFSRUp", "RelativeStatFSRDown", "RelativeStatECUp", "RelativeStatECDown", "RelativeStatHFUp", "RelativeStatHFDown", "PileUpDataMCUp", "PileUpDataMCDown", "PileUpPtRefUp", "PileUpPtRefDown", "PileUpPtBBUp", "PileUpPtBBDown", "PileUpPtEC1Up", "PileUpPtEC1Down", "PileUpPtEC2Up", "PileUpPtEC2Down", "PileUpPtHFUp", "PileUpPtHFDown", "PileUpMuZeroUp", "PileUpMuZeroDown", "PileUpEnvelopeUp", "PileUpEnvelopeDown", "SubTotalPileUpUp", "SubTotalPileUpDown", "SubTotalRelativeUp", "SubTotalRelativeDown", "SubTotalPtUp", "SubTotalPtDown", "SubTotalScaleUp", "SubTotalScaleDown", "SubTotalAbsoluteUp", "SubTotalAbsoluteDown", "SubTotalMCUp", "SubTotalMCDown", "TotalUp", "TotalDown", "FlavorZJetUp", "FlavorZJetDown", "FlavorPhotonJetUp", "FlavorPhotonJetDown", "FlavorPureGluonUp", "FlavorPureGluonDown", "FlavorPureQuarkUp", "FlavorPureQuarkDown", "FlavorPureCharmUp", "FlavorPureCharmDown", "FlavorPureBottomUp", "FlavorPureBottomDown", "TimeRunBCDUp", "TimeRunBCDDown", "TimeRunEFUp", "TimeRunEFDown", "TimeRunGUp", "TimeRunGDown", "TimeRunHUp", "TimeRunHDown", "CorrelationGroupMPFInSituUp", "CorrelationGroupMPFInSituDown", "CorrelationGroupIntercalibrationUp", "CorrelationGroupIntercalibrationDown", "CorrelationGroupbJESUp", "CorrelationGroupbJESDown", "CorrelationGroupFlavorUp", "CorrelationGroupFlavorDown", "CorrelationGroupUncorrelatedUp", "CorrelationGroupUncorrelatedDown" };

    bool foundVar = false;
    for (int i=0; i < 102; i++) {
        if (variation == corrsToCalc[i].c_str()) foundVar = true;
    }

    if (!foundVar) return;

    for (int j=0; j < *in["nJet"]; j++) {
        float Jet_corr_var = f["Jet_corr_" + variation][j];
        float corr_nominal = f["Jet_corr"][j];
        float scaleShift = Jet_corr_var / corr_nominal;
        f["Jet_corr_"+variation+"_ratio"][j] = scaleShift;
        float Jet_pt_nom = f["Jet_pt"][j];
        float Jet_pt_reg_nom = evaluateRegression(j); // probably safer to re-evaluate the nominal in case there is any residual bias in re-implementating of reg.
        f["Jet_pt"][j] = f["Jet_pt"][j] * scaleShift;
        float Jet_pt_reg_var = evaluateRegression(j); // re-evaluate regression on top of variation
        f["Jet_pt"][j] = Jet_pt_nom;
        f["Jet_pt_reg_corr"+variation+"_ratio"][j] = Jet_pt_reg_var / Jet_pt_reg_nom;
    }

}
