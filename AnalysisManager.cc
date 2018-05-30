#define AnalysisManager_cxx
#include "AnalysisManager.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <TFile.h>
#include <TMath.h>

AnalysisManager::AnalysisManager(){
    intL=20000; // pb^-1
    settingsTree = new TTree("settings","settings");
    debug = 0;
}

void AnalysisManager::Initialize(std::string filename) {
// used to generate this class and read the Tree.
    InitChain(filename);

    ui.clear();
    in.clear();
    f.clear();
    d.clear();
    b.clear();
    uc.clear();

    branches.clear();
    branchInfos.clear();
    existingBranches.clear();

    debug=0;
    if(outputTreeName==""){
        outputTreeName="condensed_tree";
    }
    safemode=1;

    return;
}


AnalysisManager::~AnalysisManager()
{
    if (!fChain) return;
    delete fChain->GetCurrentFile();

    if(debug>10000) std::cout<<"uints"<<std::endl;
    for(std::map<std::string,unsigned int*>::iterator uiit=ui.begin();
            uiit!=ui.end();  ++uiit){
        if(debug>10000) std::cout<<"I'm deleting "<<uiit->first<<std::endl;
        delete uiit->second;
    }

    if(debug>10000) std::cout<<"ints"<<std::endl;
    for(std::map<std::string,int*>::iterator iit=in.begin();
            iit!=in.end();  ++iit){
        if(debug>10000) std::cout<<"I'm deleting "<<iit->first<<std::endl;
        delete iit->second;
    }

    if(debug>10000) std::cout<<"floats"<<std::endl;
    for(std::map<std::string,float*>::iterator fit=f.begin();
            fit!=f.end();  ++fit){
        if(debug>10000) std::cout<<"I'm deleting "<<fit->first<<std::endl;
        delete fit->second;
    }

    if(debug>10000) std::cout<<"doubles"<<std::endl;
    for(std::map<std::string,double*>::iterator dit=d.begin();
            dit!=d.end();  ++dit){
        if(debug>10000) std::cout<<"I'm deleting "<<dit->first<<std::endl;
        delete dit->second;
    }

    if(debug>10000) std::cout<<"bools"<<std::endl;
    for(std::map<std::string,bool*>::iterator bit=b.begin();
            bit!=b.end();  ++bit){
        if(debug>10000) std::cout<<"I'm deleting "<<bit->first<<std::endl;
        delete bit->second;
    }

    //if(debug>1000) std::cout<<"deleting settingsTree"<<std::endl;
    //delete settingsTree;
    for(std::map<std::string,BDTInfo*>::iterator iterBDT=bdtInfos.begin();
           iterBDT!=bdtInfos.end(); iterBDT++){
        if(debug>10000) std::cout<<"I'm deleting "<<iterBDT->first<<std::endl;
        delete iterBDT->second;;
    }
}


void AnalysisManager::AddSample(SampleContainer sample){
    samples.push_back(sample);
}

void AnalysisManager::AddSystematic(SystematicContainer syst){
    systematics.push_back(syst);
}

void AnalysisManager::AddScaleFactor(SFContainer sf) {
    scaleFactors.push_back(sf);
    std::cout<<"added sf: "<<sf.name<<", now going to test it"<<std::endl;
    float sf_err = 1.0;
    float sfac = sf.getScaleFactor(100.,1.2,sf_err);
    std::cout<<sfac<<std::endl;
}

void AnalysisManager::InitializeBTagSF(const std::string & bTagCalibFile) {

    BTagCalibration calib("taggernamedoesntmatter", bTagCalibFile);
    bTagCalibReader.reset(new BTagCalibrationReader(BTagEntry::OP_RESHAPING, "central",{"up_jes","down_jes","up_lf","down_lf","up_hf","down_hf","up_hfstats1","down_hfstats1","up_hfstats2","down_hfstats2","up_lfstats1","down_lfstats1","up_lfstats2","down_lfstats2","up_cferr1","down_cferr1","up_cferr2","down_cferr2"}));

    bTagCalibReader->load(calib, BTagEntry::FLAV_B,     "iterativeFit");
    bTagCalibReader->load(calib, BTagEntry::FLAV_C,     "iterativeFit");
    bTagCalibReader->load(calib, BTagEntry::FLAV_UDSG,  "iterativeFit");
}

void AnalysisManager::AddBDT(std::string settingName, BDTInfo bdt) {
    bdtInfos[settingName]=new BDTInfo(bdt);
    SetupBDT(bdtInfos[settingName]);
}


void AnalysisManager::PrintBDTInfoValues(BDTInfo* bdt) {
    std::cout<<"Printing information for BDT "<<bdt->bdtname<<"..."<<std::endl;
    for (unsigned int i=0; i < bdt->bdtVars.size(); i++) {
        BDTVariable bdtvar = bdt->bdtVars[i];
        std::cout<<"Input variable: "<<bdtvar.varName.c_str()<<", reference in tree: "<<bdtvar.localVarName.c_str()<<", current value: "<<*f[bdtvar.localVarName]<<", isSpec: "<<bdtvar.isSpec<<std::endl;
    }
}

Int_t AnalysisManager::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}


Long64_t AnalysisManager::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (!fChain->InheritsFrom(TChain::Class()))  return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
   }
   return centry;
}


void AnalysisManager::InitChain(std::string filename)
{
    // The InitChain() function is called when the selector needs to initialize
    // a new tree or chain. Typically here the branch addresses and branch
    // pointers of the tree will be set.
    //fChain = new TChain("tree");
    fChain = new TChain("Events");
    //std::cout<<"opening "<<filename.c_str()<<std::endl;
    //TFile* tf = TFile::Open(filename.c_str());
    //std::cout<<"adding to chain"<<std::endl;
    //fChain->Add(tf);
    fChain->Add(filename.c_str());
    fCurrent = -1;
    fChain->SetMakeClass(1);

    if(debug>100) std::cout<<"Resetting branches"<<std::endl;
    ResetBranches();

    // Set special branches
    //fChain->SetBranchAddress("H", &H, &b_H);
    //fChain->SetBranchAddress("V", &V, &b_V);
    //fChain->SetBranchAddress("METtype1corr", &METtype1corr, &b_METtype1corr);
}


void AnalysisManager::SetupBranch(std::string name, int type, int length, int onlyMC, std::string prov, std::string lengthBranch, int allowMissingBranch, int dropBranch){
    branches[name] = new TBranch;
    if(debug>10000) std::cout<<"new TBranch"<<std::endl;
    branchInfos[name] = new BranchInfo(name,type,length,onlyMC,prov,-999.,lengthBranch, allowMissingBranch, dropBranch);
    if(debug>10000) std::cout<<"new BranchInfo"<<std::endl;

    // Only 0-11 are setup with types for the moment.
    if(type>11 || type<0) {
        std::cout<<"Branch "<<name<<" cannot be set to type "<<type<<std::endl;
        return;
    }
    if(debug>10000) std::cout<<"checking type and setting branch address.  type "<<type<<std::endl;

    if(type==0) {
        ui[name] = new unsigned int();
        fChain->SetBranchAddress(name.c_str(), ui[name], &branches[name]);
    } else if(type==1) {
        in[name] = new int();
        fChain->SetBranchAddress(name.c_str(), in[name], &branches[name]);
    } else if(type==2) {
        f[name] = new float();
        fChain->SetBranchAddress(name.c_str(), f[name], &branches[name]);
    } else if(type==3) {
        d[name] = new double();
        fChain->SetBranchAddress(name.c_str(), d[name], &branches[name]);
    } else if(type==4) {
        b[name] = new bool();
        fChain->SetBranchAddress(name.c_str(), b[name], &branches[name]);
    } else if(type==5) {
        uc[name] = new char();
        fChain->SetBranchAddress(name.c_str(), uc[name], &branches[name]);
    }


    if(type>5 && type<11 && length<0) {
        std::cout<<"Types 6-11 are arrays and need a length greater than 0... not "
            <<length<<" name "<<name.c_str()<<std::endl;
        return;
    }

    if(type==6) {
        ui[name] = new unsigned int[length]();
        fChain->SetBranchAddress(name.c_str(), ui[name], &branches[name]);
    } else if(type==7) {
        in[name] = new int[length]();
        fChain->SetBranchAddress(name.c_str(), in[name], &branches[name]);
    } else if(type==8) {
        f[name] = new float[length]();
        fChain->SetBranchAddress(name.c_str(), f[name], &branches[name]);
    } else if(type==9) {
        d[name] = new double[length]();
        fChain->SetBranchAddress(name.c_str(), d[name], &branches[name]);
    } else if(type==10) {
        b[name] = new bool[length]();
        fChain->SetBranchAddress(name.c_str(), b[name], &branches[name]);
    } else if(type==11) {
        uc[name] = new char[length]();
        fChain->SetBranchAddress(name.c_str(), uc[name], &branches[name]);
    }       

    return;
}

// call after adding all pre-existing branches, but before adding any new branches
void AnalysisManager::ConfigureOutputTree() {
    outputTree = fChain->CloneTree(0);
    if(debug>100000) std::cout<<"getting entries of outputTree"<<std::endl;
    if(debug>100000) std::cout<<"entries "<<outputTree->GetEntries()<<std::endl;
}

void AnalysisManager::SetupNewBranch(std::string name, int type, int length, bool newmem, std::string treetype, float val){
    if(debug>1000) {
        std::cout<<"SetupNewBranch "<<name<<std::endl;
        std::cout<<"treetype name type length val \t"<<treetype<<" "<<name<<" "<<type<<" "<<length<<" "<<val<<std::endl;
    }

    TTree* treeptr;
    if(treetype=="output") { // outputtree
        treeptr=outputTree;
    } else if(treetype=="settings") { // settingstree
        treeptr=settingsTree;
    } else {
        std::cout<<"treetype "<<treetype<<" is unknown.  Not setting up "<<name<<std::endl;
        return;
    }

    if(newmem) {
        if(treetype=="output") { // outputtree
            branchInfos[name] = new BranchInfo(name,type,length,false,"new");
        } else if(treetype=="settings") { // settingstree
            branchInfos[name] = new BranchInfo(name,type,length,false,"settings",val);
        }
    }
    if(debug>1000) std::cout<<"BranchInfo instaniated"<<std::endl;


    if(type>11 || type<0) {
        std::cout<<"New Branch "<<name<<" cannot be set to type "<<type<<std::endl;
        return;
    }

    if(type==0) {
        if(newmem) ui[name] = new unsigned int();
        branches[name] = treeptr->Branch(name.c_str(), ui[name], Form("%s/i",name.c_str()));
    } else if(type==1) {
        if(newmem) in[name] = new int();
        branches[name] = treeptr->Branch(name.c_str(), in[name], Form("%s/I",name.c_str()));
    } else if(type==2) {
        if(newmem) f[name] = new float();
        branches[name] = treeptr->Branch(name.c_str(), f[name], Form("%s/F",name.c_str()));
    } else if(type==3) {
        if(newmem) d[name] = new double();
        branches[name] = treeptr->Branch(name.c_str(), d[name], Form("%s/D",name.c_str()));
    } else if(type==4) {
        if(newmem) b[name] = new bool();
        branches[name] = treeptr->Branch(name.c_str(), b[name], Form("%s/O",name.c_str()));
    } /*else if(type==5) {
        if(newmem) uc[name] = new char;
        branches[name] = treeptr->Branch(name.c_str(), uc[name], Form("%s/C",name.c_str()));
    }*/


    if(type>5 && type<12 && length<0) {
        std::cout<<"Types 6-11 are arrays and need a length greater than 0... not "
            <<length<<" name "<<name.c_str()<<std::endl;
        return;
    }

    if(type==6) {
        if(newmem) ui[name] = new unsigned int[length]();
        branches[name] = treeptr->Branch(Form("%s",name.c_str()), ui[name], Form("%s[%i]/i",name.c_str(),length));
    } else if(type==7) {
        if(newmem) in[name] = new int[length]();
        branches[name] = treeptr->Branch(Form("%s",name.c_str()), in[name], Form("%s[%i]/I",name.c_str(),length));
    } else if(type==8) {
        if(newmem) f[name] = new float[length]();
        branches[name] = treeptr->Branch(Form("%s",name.c_str()), f[name], Form("%s[%i]/F",name.c_str(),length));
    } else if(type==9) {
        if(newmem) d[name] = new double[length]();
        branches[name] = treeptr->Branch(Form("%s",name.c_str()), d[name], Form("%s[%i]/D",name.c_str(),length));
    } else if(type==10) {
        if(newmem) b[name] = new bool[length]();
        branches[name] = treeptr->Branch(Form("%s",name.c_str()), b[name], Form("%s[%i]/O",name.c_str(),length));
    } /*else if(type==11) {
        if(newmem) uc[name] = new char[length]();
        branches[name] = treeptr->Branch(Form("%s",name.c_str()), uc[name], Form("%s[%i]/C",name.c_str(),length));
    }*/

    return;
}

void AnalysisManager::ResetBranches(){

    for(std::map<std::string,BranchInfo*>::iterator biit=branchInfos.begin();
            biit!=branchInfos.end(); ++biit) {
        if(debug>100) std::cout<<"Branch "<<biit->second->name<<" of type, prov "<<biit->second->type<<" "<<biit->second->prov<<std::endl;
        if(biit->second->prov == "existing" || biit->second->prov == "early"){
            std::string name(biit->first);
            if(biit->second->type>11 || biit->second->type<0){
                std::cout<<"Branch "<<name<<" of unknown type "<<biit->second->type<<std::endl;
                continue;
            }
            if(biit->second->type%6==0) {
                fChain->SetBranchAddress(name.c_str(), ui[name], &branches[name]);
            } else if(biit->second->type%6==1) {
                fChain->SetBranchAddress(name.c_str(), in[name], &branches[name]);
            } else if(biit->second->type%6==2) {
                fChain->SetBranchAddress(name.c_str(), f[name], &branches[name]);
            } else if(biit->second->type%6==3) {
                fChain->SetBranchAddress(name.c_str(), d[name], &branches[name]);
            } else if(biit->second->type%6==4) {
                fChain->SetBranchAddress(name.c_str(), b[name], &branches[name]);
            } else if(biit->second->type%6==5) {
                fChain->SetBranchAddress(name.c_str(), uc[name], &branches[name]);
            }
        }
    }
}


void AnalysisManager::PrintBranches(){
    std::cout<<"Branches in branch map"<<std::endl;
    for(std::map<std::string,TBranch*>::iterator ibranch=branches.begin();
            ibranch!=branches.end(); ++ibranch){
        std::cout<<ibranch->first<<" "<<branchInfos[ibranch->first]->prov<<std::endl;
    }
}


void AnalysisManager::SetBranches(){
    // only copy branches specified
    if(debug>10) std::cout<<"SetBranchStatus of existing branches"<<std::endl;
    fChain->SetBranchStatus("*", 0);
    for(std::map<std::string,BranchInfo*>::iterator ibranch=branchInfos.begin();
            ibranch!=branchInfos.end(); ++ibranch){
        if(ibranch->second->prov == "existing" || ibranch->second->prov == "early") {
            if(debug>100) std::cout<<"fChain->SetBranchStatus("<<ibranch->first.c_str()<<", 1);"<<std::endl;
            fChain->SetBranchStatus(ibranch->first.c_str(), 1);
        }
        if(ibranch->second->prov == "new") {

        }
    }
}

void AnalysisManager::SetNewBranches(){
    //
    if(debug>10) std::cout<<"SetupNewBranches"<<std::endl;
    for(std::map<std::string,BranchInfo*>::iterator ibranch=branchInfos.begin();
            ibranch!=branchInfos.end(); ++ibranch){
        if(ibranch->second->prov == "new") {
            if(debug>100) std::cout<<"SetupNewBranch "<<ibranch->first.c_str()<<std::endl;
            SetupNewBranch(ibranch->first, ibranch->second->type, ibranch->second->length, false); //newmem is false
        }
        if(ibranch->second->prov == "settings") {
            if(debug>100) std::cout<<"Setting value "<<ibranch->first.c_str()<<std::endl;
            // Branch is already setup properly
            //SetupNewBranch(ibranch->first, ibranch->second->type, ibranch->second->length, false, "settings", ibranch->second->val); //newmem is false
            // Set the value to val which is read from settings.txt
            *f[ibranch->second->name]=ibranch->second->val;
            if(debug>100) std::cout<<"Setting value done "<<ibranch->first.c_str()<<std::endl;
        }
    }
}

//CheckBranchLengths reads in all "early" and "existing" integer branches (branches which could be a lengthBranch),
//then checks that the maximum length set for "early" and "existing" branches is greater than the value of a specified lengthBranch.
//If not, the program exits cleanly.
void AnalysisManager::CheckBranchLengths(Long64_t entry, bool isData){
    for(std::map<std::string,BranchInfo*>::iterator ibranch=branchInfos.begin();
            ibranch!=branchInfos.end(); ++ibranch){
        if( !(isData && ibranch->second->onlyMC) && (ibranch->second->prov=="early"||ibranch->second->prov=="existing") && ibranch->second->type<2){
            if(BranchExists(ibranch->first)){
                branches[ibranch->first]->GetEntry(entry);
            } else if (!(ibranch->second->allowMissingBranch)){
               std::cout<<"Branch "<<ibranch->first<<" is missing and allowMissingBranch is not set. Exiting..."<<std::endl;
               std::exit(0);
            }
        }
    }
    for(std::map<std::string,BranchInfo*>::iterator ibranch=branchInfos.begin();
            ibranch!=branchInfos.end(); ++ibranch){
        if( !(isData && ibranch->second->onlyMC) && (ibranch->second->prov=="early"||ibranch->second->prov=="existing") && ibranch->second->type>5){
            if((ibranch->second->lengthBranch)!=""){
                if(ibranch->second->length < mInt((ibranch->second->lengthBranch).c_str())){
                    std::cout<<"Branch "<<ibranch->first<< " has max length "<<ibranch->second->length<< " but lengthBranch " <<ibranch->second->lengthBranch <<" is "<<mInt((ibranch->second->lengthBranch).c_str())<<std::endl;
                    std::cout<<"Exiting...."<<std::endl;
                    std::exit(0);
                 }
            }
        }
     }
}


void AnalysisManager::GetEarlyEntries(Long64_t entry, bool isData){
    for(std::map<std::string,BranchInfo*>::iterator ibranch=branchInfos.begin();
            ibranch!=branchInfos.end(); ++ibranch){
        if(ibranch->second->prov == "early" && !(isData && ibranch->second->onlyMC)) {
            if(BranchExists(ibranch->first)){
                if(debug>100000) std::cout<<"Getting entry for early branch "<<ibranch->first<<std::endl;
                branches[ibranch->first]->GetEntry(entry);
                if(debug>100000) std::cout<<"Got entry for early branch "<<ibranch->first<<std::endl;
            } else if (!(ibranch->second->allowMissingBranch)){
               std::cout<<"Branch "<<ibranch->first<<" is missing and allowMissingBranch is not set. Exiting..."<<std::endl;
               std::exit(0);
            }
        }
    }
}

bool AnalysisManager::BranchExists(std::string branchName){
    if( existingBranches.find(branchName) != existingBranches.end()){
        return existingBranches.find(branchName)->second;
    } else {
        bool branchstatus = fChain->GetBranchStatus(branchName.c_str());
        existingBranches[branchName] = branchstatus;
        return branchstatus;
    }
}

// Return a std::vector of the names of the samples
std::vector<std::string> AnalysisManager::ListSampleNames() {
    std::vector<std::string> snlist = std::vector<std::string>();
    for(int i=0; i < (int)samples.size(); i++) {
        snlist.push_back(samples[i].sampleName);
    }
    return snlist;
}


// Process all input samples and all events
void AnalysisManager::Loop(std::string sampleName, std::string filename, std::string ofilename, bool doSkim, float startFrac, float endFrac){
    // Specify sample name if we want to run on only a particular sample, specify
    // filenames if we want to run only on specific files from that sample.
    std::vector<std::string> filenames;
    if (!filename.empty()) {
        std::stringstream ss(filename);
        std::string file;
        while (std::getline(ss, file, ','))
        {
            filenames.push_back(file);
        }
    }
    if(!sampleName.empty()){
        bool sampleFound=false;
        for(int i=0; i<(int)samples.size(); i++) {
            if(sampleName == samples[i].sampleName) {
                ofile = new TFile(ofilename.c_str(), "recreate"); // moved here by Jan
                sampleFound=true;
                SampleContainer* onlySample = new SampleContainer(samples[i]);
                if (filenames.size() > 0) {
                    // keep track of total number of processed events for the sample
                    // so we still get the weight right
                    //int processedEvents = onlySample->processedEvents;
                    std::vector<std::string> sampleFiles = onlySample->files;
                    onlySample->files.clear();
                    //onlySample->sampleChain = new TChain("tree");
                    onlySample->sampleChain = new TChain("Events");
                    for (int j=0; j<(int)filenames.size(); j++) {
                        if (std::find(sampleFiles.begin(), sampleFiles.end(), filenames[j]) != sampleFiles.end() ) {
                            onlySample->AddFile(filenames[j].c_str(),1);
                        }
                        else {
                            std::cout<<"Analysis Manager tried to run on file "<<filenames[j]<<" in sample "<<sampleName<<", but this file is not in the sample's list of files. Skipping..."<<std::endl;
                            std::cout<<"Let's print the full list of files for this sample..."<<std::endl;
                            for (int k=0; k<(int)sampleFiles.size(); k++) {
                                std::cout<<sampleFiles[k]<<std::endl;
                            }
                        }
                    }
                    //onlySample->processedEvents = processedEvents;
                    //ofile = new TFile(Form("%s_%s_%i.root",outputTreeName.c_str(),samples[0].sampleName.c_str(),fNum),"recreate");
                }
                else {
                    //ofile = new TFile(Form("%s_%s.root",outputTreeName.c_str(),samples[0].sampleName.c_str()),"recreate");
                }
                //ofile = new TFile(ofilename.c_str(), "recreate");
                samples.clear();
                samples.push_back(*onlySample);
                break;
            }
        }
        if(!sampleFound) {
            std::cout<<"Analysis Manager tried to loop over the individual sample "<<sampleName<<", but this sample is not in the samples list. Skipping..."<<std::endl;
        }
        if(debug>10) std::cout<<"Starting Loop over individual sample %s"<<sampleName<<std::endl;
    }
    else {
        ofile = new TFile(Form("%s.root",outputTreeName.c_str()),"recreate");
    }

    if(debug>10) std::cout<<"Starting Loop"<<std::endl;
    SetBranches();

    // add new branches
    //ofile->cd();
    tempfile = new TFile("temp.root","recreate");
    tempfile->cd();
    //TTree *newtree = fChain->CloneTree(0);
    // for now we will use a default file to set the structure for the output tree
    outputTree = new TTree();
    outputTree = fChain->CloneTree(0); // need this one to only keep the branches you want

    // it would probably be smarter to check also if SetupBranch() has been called since the last time outputTree was initialized
    //if (!outputTree) {
    //    std::cout<<"Did not call ConfigureOutputTree() before Loop(). Assuming there are no new branches..."<<std::endl;
    //    outputTree = fChain->CloneTree(0);
    //}
    // let's add some of our own branches
    SetNewBranches();
    if(!doSkim){
        SetupSystematicsBranches();
    }
    // FIXME add branches to settings regarding splitting
    //SetupNewBranch("jobNum", 3, -1, true, "settings", jobNum);
    ofile->cd();
    settingsTree->Fill();
    settingsTree->Write();
    delete settingsTree;


    if(debug>10) std::cout<<"Done setting up branches; about to Init"<<std::endl;
    InitAnalysis();


    if(debug>10) std::cout<<"About to loop over samples"<<std::endl;
    // loop through one sample at a time
    for (int i = 0; i < (int)samples.size(); i++) {
        cursample = &samples[i];
        cursample->ComputeWeight(*f["intL"]);

        if(cursample->InputPU != NULL){
            TH1D * puhist_to_save = cursample->InputPU;
            ofile->cd();
            puhist_to_save->Write();
            EvaluatePUReweighting(cursample->InputPU);
            if(globalPUTargetUP!= NULL){
                EvaluatePUReweighting(cursample->InputPU,1);
            }
            if(globalPUTargetDOWN!= NULL){
                EvaluatePUReweighting(cursample->InputPU,-1);
            }
        }

        for(int ifile=0; ifile<(int)(cursample->files.size()); ifile++){
            if(std::find(readFiles.begin(), readFiles.end(), cursample->files[ifile]) != readFiles.end() ) {
                // file has already been processed, skip it
                continue;
            }
            readFiles.push_back(cursample->files[ifile]);
            // set fChain to the TChain for the current sample
            InitChain(cursample->files[ifile]);

            // FIXME should have a sample name, but doesn't right now
            if(debug>0) std::cout<<"About to loop over events in "<<cursample->files[ifile]<<std::endl;
            // loop through the events
            Long64_t nentries = round(fChain->GetEntries()*cursample->procEff);
            Long64_t nentry_begin=0;
            Long64_t nentry_end=nentries;
            if (startFrac!=0){
                nentry_begin=ceil(startFrac*nentries);
            }
            if (endFrac!=1){
                nentry_end=ceil(endFrac*nentries);
            }
            if(debug>1) std::cout<<"looping over "<<nentry_end-nentry_begin<<std::endl;
            Long64_t nbytes = 0, nb = 0;
            int saved=0;
            // FIXME need a loop over systematics
            for (Long64_t jentry=nentry_begin; jentry<nentry_end;jentry++) {
            //for (Long64_t jentry=0; jentry<nentries;jentry++) {
                if((jentry%1000==0 && debug>0) || debug>100000)  std::cout<<"entry saved weighted "<<jentry<<" "<<saved<<" "<<saved*cursample->intWeight<<std::endl;
                //if((jentry%10000==0 && debug>0) || debug>100000)  std::cout<<"entry saved weighted "<<jentry<<" "<<saved<<" "<<saved*cursample->intWeight<<std::endl;
                CheckBranchLengths(jentry, cursample->sampleNum==0);
                bool anyPassing=false;
                for(unsigned iSyst=0; iSyst<systematics.size(); iSyst++){
                    GetEarlyEntries(jentry, cursample->sampleNum==0);
                    cursyst=&(systematics[iSyst]);
                    if (cursample->sampleNum == 0 && cursyst->name != "nominal") continue;
                    //ApplySystematics(true);

                    if(debug>100000) std::cout<<"checking preselection"<<std::endl;
                    bool presel = Preselection();
                    if(presel) anyPassing=true;
                }

                if(!anyPassing) continue;


                if(debug>1000) std::cout<<"passed presel; Loading tree"<<std::endl;
                Long64_t ientry = LoadTree(jentry);
                if (ientry < 0) break;
                //nb = fChain->GetEntry(jentry);
                //nbytes += nb;

                anyPassing=false;
                if(!doSkim){
                    for(unsigned iSyst=0; iSyst<systematics.size(); iSyst++){
                        nb = fChain->GetEntry(jentry);
                        nbytes += nb;

                        cursyst=&(systematics[iSyst]);
                        // set bdt inputs and output to std value before evaluating/saving the event
                        for(std::map<std::string,BDTInfo*>::iterator itBDTInfo=bdtInfos.begin(); itBDTInfo!=bdtInfos.end(); itBDTInfo++){
                            InitializeBDTVariables(itBDTInfo->second);
                            std::string bdt_syst_name = itBDTInfo->second->bdtname;
                            if(cursyst->name != "nominal"){
                                bdt_syst_name.append("_");
                                bdt_syst_name.append(cursyst->name);
                                *f[bdt_syst_name] = -99;
                            }
                        }


                        if (cursample->sampleNum == 0 && cursyst->name != "nominal") continue;
                        ApplySystematics();
                        if(debug>1000) std::cout<<"running analysis"<<std::endl;
                        bool select = Analyze();
                        *b[Form("Pass_%s",cursyst->name.c_str())] = false;
                        if(select) {
                            anyPassing=true;
                            *b[Form("Pass_%s",cursyst->name.c_str())] = true;
                        }
                        if(select || (cursyst->name=="nominal" && anyPassing)){
                            if(debug>1000) std::cout<<"selected event; Finishing"<<std::endl;
                            for (unsigned i=0; i < scaleFactors.size(); i++) {
                                SFContainer sf = scaleFactors[i];
                                float sf_err = 0.0;
                                if (cursample->sampleNum != 0) {
                                    for (int j=0; j < mInt(sf.length); j++) {
                                        if (sf.binning.find("abs") < 0) {  
                                            f[sf.branchname][j] = sf.getScaleFactor(m(sf.branches[0],j) , m(sf.branches[1],j), sf_err);
                                        }else {
                                            f[sf.branchname][j] = sf.getScaleFactor(fabs(m(sf.branches[0],j)) , fabs(m(sf.branches[1],j)), sf_err);
                                        }
                                        f[Form("%s_err",sf.branchname.c_str())][j] = sf_err;
                                    }
                                }
                                else {
                                    // data event, scale factor should just be 1.0
                                    for (int j=0; j < mInt(sf.length); j++) {
                                        f[sf.branchname][j] = 1.0;
                                        f[Form("%s_err",sf.branchname.c_str())][j] = 0.0;
                                    }
                                }
                            
                            }
                            FinishEvent();
                            if(cursyst->name=="nominal") saved++;
                        }
                    }
                 } else {
                     nb = fChain->GetEntry(jentry);
                     nbytes += nb;
                     // running skim
                     //ofile->cd();
                     tempfile->cd();
                     outputTree->Fill();
                     saved++;
                }
            } // end event loop
        } // end file loop
        ofile->cd();
        if (!doSkim) {
            cursample->CountWeightedLHEWeightScale->Write(Form("CountWeightedLHEWeightScale_%s",cursample->sampleName.c_str()));
            cursample->CountWeightedLHEWeightPdf->Write(Form("CountWeightedLHEWeightPdf_%s",cursample->sampleName.c_str()));
            cursample->CountWeighted->Write(Form("CountWeighted_%s",cursample->sampleName.c_str()));
            //cursample->CountFullWeighted->Write(Form("CountFullWeighted_%s",cursample->sampleName.c_str()));
        }
        else {
            cursample->CountWeightedLHEWeightScale->Write();
            cursample->CountWeightedLHEWeightPdf->Write();
            cursample->CountWeighted->Write();
            //cursample->CountFullWeighted->Write();
        }
    } // end sample loop
    std::cout<<"Finished looping"<<std::endl;


    TermAnalysis(doSkim);
}


void AnalysisManager::InitAnalysis(){
    if(debug>100) std::cout<<"InitAnalysis"<<std::endl;
}


bool AnalysisManager::Preselection(){
    bool sel=false;
    if(1) sel=true;
    return sel;
}


bool AnalysisManager::Analyze(){
    bool sel=false;
    return sel;
}


void AnalysisManager::FinishEvent(){
    //need to fill the tree, hist, or whatever here.

    // FIXME nominal must be last!
    if(cursyst->name=="nominal"){
        //ofile->cd();
        tempfile->cd();
        outputTree->Fill();
    }
    return;
}


void AnalysisManager::TermAnalysis(bool skim) {
    // save tree here
    ofile->cd();
    outputTreeSlim = outputTree->CloneTree();
    outputTreeSlim->Write();
    ofile->Close();
}


// Set up all the BDT branches and configure the BDT's with the same input variables as used in training. Run before looping over events.
void AnalysisManager::SetupBDT(BDTInfo* bdtInfo) {

    for (unsigned int i=0; i < bdtInfo->bdtVars.size(); i++) {
        BDTVariable bdtvar = bdtInfo->bdtVars[i];
        std::string thisVar("bdtInput_");
        thisVar.append(bdtInfo->bdtname);
        thisVar.append("_");
        thisVar.append(bdtvar.localVarName);
        if (!bdtvar.isSpec) {
            bdtInfo->reader->AddVariable( bdtvar.varName, f[thisVar]);
        } else {
            bdtInfo->reader->AddSpectator(bdtvar.varName, f[thisVar]);
        }
    }

    std::cout<<"booking MVA for bdt with name...  "<<bdtInfo->bdtname<<std::endl;
    //bdtInfo->reader->BookMVA(bdtInfo->bdtmethod, bdtInfo->xmlFile);
    bdtInfo->BookMVA();
}

void AnalysisManager::InitializeBDTVariables(BDTInfo* bdtInfo){
    for (unsigned int i=0; i < bdtInfo->bdtVars.size(); i++) {
        BDTVariable bdtvar = bdtInfo->bdtVars[i];
        std::string thisVar("bdtInput_");
        thisVar.append(bdtInfo->bdtname);
        thisVar.append("_");
        thisVar.append(bdtvar.localVarName);
        *f[thisVar]=-99;
    }
    *f[bdtInfo->bdtname] = -99;
}

void AnalysisManager::SetBDTVariables(BDTInfo* bdtInfo){
    for (unsigned int i=0; i < bdtInfo->bdtVars.size(); i++) {
        BDTVariable bdtvar = bdtInfo->bdtVars[i];
        std::string thisVar("bdtInput_");
        thisVar.append(bdtInfo->bdtname);
        thisVar.append("_");
        thisVar.append(bdtvar.localVarName);
        *f[thisVar]=*f[bdtvar.localVarName];
    }
}

float AnalysisManager::EvaluateMVA(BDTInfo* bdtInfo){
    SetBDTVariables(bdtInfo);
    return bdtInfo->reader->EvaluateMVA(bdtInfo->methodName);
}

float AnalysisManager::EvaluateRegression(BDTInfo* bdtInfo){
    SetBDTVariables(bdtInfo);
    return bdtInfo->reader->EvaluateRegression(bdtInfo->methodName)[0];
}

void AnalysisManager::SetupSystematicsBranches(){
    std::cout<<"loop through systematics "<<systematics.size()<<std::endl;
    for(unsigned iSyst=0; iSyst<systematics.size(); iSyst++){
        std::cout<<"syst name "<<systematics[iSyst].name<<std::endl;
        for(unsigned iBrnch=0; iBrnch<systematics[iSyst].branchesToEdit.size(); iBrnch++){
            std::string newName(systematics[iSyst].branchesToEdit[iBrnch]);
            newName.append("_");
            newName.append(systematics[iSyst].name);
            std::cout<<"\tbranch name "<<newName<<std::endl;
            std::cout<<"\t\ttype, length,  "
                <<branchInfos[systematics[iSyst].branchesToEdit[iBrnch]]->type<<", "
                <<branchInfos[systematics[iSyst].branchesToEdit[iBrnch]]->length<<std::endl;
                //<<branchInfos[systematics[iSyst].branchesToEdit[iBrnch]]->length<<" "
            SetupNewBranch(newName,
                branchInfos[systematics[iSyst].branchesToEdit[iBrnch]]->type,
                branchInfos[systematics[iSyst].branchesToEdit[iBrnch]]->length);
        }
        for(std::map<std::string,BDTInfo*>::iterator iterBDT=bdtInfos.begin();
              iterBDT!=bdtInfos.end(); iterBDT++){
            std::string bdtname(iterBDT->second->bdtname);
            if (systematics[iSyst].name != "nominal") {
                bdtname.append("_");
                bdtname.append(systematics[iSyst].name);
                SetupNewBranch(bdtname, 2);
            }
        }
        if (systematics[iSyst].name != "nominal") {
            SetupNewBranch(Form("H_mass_%s", systematics[iSyst].name.c_str()), 2);
            SetupNewBranch(Form("H_pt_%s", systematics[iSyst].name.c_str()), 2);
            SetupNewBranch(Form("V_pt_%s", systematics[iSyst].name.c_str()), 2);
            SetupNewBranch(Form("controlSample_%s", systematics[iSyst].name.c_str()), 1);
            SetupNewBranch(Form("weight_%s", systematics[iSyst].name.c_str()), 2);
            SetupNewBranch(Form("Jet_btagCSV_%s", systematics[iSyst].name.c_str()), 8, 250);
            SetupNewBranch(Form("nAddJets252p9_puid_%s", systematics[iSyst].name.c_str()), 1);
        }
    }
}

// applies scales to floats and doubles (and arrays)
// smearing can be added; other types can be added.
void AnalysisManager::ApplySystematics(bool early){

    for(unsigned iBrnch=0; iBrnch<cursyst->branchesToEdit.size(); iBrnch++){
        //std::cout<<"iBrnch "<<iBrnch<<std::endl;
        std::string oldBranchName(cursyst->branchesToEdit[iBrnch]);
        BranchInfo* existingBranchInfo = branchInfos[oldBranchName.c_str()];
        //std::cout<<"early? "<< branchInfos[cursyst->branchesToEdit[iBrnch]]->prov<<std::endl;
        if(!early || (early && branchInfos[cursyst->branchesToEdit[iBrnch]]->prov == "early")){
            std::string systBranchName(oldBranchName);
            systBranchName.append("_");
            systBranchName.append(cursyst->name);
            //std::cout<<"new branch "<<systBranchName.c_str()<<std::endl;
            int thisType=existingBranchInfo->type;
            // just doing floats and doubles
            // smearing can be added at any time
            float jetPtSplit = 100.; // classify high/low jets by this pT threshold
            float jetEtaSplit = 1.4; // classify central/forward jets by this eta threshold
            int ptSplit = 0; // if 0 don't split, if 1 split low, if 2 splight high
            int etaSplit = 0; // if 0 don't split, if 1 split low, if 2 splight high
            if (cursyst->name.find("High")!=std::string::npos) { ptSplit = 2; }
            if (cursyst->name.find("Low")!=std::string::npos) { ptSplit = 1; }
            if (cursyst->name.find("Forward")!=std::string::npos) { etaSplit = 2; }
            if (cursyst->name.find("Central")!=std::string::npos) { etaSplit = 1; }
            //std::cout<<cursyst->name<<std::endl;
            //std::cout<<ptSplit<<", "<<etaSplit<<std::endl;
            if(thisType==2){
                // scale the current branch
                if (cursyst->scaleVar[iBrnch] == "") {
                    // flat scaling
                    *f[systBranchName.c_str()]=*f[oldBranchName.c_str()] * cursyst->scales[iBrnch];
                }
                else {
                    if (cursyst->scaleVarRef[iBrnch] == ""){
                        // dynamic scaling
                        *f[systBranchName.c_str()]=*f[oldBranchName.c_str()] * *f[cursyst->scaleVar[iBrnch]];
                    } else if (cursyst->scaleVarRef[iBrnch] == "None"){
                        *f[systBranchName.c_str()] = *f[cursyst->scaleVar[iBrnch]];
                    } else {
                        *f[systBranchName.c_str()] = *f[oldBranchName.c_str()] * (*f[cursyst->scaleVar[iBrnch]]/ *f[cursyst->scaleVarRef[iBrnch]]);
                    }
                }
                // copy the value to the new branch
                //*f[systBranchName.c_str()]=*f[oldBranchName.c_str()];
            } else if(thisType==3){
                // scale the current branch
                if (cursyst->scaleVar[iBrnch] == "") {
                    // flat scaling
                    *d[systBranchName.c_str()]=*d[oldBranchName.c_str()] * cursyst->scales[iBrnch];
                }
                else {
                    if (cursyst->scaleVarRef[iBrnch] == ""){
                        // dynamic scaling
                        *d[systBranchName.c_str()]=*d[oldBranchName.c_str()] * *d[cursyst->scaleVar[iBrnch]];
                    } else if (cursyst->scaleVarRef[iBrnch] == "None"){
                        *d[systBranchName.c_str()] = *d[cursyst->scaleVar[iBrnch]];
                    } else {
                        *d[systBranchName.c_str()] = *d[oldBranchName.c_str()] * (*d[cursyst->scaleVar[iBrnch]]/ *d[cursyst->scaleVarRef[iBrnch]]);
                    }
                }
                // copy the value to the new branch
               // *d[systBranchName.c_str()]=*d[oldBranchName.c_str()];
            } else if(thisType==8){
                //std::cout<<"length branch "<<existingBranchInfo->lengthBranch<<std::endl;
                //std::cout<<"length "<<*in[existingBranchInfo->lengthBranch]<<std::endl;
                for(int ind=0; ind<*in[existingBranchInfo->lengthBranch]; ind++){// scale the current branch
                    //std::cout<<"Jet pt, eta: "<<f["Jet_pt_reg"][ind]<<", "<<f["Jet_eta"][ind]<<std::endl;
                    //std::cout<<"old val "<<f[oldBranchName.c_str()][ind]<<std::endl;
                    // FIXME: we are assuming that the only array variables we scale are jet variables, we should in principle make this more generic
                        double jet_reg = m("Jet_bReg",ind);
                    if (ptSplit==0 || (ptSplit==1 && jet_reg<jetPtSplit) || (ptSplit==2 && jet_reg>jetPtSplit)) {
                        if (etaSplit==0 || (etaSplit==1 && fabs(f["Jet_eta"][ind])<jetEtaSplit) || (etaSplit==2 && fabs(f["Jet_eta"][ind])>jetEtaSplit)) {
                            //std::cout<<"got through"<<std::endl;
                            if (cursyst->scaleVar[iBrnch] == "") {
                                // flat scaling
                                f[systBranchName.c_str()][ind]=f[oldBranchName.c_str()][ind] * cursyst->scales[iBrnch];
                            }
                            else {
                                if (cursyst->scaleVarRef[iBrnch] == ""){
                                    // dynamic scaling
                                    f[systBranchName.c_str()][ind]=f[oldBranchName.c_str()][ind] * f[cursyst->scaleVar[iBrnch]][ind];
                                } else if (cursyst->scaleVarRef[iBrnch] == "None"){
                                    f[systBranchName.c_str()][ind] = f[cursyst->scaleVar[iBrnch]][ind];
                                } else {
                                    f[systBranchName.c_str()][ind] = f[oldBranchName.c_str()][ind] * (f[cursyst->scaleVar[iBrnch]][ind]/ f[cursyst->scaleVarRef[iBrnch]][ind]);
                                }
                            }
                            //std::cout<<"new val "<<f[oldBranchName.c_str()][ind]<<std::endl;
                            // copy the value to the new branch
                            //f[systBranchName.c_str()][ind]=f[oldBranchName.c_str()][ind];
                            //std::cout<<"new branch "<<f[systBranchName.c_str()][ind]<<std::endl;
                        }
                    }
                }
            } else if(thisType==9){
                for(int ind=0; ind<*in[existingBranchInfo->lengthBranch]; ind++){// scale the current branch
                    if (cursyst->scaleVar[iBrnch] == "") {
                         // flat scaling
                        d[systBranchName.c_str()][ind]=d[oldBranchName.c_str()][ind] * cursyst->scales[iBrnch];
                    }
                    else {
                        if (cursyst->scaleVarRef[iBrnch] == ""){
                            // dynamic scaling
                            d[systBranchName.c_str()][ind]=d[oldBranchName.c_str()][ind] * d[cursyst->scaleVar[iBrnch]][ind];
                        } else if (cursyst->scaleVarRef[iBrnch] == "None"){
                            d[systBranchName.c_str()][ind] = d[cursyst->scaleVar[iBrnch]][ind];
                        } else {
                            d[systBranchName.c_str()][ind] = d[oldBranchName.c_str()][ind] * (d[cursyst->scaleVar[iBrnch]][ind]/ d[cursyst->scaleVarRef[iBrnch]][ind]);
                        }
                    }
                    // copy the value to the new branch
                    //d[systBranchName.c_str()][ind]=d[oldBranchName.c_str()][ind];
                }
            }
        }
    }

    //Loop again to set the nominal value of the branch to the updated value just for analyze() to work correctly
    for(unsigned iBrnch=0; iBrnch<cursyst->branchesToEdit.size(); iBrnch++){
        std::string oldBranchName(cursyst->branchesToEdit[iBrnch]);
        BranchInfo* existingBranchInfo = branchInfos[oldBranchName.c_str()];
        if(!early || (early && branchInfos[cursyst->branchesToEdit[iBrnch]]->prov == "early")){
            std::string systBranchName(oldBranchName);
            systBranchName.append("_");
            systBranchName.append(cursyst->name);
            int thisType=existingBranchInfo->type;
            if(thisType==2){
                *f[oldBranchName.c_str()]=*f[systBranchName.c_str()]; 
            } else if(thisType==2){
                *d[oldBranchName.c_str()]=*d[systBranchName.c_str()]; 
            } else if(thisType==8 ){
                //For Jet_bReg we have already correctly organised the various branches above
                for(int ind=0; ind<*in[existingBranchInfo->lengthBranch]; ind++){
                    f[oldBranchName.c_str()][ind]=f[systBranchName.c_str()][ind];
                }
            } else if(thisType==9){
                for(int ind=0; ind<*in[existingBranchInfo->lengthBranch]; ind++){
                    d[oldBranchName.c_str()][ind]=d[systBranchName.c_str()][ind];
                }
            }
        }
    }
}

double   AnalysisManager::GetPUWeight(int thisPU, int puType){
    float thisWeight=1.0;
    if(thisPU<0){
        if(debug>10){
            std::cout<<"thisPU is negative..."<<thisPU<<std::endl;
        }
        return thisWeight;
    }
    
    TH1D** pointerToReweightingPointer=NULL;

    if(puType==0){
        pointerToReweightingPointer=&PUReWeighting;
    }else if(puType==-1){
        pointerToReweightingPointer=&PUReWeightingDOWN;
    }else if(puType==1){
        pointerToReweightingPointer=&PUReWeightingUP;
    }else{
        std::cout<<"GetPUWeight not set for puType "<<puType<<std::endl;
        return 1.0;
    }

    if((*pointerToReweightingPointer)==NULL){
        if(debug>10){
            std::cout<<"PU reweighting histogram is NULL."<<std::endl;
        }
        return thisWeight;
    }

    unsigned int nBins = (*pointerToReweightingPointer)->GetNbinsX();
    float thisLowEdge=-1;
    float nextLowEdge=-1;
    int bestBin=-1;
    for(unsigned int iBin=0;iBin<nBins; iBin++){
        thisLowEdge=(*pointerToReweightingPointer)->GetBinLowEdge(iBin);
        nextLowEdge=(*pointerToReweightingPointer)->GetBinLowEdge(iBin+1);
        if(thisPU>=thisLowEdge && thisPU<nextLowEdge){
            bestBin=iBin;
            break;
        }
    }
    if(bestBin!=-1){
        if(debug>100000) std::cout<<"thisPU last thisLowEdge nextLowEdge "<<thisPU<<" "<<thisLowEdge<<" "<<nextLowEdge<<std::endl;
        thisWeight=(*pointerToReweightingPointer)->GetBinContent(bestBin);
    } else {
        if(debug>100000) std::cout<<"bestBin is -1.  PU is "<<thisPU<<std::endl;
    }
    return thisWeight;
}

bool     AnalysisManager::EvaluatePUReweighting(TH1D* inputPU, int puType){
    bool success=false;
    TH1D** pointerToTargetPointer=NULL;

    if(puType==0){
        pointerToTargetPointer=&globalPUTarget;
    }else if(puType==-1){
        pointerToTargetPointer=&globalPUTargetDOWN;
    }else if(puType==1){
        pointerToTargetPointer=&globalPUTargetUP;
    }else{
        std::cout<<"EvaluatePUReweighting not set for puType "<<puType<<std::endl;
        return success;
    }

    if((*pointerToTargetPointer)==NULL || inputPU==NULL){
        if(debug>10) std::cout<<"(*pointerToTargetPointer) or inputPU is NULL"<<std::endl;
    } else {
        inputPU->Scale(1./inputPU->Integral());
       
        if((*pointerToTargetPointer)->GetNbinsX()!=inputPU->GetNbinsX()){
            int maxBins=std::max((*pointerToTargetPointer)->GetNbinsX(),inputPU->GetNbinsX());
            float minPU=std::min((*pointerToTargetPointer)->GetBinLowEdge(1),inputPU->GetBinLowEdge(1));
            float maxPU=std::max((*pointerToTargetPointer)->GetBinLowEdge((*pointerToTargetPointer)->GetNbinsX()+1),inputPU->GetBinLowEdge(inputPU->GetNbinsX()+1));
            if(debug>10) std::cout<<"maxBins "<<maxBins<<"  maxPU  "<<maxPU<<"  minPU  "<<minPU<<std::endl;
            if((*pointerToTargetPointer)->GetNbinsX()!=maxBins || (*pointerToTargetPointer)->GetBinLowEdge(1)!=minPU || (*pointerToTargetPointer)->GetBinLowEdge((*pointerToTargetPointer)->GetNbinsX()+1)!=maxPU){
                TH1D * PUTargetReplacement= new TH1D("PUTarget",";nPU:Fraction",maxBins,minPU,maxPU);
                for(int iBin=0;iBin<maxBins+1; iBin++){
                    if(iBin<=(*pointerToTargetPointer)->GetNbinsX()+1){
                        PUTargetReplacement->Fill(iBin,(*pointerToTargetPointer)->GetBinContent(iBin));
                    } else {
                        PUTargetReplacement->Fill(iBin,(*pointerToTargetPointer)->GetBinContent((*pointerToTargetPointer)->GetNbinsX()+1));
                    }
                }
                (*pointerToTargetPointer)=(TH1D*)PUTargetReplacement->Clone();
                (*pointerToTargetPointer)->SetDirectory(0);
                (*pointerToTargetPointer)->Scale(1./(*pointerToTargetPointer)->Integral());
            }
            
            if(inputPU->GetNbinsX()!=maxBins || inputPU->GetBinLowEdge(1)!=minPU || inputPU->GetBinLowEdge(inputPU->GetNbinsX()+1)!=maxPU){
                TH1D * PUInputReplacement= new TH1D("PUInput",";nPU:Fraction",maxBins,minPU,maxPU);
                for(int iBin=0;iBin<maxBins+1; iBin++){
                    if(iBin<=inputPU->GetNbinsX()+1){
                        PUInputReplacement->Fill(iBin,inputPU->GetBinContent(iBin));
                    } else {
                        PUInputReplacement->Fill(iBin,inputPU->GetBinContent(inputPU->GetNbinsX()+1));
                    }
                }
                inputPU=(TH1D*)PUInputReplacement->Clone();
                inputPU->SetDirectory(0);
                inputPU->Scale(1./inputPU->Integral());
            }
        }
        TH1D** pointerToReweightingPointer=NULL;

        if(puType==0){
            pointerToReweightingPointer=&PUReWeighting;
        }else if(puType==-1){
            pointerToReweightingPointer=&PUReWeightingDOWN;
        }else if(puType==1){
            pointerToReweightingPointer=&PUReWeightingUP;
        }else{
            std::cout<<"EvaluatePUReweighting not set for puType "<<puType<<std::endl;
            return success;
        }

        (*pointerToReweightingPointer)=(TH1D*)(*pointerToTargetPointer)->Clone("PUReWeighting");
        (*pointerToReweightingPointer)->Divide(inputPU);
       
        if(debug>100){ 
            for(int iBin=0;iBin<(*pointerToTargetPointer)->GetNbinsX()+1; iBin++){
                std::cout<<"iBin "<<iBin<<" input "<<inputPU->GetBinContent(iBin)<<" target "<<(*pointerToTargetPointer)->GetBinContent(iBin)<<" weight "<<PUReWeighting->GetBinContent(iBin)<<std::endl; 
            }
        }
        success=true;
    }
    return success;
}

void     AnalysisManager::SetGlobalPUTarget(TH1D targetHist, int puType){
    TH1D** pointerToTargetPointer=NULL;

    if(puType==0){
        pointerToTargetPointer=&globalPUTarget;
    }else if(puType==-1){
        pointerToTargetPointer=&globalPUTargetDOWN;
    }else if(puType==1){
        pointerToTargetPointer=&globalPUTargetUP;
    }else{
        std::cout<<"SetGlobalPUTarget not set for puType "<<puType<<std::endl;
        return;
    }

    (*pointerToTargetPointer)=(TH1D*)targetHist.Clone("globalPUTarget"); 
    (*pointerToTargetPointer)->SetDirectory(0);
    (*pointerToTargetPointer)->Scale(1./(*pointerToTargetPointer)->Integral());
}

void     AnalysisManager::SetGlobalPUInput(TH1D inputHist){
    globalPUInput=(TH1D*)inputHist.Clone("globalPUInput"); 
    globalPUInput->SetDirectory(0);
    globalPUInput->Scale(1./globalPUInput->Integral());
}


double AnalysisManager::m(std::string key, int index){
    if(debug>100000) std::cout<<"looking for key "<<key<<" with index "<<index<<std::endl;
    if(branchInfos.count(key)==0){
        std::cout<<"There is no branch with name "<<key<<std::endl;
        if(safemode){
            std::cout<<"The program must be terminated..."<<std::endl;
            std::exit(0);
        } else {
            if(debug>1) std::cout<<"Returning -999 and hoping for the best."<<std::endl;
            return -999;
        }
    } else {
        if(debug>100000) std::cout<<"here is where I should return the right branch value"<<std::endl;
        if(branchInfos[key]->type%6==5){
            std::cout<<"Key "<<key<<" of type "<<branchInfos[key]->type<<", please use AnalysisManager::mInt() instead."<<std::endl;
            std::cout<<"Exiting"<<std::endl;
            std::exit(0); 
        }
        if(branchInfos[key]->type>5&&index<0){
            std::cout<<"No valid index ("<<index<<") specified for branch "<<key<<std::endl;
            std::cout<<"Exiting..."<<std::endl;
            std::exit(0);
        }
        switch(branchInfos[key]->type)
        {
        case 0:
            if(debug>100000) std::cout<<"unsigned int "<<*ui[key]<<std::endl;
            return (double)*ui[key];
        case 1:
            if(debug>100000) std::cout<<"int "<<*in[key]<<std::endl;
            return (double)*in[key];
        case 2:
            if(debug>100000) std::cout<<"float "<<*f[key]<<std::endl;
            return (double)*f[key];
        case 3:
            if(debug>100000) std::cout<<"double "<<*d[key]<<std::endl;
            return (double)*d[key];
        case 4:
            if(debug>100000) std::cout<<"boolean "<<*b[key]<<std::endl;
            return (double)*b[key];
        case 6:
            if(debug>100000) std::cout<<"unsigned int "<<ui[key][index]<<std::endl;
            return (double)ui[key][index];
        case 7:
            if(debug>100000) std::cout<<"int "<<in[key][index]<<std::endl;
            return (double)in[key][index];
        case 8:
            if(debug>100000) std::cout<<"float "<<f[key][index]<<std::endl;
            return (double)f[key][index];
        case 9:
            if(debug>100000) std::cout<<"double "<<d[key][index]<<std::endl;
            return (double)d[key][index];
        case 10:
            if(debug>100000) std::cout<<"boolean "<<b[key][index]<<std::endl;
            return (double)b[key][index];
        default:
            if(debug>10) std::cout<<"I don't know type "<<branchInfos[key]->type<<" yet..."<<std::endl;
            return -999;
        }
    }
}

int AnalysisManager::mInt(std::string key, int index){
    if(debug>100000) std::cout<<"looking for key "<<key<<" with index "<<index<<std::endl;
    if(branchInfos.count(key)==0){
        std::cout<<"There is no branch with name "<<key<<std::endl;
        if(safemode){
            std::cout<<"The program must be terminated..."<<std::endl;
            std::exit(0);
        } else {
            if(debug>1) std::cout<<"Returning -999 and hoping for the best."<<std::endl;
            return -999;
        }
    } else {
        if(debug>100000) std::cout<<"here is where I should return the right branch value"<<std::endl;
        if(branchInfos[key]->type%6==2 || branchInfos[key]->type%6==3){
            std::cout<<"Key "<<key<<" of type "<<branchInfos[key]->type<<", please use AnalysisManager::m() instead."<<std::endl;
            std::cout<<"Exiting"<<std::endl;
            std::exit(0); 
        }
        if(branchInfos[key]->type>5&&index<0){
            std::cout<<"No valid index specified for branch "<<key<<std::endl;
            std::cout<<"Exiting..."<<std::endl;
            std::exit(0);
        }
        switch(branchInfos[key]->type)
        {
        case 0:
            if(debug>100000) std::cout<<"unsigned int "<<*ui[key]<<std::endl;
            return (int)*ui[key];
        case 1:
            if(debug>100000) std::cout<<"int "<<*in[key]<<std::endl;
            return (int)*in[key];
        case 4:
            if(debug>100000) std::cout<<"boolean "<<*b[key]<<std::endl;
            return (int)*b[key];
        case 5:
            if(debug>100000) std::cout<<"char "<<*uc[key]<<std::endl;
            return (int)*uc[key]; //character form of an int always has '0' appended
        case 6:
            if(debug>100000) std::cout<<"unsigned int "<<ui[key][index]<<std::endl;
            return (int)ui[key][index];
        case 7:
            if(debug>100000) std::cout<<"int "<<in[key][index]<<std::endl;
            return (int)in[key][index];
        case 10:
            if(debug>100000) std::cout<<"boolean "<<b[key][index]<<std::endl;
            return (int)b[key][index];
        case 11:
            if(debug>100000) std::cout<<"char "<<uc[key][index]<<std::endl;
            return (int)uc[key][index];
        default:
            if(debug>10) std::cout<<"I don't know type "<<branchInfos[key]->type<<" yet..."<<std::endl;
            return -999;
        }
    }
}


double AnalysisManager::EvalDeltaPhi(double phi0, double phi1){
    double dPhi = fabs(phi0-phi1);
    //std::cout<<"dPhi PI "<<dPhi<<" "<<PI<<std::endl;

    if(dPhi > PI)
        dPhi = 2.0*PI - dPhi;

    return dPhi;
}


double AnalysisManager::EvalDeltaR(double eta0, double phi0, double eta1, double phi1)
{
  double dEta = fabs(eta0-eta1);
  double dPhi = EvalDeltaPhi(phi0, phi1);


  return TMath::Sqrt(TMath::Power(dEta,2)+TMath::Power(dPhi,2));
}
