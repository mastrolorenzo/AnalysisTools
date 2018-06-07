// code template from https://github.com/h2gglobe/h2gglobe

#include "SampleContainer.h"
#include <utility>
#include <iostream>

inline SampleContainer::~SampleContainer() 
{}

//float SampleContainer::defaultextw=1.;

//SampleContainer::SampleContainer(const float * extw) :
//	extweight(extw)
inline SampleContainer::SampleContainer() 
{
	/*intWeight = 1;
	itype = 0;
	ind = 0;
	histoplotit = 1;
	filesshortnam = "";
	ntot = 0;
	nred = 0;
	lumi = 0;
	lumireal = 1;
	hasLumiSelection = false;
	hasEventList = false;
	pileup = "";
	forceVersion = 0;*/
    sampleName = "";
    sampleNum = 0;
    sampleChain = new TChain("Events");
    //sampleChain = new TChain("tree");
    xsec = 0;
    kFactor= 1;
    scale = 1;
    files.clear();
    processedEvents=0.;
    intWeight=1;
    nProFromFile=false;
    doJetFlavorSplit = false;
    procEff = 1.;
    CountWeightedLHEWeightScale = new TH1F("CountWeightedLHEWeightScale","CountWeightedLHEWeightScale",9,-0.5,8.5);
    CountWeightedLHEWeightPdf = new TH1F("CountWeightedLHEWeightPdf","CountWeightedLHEWeightPdf",120,-0.5,119.5);
    CountWeighted = new TH1F("CountWeighted","CountWeighted",1,0.,2.0);
    CountFullWeighted = new TH1F("CountFullWeighted","CountFullWeighted",1,0.,2.0);
    InputPU = NULL;
    PUHistName = "";
    lepFlav = -1;
    forceOverWriteExtFile = false;
    externFileName = "";
    externFileExists = false;
}

inline void SampleContainer::AddFile(const char* fname,int isBatch, int doSkim) {
    files.push_back(fname);
    if( isBatch==1 ) return;
    
    sampleChain->Add(fname);
    //std::cout<<"nProFromFile "<<nProFromFile<<" doSkim "<<doSkim<<std::endl; 
    if(forceOverWriteExtFile || !externFileExists){
        bool openTFile=(nProFromFile || !(PUHistName.empty()));
        TFile *file;
   
        if(openTFile){
            file = TFile::Open(fname);
            if (file->IsZombie()) return;
        }
   
        if(!PUHistName.empty() && sampleNum != 0){
            if(doSkim==2) PUHistName+="skim";
            TH1D* thisFilePUHist=(TH1D*)((TH1D*)file->Get(PUHistName.c_str()))->Clone("thisFilePUHist");
            if(InputPU == NULL){
                if(doSkim==1){
                    PUHistName+="skim";
                    InputPU=(TH1D*)(thisFilePUHist->Clone(PUHistName.c_str()));
                } else {
                    InputPU=(TH1D*)(thisFilePUHist->Clone("PUTarget"));
                }
            } else {
                InputPU->Add(thisFilePUHist);
            }
            InputPU->SetDirectory(0);
            delete thisFilePUHist;
        }
        
        if(nProFromFile) {
            if (doSkim == 2 && files.size() > 1) return; // skimmed files already have the summed count histograms
            if (doSkim == 2) {
                // doSkim == 2 -> running on skimmed samples by AT
                // for now our skimmed ntuples we preserve the counting by histogram structure
                // from Heppy, since this is how the datacard maker is set up. We may want to eventually
                // change this.
                
                TH1F* counter = (TH1F*)file->Get("CountWeighted");
                float nEffective = counter->GetBinContent(1);
                std::cout<<"nEffective "<<nEffective<<std::endl;
                //if(sampleNum==49 or sampleNum==491) {
                //    //special prescription for WJets_BGenFilter sample weighting
                //    //TH1F* counterFullWeight = (TH1F*)file->Get("CountFullWeighted");
                //    nEffective = counterFullWeight->GetBinContent(1); 
                //}
                CountWeighted->Add(counter);
                std::cout<<"pe = "<<processedEvents<<std::endl;
                processedEvents += nEffective;
                std::cout<<"pe = "<<processedEvents<<std::endl;
                TH1F* CountWeightedLHEWeightScale_thisfile = (TH1F*)file->Get("CountWeightedLHEWeightScale");
                TH1F* CountWeightedLHEWeightPdf_thisfile = (TH1F*)file->Get("CountWeightedLHEWeightPdf");
                //std::cout<<"lhe = "<<CountWeightedLHEWeightPdf->GetBinContent(1)<<std::endl;
                CountWeightedLHEWeightScale->Add(CountWeightedLHEWeightScale_thisfile);
                CountWeightedLHEWeightPdf->Add(CountWeightedLHEWeightPdf_thisfile);
                //std::cout<<"lhe = "<<CountWeightedLHEWeightPdf->GetBinContent(1)<<std::endl;
            } else {
                // totally different setup for grabbing event count in nanoAOD
                TTree *Runs = (TTree*) file->Get("Runs");
                Double_t genEventSumw = 0;
                Double_t lheScales[9];
                Double_t lhePdfs[120];
                Runs->SetBranchAddress("genEventSumw",&genEventSumw);
                Runs->SetBranchAddress("LHEScaleSumw",&lheScales);
                Runs->SetBranchAddress("LHEPdfSumw",&lhePdfs);
                // one entry in Runs tree per input NanoAOD file
                int nNanoInputFiles = Runs->GetEntries();
                std::cout<<"From Runs tree, processedEvents is (before) "<<processedEvents<<std::endl;
                for (int i=0; i < nNanoInputFiles; i++) {
                    Runs->GetEntry(i);
                    std::cout<<fname<<" genEventSumw: "<<genEventSumw<<std::endl;
                    CountWeighted->SetBinContent(1,CountWeighted->GetBinContent(1)+genEventSumw);
                    processedEvents += genEventSumw;
                    for (int j=0; j<9; j++) {
                        CountWeightedLHEWeightScale->SetBinContent(j+1,CountWeightedLHEWeightScale->GetBinContent(j+1)+lheScales[j]);
                    }
                    for (int j=0; j<120; j++) {
                        CountWeightedLHEWeightPdf->SetBinContent(j+1,CountWeightedLHEWeightPdf->GetBinContent(j+1)+lhePdfs[j]);
                    }
                }
                std::cout<<"From Runs tree, processedEvents is (after)  "<<processedEvents<<std::endl;
            }
        }
        if(openTFile){
            file->Close();
        }
    }
}

inline void SampleContainer::ComputeWeight(float intL) {
    if(sampleNum==0) { //this is data
        intWeight = 1; 
    } else {
        //std::cout << "Computing Weight for type - " << sampleNum << ", Using " << processedEvents << " Processed Events" << std::endl;
        intWeight = (kFactor*scale*xsec*intL)/(processedEvents*procEff);
    }
}

inline void SampleContainer::CreateSampleInfoFile(){
    //TFile* sampleInfoFile= new TFile(externFileName.c_str(),"RECREATE");
    std::string sampleInfoFileName(sampleName);
    sampleInfoFileName.append("_sampleInfo.root");
    TFile* sampleInfoFile= new TFile(sampleInfoFileName.c_str(),"RECREATE");
    sampleInfoFile->cd();

    CountWeightedLHEWeightScale->Write();
    CountWeightedLHEWeightPdf->Write();
    CountWeighted->Write();
    InputPU->Write();
   
    // CountWeighted does this already.  Keep for reference...
    // Maybe a tree will be useful in the future.
    //TTree* sampleInfoTree = new TTree("sampleInfo","Sample Info Tree");
    //sampleInfoTree->Branch("processedEvents",&processedEvents,"processedEvents/f");
    //sampleInfoTree->Fill();
    //sampleInfoTree->Write();

    sampleInfoFile->Close();
    delete sampleInfoFile;
}

inline void SampleContainer::ReadSampleInfoFile(){
    TFile* sampleInfoFile = TFile::Open(externFileName.c_str());

    TH1F* CountWeightedLHEWeightScale_thisfile = (TH1F*)sampleInfoFile->Get("CountWeightedLHEWeightScale");
    TH1F* CountWeightedLHEWeightPdf_thisfile = (TH1F*)sampleInfoFile->Get("CountWeightedLHEWeightPdf");
    CountWeightedLHEWeightScale->Add(CountWeightedLHEWeightScale_thisfile);
    CountWeightedLHEWeightPdf->Add(CountWeightedLHEWeightPdf_thisfile);
    
    TH1F* counter = (TH1F*)sampleInfoFile->Get("CountWeighted");
    processedEvents = counter->GetBinContent(1);
            
    TH1D* thisFilePUHist=(TH1D*)((TH1D*)sampleInfoFile->Get("PUTarget"))->Clone("thisFilePUHist");
    InputPU=(TH1D*)(thisFilePUHist->Clone(PUHistName.c_str()));
    InputPU->SetDirectory(0);
    
    delete counter;
    delete thisFilePUHist;
    delete sampleInfoFile;
}

/* 
// ----------------------------------------------------------------------------------------------------------------------
void SampleContainer::addGoodLumi(int run, int lumi1, int lumi2 )
{
	hasLumiSelection = true;
	goodLumis[run].push_back( std::make_pair(lumi1,lumi2) );
}


// ----------------------------------------------------------------------------------------------------------------------
void SampleContainer::addEventToList(int run, int lumi, int event )
{
	hasEventList = true;
	eventList[run].push_back( std::make_pair(lumi,event) );
}*/
