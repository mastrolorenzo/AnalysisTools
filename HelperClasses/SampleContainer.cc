// code template from https://github.com/h2gglobe/h2gglobe

#include "HelperClasses/SampleContainer.h"
#include <utility>
#include <iostream>

SampleContainer::~SampleContainer() 
{}

//float SampleContainer::defaultextw=1.;

//SampleContainer::SampleContainer(const float * extw) :
//	extweight(extw)
SampleContainer::SampleContainer() 
{
	/*intweight = 1;
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
    sampleChain = new TChain("tree");
    xsec = 0;
    kfactor= 1;
    scale = 1;
    files.clear();
}

void SampleContainer::AddFile(char* fname) {
    //std::cout<<"fname "<<fname<<std::endl;
    sampleChain->Add(fname);
    files.push_back(fname);
}

/*void SampleContainer::computeWeight(float intL) {
  if(itype==0) { //this is data
    intweight = 1; 
  } else {
    std::cout << "Computing Weight for type - " << itype << ", Using " << ntot << " Processed Events" << std::endl;
    intweight = kfactor*scale*xsec*intL/ntot;
  }
}
 
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
