// code template from https://github.com/h2gglobe/h2gglobe

#ifndef SYSTEMATICCONTAINER
#define SYSTEMATICCONTAINER

#include <string>
#include <map>
#include <vector>

#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>

class SystematicContainer {
	
 public:	
    SystematicContainer();
    ~SystematicContainer();
 
    bool apply;
    std::string name;
    std::vector<std::string> branchesToEdit;
    std::vector<float> scales;  // scale value by
    std::vector<float> smears;  // smear value (after scaling)
    std::vector<std::string> scaleVar; // dynamic scaling, JER's for example
    std::vector<std::string> scaleVarRef; // Reference for calculating dynamic scaling
    std::map<std::string,float> gtcuts; // only apply systematic effect if branchname > cutvalue
    std::map<std::string,float> ltcuts; // only apply systematic effect if branchname < cutvalue

    // anticipating tables of eta/pt
    TH2F scaleTable;
    TH2F smearTable;

    void TurnOff();

    void AddBranchName(std::string brnchName);
    void AddScale(float newScale);
    void AddScaleVar(std::string newScaleVar);
    void AddScaleVarRef(std::string newScaleVarRef);
    void AddSmear(float newSmear);
    
};

#endif
