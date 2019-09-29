#include "BDTInfo.h"

BDTInfo::BDTInfo(std::string _methodName, std::string _bdtname, std::string _xmlFile, std::string _mvaType, unsigned int _nOutputs) {
    mvaType = _mvaType;
    bdtname = _bdtname;
    methodName = _methodName;
    bdtVars = std::vector<BDTVariable>();
    //method = _method;
    xmlFile = _xmlFile;
    reader = std::make_unique<TMVA::Reader>( "!Color:!Silent:Error" );
    reader->SetVerbose(kTRUE);
    mostProbIndex= _bdtname+"MaxInd";
    nOutputs=_nOutputs;
}

BDTInfo::BDTInfo(BDTInfo& _bdtInfo) {
    mvaType = _bdtInfo.mvaType;
    bdtname = _bdtInfo.bdtname;
    methodName = _bdtInfo.methodName;
    mostProbIndex = _bdtInfo.mostProbIndex;
    nOutputs = _bdtInfo.nOutputs;
    bdtVars = std::vector<BDTVariable>();
    for(unsigned int iVar=0; iVar<_bdtInfo.bdtVars.size(); iVar++){
        bdtVars.push_back(_bdtInfo.bdtVars[iVar]);
    }
    xmlFile = _bdtInfo.xmlFile;
    reader = std::make_unique<TMVA::Reader>( "!Color:!Silent:Error" );
    reader->SetVerbose(kTRUE);
}

BDTInfo::BDTInfo() {
    BDTInfo("", "", "");
}

void BDTInfo::AddVariable(std::string varName, std::string localVarName, bool isExisting, bool isSpec) {
    bdtVars.push_back(BDTVariable(varName, localVarName, isExisting, isSpec));
}

void BDTInfo::BookMVA(){
    reader->BookMVA(methodName,xmlFile);
}
