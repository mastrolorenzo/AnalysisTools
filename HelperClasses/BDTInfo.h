#ifndef BDTInfo_h
#define BDTInfo_h

#include "TMVA/Reader.h"

#include "BDTVariable.h"

class BDTInfo{

public:
  std::string bdtname;
  std::string methodName;
  std::vector<BDTVariable> bdtVars;
  std::string xmlFile;
  TMVA::Reader *reader;
  std::string mvaType;

  BDTInfo(std::string, std::string, std::string, std::string="BDT");
  BDTInfo(BDTInfo&);
  BDTInfo();
  void AddVariable(std::string, std::string, bool, bool);
  void BookMVA();

};

#endif
