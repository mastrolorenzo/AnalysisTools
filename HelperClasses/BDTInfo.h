#ifndef BDTInfo_h
#define BDTInfo_h

#include <memory>

#include "TMVA/Reader.h"

#include "BDTVariable.h"

class BDTInfo{

public:
  std::string bdtname;
  std::string methodName;
  std::vector<BDTVariable> bdtVars;
  std::string xmlFile;
  std::unique_ptr<TMVA::Reader> reader;
  std::string mvaType;
  std::string mostProbIndex;
  unsigned int nOutputs;

  BDTInfo(std::string, std::string, std::string, std::string="BDT", unsigned int=1);
  BDTInfo(BDTInfo&);
  BDTInfo();
  void AddVariable(std::string, std::string, bool, bool);
  void BookMVA();
  void DeleteReader() {delete reader.release();}

};

#endif
