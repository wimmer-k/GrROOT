#ifndef __SIMSETTINGS_HH
#define __SIMSETTINGS_HH

#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <vector>
#include "Settings.hh"

#include "TSystem.h"
#include "TEnv.h"

using namespace std;

class SimSettings : public Settings {
public:
  SimSettings();//default ctor
  SimSettings(const char* filename);
  SimSettings(vector<char*> files);
  ~SimSettings();

  void ReadSimSettings(TEnv* set);
  void PrintSimSettings();
  const char* SimResolutionFile(){return fResFile.c_str();}
  const char* SimThresholdFile(){return fThreshFile.c_str();}
  double SimGretinaPositionResolution(){return fGretPosRes;}
  double SimS800AngleResolution(){return fS800AngleRes;}
  double SimS800YTAResolution(){return fS800YTARes;}

protected:
  string fResFile;
  string fThreshFile;
  double fGretPosRes;
  double fS800AngleRes;
  double fS800YTARes;

  ClassDef(SimSettings, 1)
};

#endif
