////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////
////////////                       GrROOT
////////////
////////////          Purpose:
////////////                   To assist in the analysis of data from
////////////                 the gretina/S800 experimental setup.
////////////                          
////////////          Current Maintainers:
////////////                 Kathrin Wimmer  (wimmer@phys.s.u-tokyo.ac.jp)
////////////                 Eric Lunderberg (lunderberg@nscl.msu.edu)
////////////
////////////          Distribution:
////////////                   Please do not redistribute this software directly.
////////////                   If someone new wants a copy of this software,
////////////                 email one of the maintainers for the download link.
////////////                   This allows us to keep track of who has the software,
////////////                 and properly distribute updates and bug fixes.
////////////                 
////////////          Suggestions:
////////////                   We view the development of the software as a collaborative
////////////                 effort, and as such, welcome and appreciate any suggestions
////////////                 for bug fixes and improvements.
////////////
////////////          Disclaimer:
////////////                 This software is provided as-is, with no warranty.
////////////                 No current or future support is guaranteed for this software.
////////////
////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
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
