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
#ifndef __TRACKSETTINGS_HH
#define __TRACKSETTINGS_HH

#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include "Settings.hh"

#include "TSystem.h"
#include "TEnv.h"

using namespace std;

class TrackSettings : public Settings {
public:
  TrackSettings();
  TrackSettings(const char* filename);
  TrackSettings(vector<char*> files);
  ~TrackSettings();
  void ReadTrackSettings(TEnv* set);
  void PrintTrackSettings();
  int MaxBad(){
    return fmaxBad;
  }
  bool RedoMap(){
    return fredoMap;
  }
  double JumpFOM(){
    return fjumpFOM;
  }
  int Probabilities(){
    return fprob;
  }
  double RadiusGE(){
    return fradiusGE;
  }
protected:
  float fjumpFOM;				
  int fmaxBad;
  int fredoMap;
  int fprob;
  float fradiusGE;				

  ClassDef(TrackSettings, 1)
};

#endif
