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
#include "SimSettings.hh"

using namespace std;
SimSettings::SimSettings(){
}
SimSettings::SimSettings(const char* filename): Settings(filename){
  TEnv set(filename);
  ReadSimSettings(&set);
  if(fVerboseLevel>1)
    PrintSimSettings();
}
SimSettings::SimSettings(vector<char*> files): Settings(files){
  TEnv set;
  //Reverse order, because duplicate entries are ignored, instead of overwriting.
  for(vector<char*>::reverse_iterator it = files.rbegin(); it!=files.rend(); it++){
    set.ReadFile(*it,kEnvLocal);
    fInputFiles.push_back(*it);
  }
  ReadSimSettings(&set);
  if(fVerboseLevel>1)
    PrintSimSettings();
}
SimSettings::~SimSettings(){}
void SimSettings::ReadSimSettings(TEnv* set){
  char* defaultfile = (char*)"~/analysis/settings/nocal.dat";
  fResFile = set->GetValue("Sim.Resolution.File",defaultfile);
  fThreshFile = set->GetValue("Sim.Threshold.File",defaultfile);
  fS800AngleRes = set->GetValue("Sim.S800.Angle.Resolution",0.0);
  fS800YTARes = set->GetValue("Sim.S800.YTA.Resolution",0.0);
  fGretPosRes = set->GetValue("Sim.Gretina.Position.Resolution",0.0);

}
void SimSettings::PrintSimSettings(){
  cout << "Sim.Resolution.File\t" << fResFile << endl;
  cout << "Sim.Threshold.File\t" << fThreshFile << endl;
  cout << "Sim.Gretina.Position.Resolution\t" << fGretPosRes << endl;
  cout << "Sim.S800.Angle.Resolution\t" << fS800AngleRes << endl;
  cout << "Sim.S800.YTA.Resolution\t" << fS800YTARes << endl;
}
