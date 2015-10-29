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
