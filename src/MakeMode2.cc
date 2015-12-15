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
#include <iostream>
#include <iomanip>
#include <string.h>
#include <sys/time.h>
#include <signal.h>


#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH2S.h"
#include "TH1S.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TMath.h"
#include "TCutG.h"
#include "TEnv.h"
#include "TKey.h"

#include "CommandLineInterface.hh"
#include "S800.hh"
#include "Gretina.hh"
#include "Trace.hh"
#include "Scaler.hh"
#include "Settings.hh"

using namespace TMath;
using namespace std;


bool signal_received = false;
void signalhandler(int sig){
  if (sig == SIGINT){
    signal_received = true;
  }
}

double get_time(){
    struct timeval t;
    gettimeofday(&t, NULL);
    double d = t.tv_sec + (double) t.tv_usec/1000000;
    return d;
}

int main(int argc, char* argv[]){
  signal(SIGINT,signalhandler);

  double time_start = get_time();

  char* InputFile = NULL;
  char* OutputFile = NULL;
  vector<char*> SettingFile;
  int lastEvent = -1;

  //Read in the command line arguments
  CommandLineInterface* interface = new CommandLineInterface();
  interface->Add("-i", "input file", &InputFile);
  interface->Add("-o", "output file", &OutputFile);
  interface->Add("-s", "settingsfile", &SettingFile);
  interface->Add("-n", "last event", &lastEvent);
  interface->CheckFlags(argc, argv);

  //Complain about missing mandatory arguments
  if(InputFile == NULL){
    cout << "No input file given " << endl;
    return 1;
  }
  if(OutputFile == NULL){
    cout << "No output file given " << endl;
    return 2;
  }
  if(SettingFile.size() == 0){
    cout << "No settings file given " << endl;
    return 3;
  }

  Settings* set = new Settings(SettingFile);
  TFile* file = new TFile(InputFile);
  TTree* tr = (TTree*)file->Get("gtr");

  if(tr == NULL){
    cout << "could not find tree gtr in file" << endl;
    return 1;
  }

  Gretina* gr = new Gretina;
  S800* s800 = new S800;
  Mode3Event* m3r = new Mode3Event;
  tr->SetBranchAddress("gretina",&gr);
  tr->SetBranchAddress("s800",&s800);
  tr->SetBranchAddress("mode3Event",&m3r);

  Double_t nentries = tr->GetEntries();
  if(lastEvent>0 && lastEvent<nentries){
    nentries = lastEvent;
  }


  //Put before-loop actions here.
  FILE* outfile = fopen(OutputFile,"wb");
  unsigned int header[4];

  for(int i=0; i<nentries;i++){
    if(signal_received){
      break;
    }
    tr->GetEvent(i);

    //Put in-loop actions here.
    for(int crys=0; crys<m3r->GetMult(); crys++){
      crys_ips_abcd5678 obj = m3r->GetHit(crys)->MakeMode2(set);
      header[0] = GRETINA_ID;
      header[1] = sizeof(crys_ips_abcd5678);
      header[2] = obj.timestamp; //Low 32 bits
      header[3] = obj.timestamp >> 32; //High 32 bits
      fwrite(header,sizeof(unsigned int),4,outfile);
      fwrite(&obj,sizeof(crys_ips_abcd5678),1,outfile);
    }

    if(i%1000 == 0){
      double time_end = get_time();
      cout<<setw(5)<<setiosflags(ios::fixed)<<setprecision(1)<<(100.*i)/nentries<<" % done\t"<<(Float_t)i/(time_end - time_start)<<" events/s " << (nentries-i)*(time_end - time_start)/(Float_t)i<<"s to go \r"<<flush;
    }

  }
  cout << endl;

  //Put after-loop actions here.
  fclose(outfile);

  delete tr;

  double time_end = get_time();
  cout << "Program run time " << time_end - time_start << " s." << endl;

  return 0;
}
