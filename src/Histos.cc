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
#include "RawHistograms.hh"
using namespace TMath;
using namespace std;



const int raw_det2hole[7] = {28,25,29,24,19,20,9};
int raw_hole2det(int hole){
  if(hole==15) return 0;
  if(hole==6) return 1;
  if(hole==-1) return 2;
  if(hole==8) return 3;
  if(hole==9) return 4;
  if(hole==7) return 5;
  if(hole==17) return 6;

  return -1;
}

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
  vector<char*> InputFiles;
  char* OutputFile = NULL;
  vector<char*> SettingFile;
  int nmax =0;
  int vl =0;
  CommandLineInterface* interface = new CommandLineInterface();

  interface->Add("-i", "inputfiles", &InputFiles);
  interface->Add("-o", "outputfile", &OutputFile);
  interface->Add("-s", "settingsfile", &SettingFile);
  interface->Add("-n", "nmax", &nmax);
  interface->Add("-v", "verbose", &vl);
  interface->CheckFlags(argc, argv);

  if(InputFiles.size() == 0 || OutputFile == NULL){
    cerr<<"You have to provide at least one input file and the output file!"<<endl;
    exit(1);
  }
  cout<<"input file(s):"<<endl;
  for(unsigned int i=0; i<InputFiles.size(); i++){
    cout<<InputFiles[i]<<endl;
  }
  cout<<"output file: "<<OutputFile<< endl;

  TChain* tr;
  tr = new TChain("gtr");
  for(unsigned int i=0; i<InputFiles.size(); i++){
    tr->Add(InputFiles[i]);
  }
  //tr->Print("toponly");


  if(tr == NULL){
    cout << "could not find tree gtr in file " << endl;
    for(unsigned int i=0; i<InputFiles.size(); i++){
      cout<<InputFiles[i]<<endl;
    }
    return 3;
  }
  Gretina* gr = new Gretina;
  S800* s800 = new S800;
  Mode3Event* m3r = new Mode3Event;
  tr->SetBranchAddress("gretina",&gr);
  tr->SetBranchAddress("s800",&s800);
  tr->SetBranchAddress("mode3Event",&m3r);

  Double_t nentries = tr->GetEntries();

  cout << "creating outputfile " << endl;
  TFile* outfile = new TFile(OutputFile,"recreate");

  if(outfile->IsZombie()){
    return 4;
  }

  Int_t nbytes = 0;
  Int_t status;

  if(SettingFile.size() == 0){
    cout << "Please give a settings file to use" << endl;
    return 5;
  }
  RawHistograms* hists = new RawHistograms(new Settings(SettingFile),nentries);

  cout << nentries << " entries in tree " << endl;

  if(nmax>0)
    nentries = nmax;

  for(int i=0; i<nentries;i++){
    if(signal_received){
      break;
    }
    gr->Clear(); s800->Clear(); m3r->Clear();

    if(vl>2)
      cout << "getting entry " << i << endl;
    status = tr->GetEvent(i);
    if(vl>2)
      cout << "status " << status << endl;
    if(status == -1){
      cerr<<"Error occured, couldn't read entry "<<i<<" from tree "<<tr->GetName()<<" in file "<<tr->GetFile()->GetName()<<endl;
      return 5;
    }
    else if(status == 0){
      cerr<<"Error occured, entry "<<i<<" in tree "<<tr->GetName()<<" in file "<<tr->GetFile()->GetName()<<" doesn't exist"<<endl;
      return 6;
    }
    nbytes += status;

    hists->FillHistograms(gr,s800,m3r);

    if(i%1000 == 0){
      double time_end = get_time();
      cout<<setw(5)<<setiosflags(ios::fixed)<<setprecision(1)<<(100.*i)/nentries<<" % done\t"<<(Float_t)i/(time_end - time_start)<<" events/s " << (nentries-i)*(time_end - time_start)/(Float_t)i<<"s to go \r"<<flush;
    }

  }
  cout << endl;


  cout << "writing histograms to file" << endl;
  hists->Write();
  outfile->Close();
  delete tr;

  double time_end = get_time();
  cout << "Run time " << time_end - time_start << " s." << endl;

  return 0;
}
