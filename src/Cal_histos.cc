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
#include "TDirectory.h"

#include "CommandLineInterface.hh"
#include "S800Calc.hh"
#include "GretinaCalc.hh"
#include "Mode3Calc.hh"
#include "Scaler.hh"
#include "CalHistograms.hh"
#include "RunInfo.hh"
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
  vector<char*> InputFiles;
  char* OutputFile = NULL;
  char* CutFile = NULL;
  char* SettingFile = NULL;
  char* TreeName = (char*)"ctr";
  int tac = 0;
  int Cal = 0;
  int nmax =0;
  int ana =0;
  int vl =0;
  CommandLineInterface* interface = new CommandLineInterface();

  interface->Add("-i", "inputfiles", &InputFiles);
  interface->Add("-o", "outputfile", &OutputFile);
  interface->Add("-s", "settingsfile", &SettingFile);
  interface->Add("-tn", "name of the tree", &TreeName);
  interface->Add("-t", "0 for TAC, 1 for TDC (default = 0)", &tac);
  interface->Add("-c", "cutfile", &CutFile);
  interface->Add("-cal", "Histos for calibration 1 crdc cal, 2 for Ionchamber, 4 for TOF and IC corrections (and combinations i.e. 3 = crdc and ic)", &Cal);
  interface->Add("-n", "nmax", &nmax);
  interface->Add("-a", "analysis type: 0 all, 1 s800 single, 2 s800 coinc, 3 s800 single or coinc, 4 hodo only", &ana);
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
  tr = new TChain(TreeName);
  vector<RunInfo*> info;
  vector<int> events;
  info.resize(InputFiles.size());
  events.resize(InputFiles.size()+1);
  events[0] = 0;
  TFile *dummy;
  Settings *set = NULL;
  if(SettingFile != NULL)
    set = new Settings(SettingFile);
  double runtime=0;

  cout << "creating outputfile " << endl;
  TFile* outfile = new TFile(OutputFile,"recreate");
    
  if(outfile->IsZombie()){
    return 4;
  }

  for(unsigned int i=0; i<InputFiles.size(); i++){
    dummy = new TFile(InputFiles[i]);
    tr->Add(InputFiles[i]);
    info[i] = (RunInfo*)dummy->Get("runinfo");
    outfile->cd();
    info[i]->Write(Form("runinfo_%04d",info[i]->GetRunNr()),TObject::kOverwrite);
    runtime += info[i]->GetRunTime();
    events[i+1] = tr->GetEntries(); //cumulative!!
    if(set==NULL){
      set = (Settings*)dummy->Get("settings");
    }
    // else if(set != (Settings*)dummy->Get("settings")){
    //   cout << "settings have changed between runs " << InputFiles[i-1] << " and " << InputFiles[i] << " aborting" << endl;
    //   return 5;
    // }
    delete dummy;
  }
  if(set==NULL){
    cout << "couldn't find the settings in any input file" << endl;
    return 4;
  }
  cout << "total run time " << runtime << " s " << endl;

  if(tr == NULL){
    cout << "could not find tree ctr in file " << endl;
    for(unsigned int i=0; i<InputFiles.size(); i++){
      cout<<InputFiles[i]<<endl;
    }
    return 3;
  }
  GretinaCalc* gr = new GretinaCalc;
  S800Calc* s800 = new S800Calc;
  Mode3Calc* m3c = new Mode3Calc;
  tr->SetBranchAddress("gretinacalc",&gr);
  tr->SetBranchAddress("s800calc",&s800);
  tr->SetBranchAddress("mode3calc",&m3c);

  Double_t nentries = tr->GetEntries();

  Int_t nbytes = 0;
  Int_t status;
  
  CalHistograms* hists = new CalHistograms(set,tac,Cal,CutFile,nentries);

  cout << nentries << " entries in tree " << endl;

  if(nmax>0)
    nentries = nmax;
  for(int i=0; i<nentries;i++){
    if(signal_received){
      break;
    }
    gr->Clear(); s800->Clear(); m3c->Clear();

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

    hists->AddEntry();

    if(ana == 1 && !s800->IsSingle())
      continue;
    if(ana == 2 && !s800->IsCoinc())
      continue;
    if(ana == 3 && !s800->IsCoinc() && !s800->IsSingle())
      continue;
    if(ana == 4 && !s800->IsHodo())
      continue;

    hists->FillHistograms(gr,s800,m3c);

    if(i%10000 == 0){
      double time_end = get_time();
      cout << setw(5) << setiosflags(ios::fixed) << setprecision(1) << (100.*i)/nentries <<
	" % done\t" << (Float_t)i/(time_end - time_start) << " events/s " << 
	(nentries-i)*(time_end - time_start)/(Float_t)i << "s to go \r" << flush;
    }
  }
  cout << endl; 

  cout << "writing histograms to file" << endl;
  outfile->cd();
  hists->Write();
  outfile->Close();
  delete tr;

  double time_end = get_time();
  cout << "Run time " << time_end - time_start << " s." << endl;
 
  return 0;
}
