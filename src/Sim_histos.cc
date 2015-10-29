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
#include "GretinaSim.hh"
#include "Mode3Calc.hh"
#include "Scaler.hh"
#include "SimHistograms.hh"
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
  char* TreeName = (char*)"ctr";
  char* SimTreeName = (char*)"str";
  int nmax =0;
  int vl =0;
  CommandLineInterface* interface = new CommandLineInterface();

  interface->Add("-i", "inputfiles", &InputFiles);
  interface->Add("-o", "outputfile", &OutputFile);
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

  TChain* tr  = new TChain(TreeName);
  TChain* str = new TChain(SimTreeName);
  vector<RunInfo*> info;
  vector<int> events;
  info.resize(InputFiles.size());
  events.resize(InputFiles.size()+1);
  events[0] = 0;
  TFile *dummy;
  Settings *set = NULL;
  double runtime=0;
  for(unsigned int i=0; i<InputFiles.size(); i++){
    dummy = new TFile(InputFiles[i]);
    tr->Add(InputFiles[i]);
    str->Add(InputFiles[i]);
    info[i] = (RunInfo*)dummy->Get("runinfo");
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
  if(vl>0)
    set->PrintSettings();
  
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

  if(str == NULL){
    cout << "could not find tree str in file " << endl;
    for(unsigned int i=0; i<InputFiles.size(); i++){
      cout<<InputFiles[i]<<endl;
    }
    return 3;
  }
  GretinaSim* gs = new GretinaSim;
  str->SetBranchAddress("GretinaSim",&gs);

  Double_t nentries = tr->GetEntries();
  Double_t nsimentries = str->GetEntries();

  cout << "creating outputfile " << endl;
  TFile* outfile = new TFile(OutputFile,"recreate");
    
  if(outfile->IsZombie()){
    return 4;
  }

  Int_t nbytes = 0;
  Int_t status = 1;
  
  SimHistograms* hists = new SimHistograms(set,nsimentries);

  cout << nsimentries << " entries in sim tree " << endl;
  cout << nentries << " entries in cal tree " << endl;

  int j = 0; // Cal tree index
  if(nmax>0)
    nsimentries = nmax;
  for(int i=0; i<nsimentries;i++){ // Sim tree index
    if(signal_received){
      break;
    }
    gr->Clear(); s800->Clear(); m3c->Clear(); gs->Clear();

    if(vl>2)
      cout << "getting entry " << i << endl;
    status = str->GetEvent(i);
    if(vl>2)
      cout << "status " << status << endl;
    if(status == -1){
      cerr<<"Error occured, couldn't read entry "<<i<<" from tree "<<str->GetName()<<" in file "<<str->GetFile()->GetName()<<endl;
      return 5;
    }
    else if(status == 0){
      cerr<<"Error occured, entry "<<i<<" in tree "<<str->GetName()<<" in file "<<str->GetFile()->GetName()<<" doesn't exist"<<endl;
      return 6;
    }
    nbytes += status;
    if(j < nentries){
      status = tr->GetEvent(j);
      if(vl>2)
    	cout << "status " << status << endl;
      if(status == -1){
    	cerr<<"Error occured, couldn't read entry "<<j<<" from tree "<<tr->GetName()<<" in file "<<tr->GetFile()->GetName()<<endl;
    	return 5;
      }
      else if(status == 0){
    	cerr<<"Error occured, entry "<<j<<" in tree "<<tr->GetName()<<" in file "<<tr->GetFile()->GetName()<<" doesn't exist"<<endl;
    	return 6;
      }
      nbytes += status;
    }

    hists->AddEntry();
    hists->FillHistograms(gr,s800,m3c, gs);

    // Synchronize the entries in the cal and sim trees
    if( gr->GetHit(0)->GetTS() == gs->GetTimeStamp() ){
      j++;
    }

    if(i%10000 == 0){
      double time_end = get_time();
      cout << setw(5) << setiosflags(ios::fixed) << setprecision(1) << (100.*i)/nsimentries <<
	" % done\t" << (Float_t)i/(time_end - time_start) << " events/s " << 
	(nsimentries-i)*(time_end - time_start)/(Float_t)i << "s to go \r" << flush;
    }
  }
  cout << endl;  

  cout << "writing histograms to file" << endl;
  outfile->cd();
  hists->Write();
  outfile->Close();
  delete tr;
  delete str;

  double time_end = get_time();
  cout << "Run time " << time_end - time_start << " s." << endl;
  
  return 0;
}
