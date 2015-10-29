#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TROOT.h"
#include "/home/wimmer/GrROOT/src/RunInfo.hh"


void doit(){  
  //gROOT->Reset();
  char* InputFiles[5];
  InputFiles[0] = (char*)"cal/cal0095.root";
  InputFiles[1] = (char*)"cal/cal0096.root";
  InputFiles[2] = (char*)"cal/cal0097.root";
  InputFiles[3] = (char*)"cal/cal0098.root";
  InputFiles[4] = (char*)"cal/cal0099.root";
  TChain* tr = new TChain("ctr");
  RunInfo* info[5];
  int events[6];
  events[0] = 0;
  double runtime=0;
  TFile *dummy;
  for(unsigned int i=0; i<5; i++){
    dummy = new TFile(InputFiles[i]);
    tr->Add(InputFiles[i]);
    info[i] = (RunInfo*)dummy->Get("runinfo");
    runtime += info[i]->GetRunTime();
    events[i+1] = 0;
    events[i+1] = tr->GetEntries(); //cumulative!!
    cout << runtime <<"\t"<< events[i] <<"\t"<< events[i+1] << endl;
    delete dummy;
  }

}
void emptyfile(){
  ofstream file("settings/emptyevt.dat");
  for(int i=0;i<400;i++){
    file << "Start.Event.Number."<<i<<":\t0" << endl;
  }
  

}
