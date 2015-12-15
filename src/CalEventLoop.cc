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

#include "S800Calc.hh"
#include "GretinaCalc.hh"
#include "Mode3Calc.hh"

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
  char InputFile[100];
  strcpy(InputFile,argv[1]);

  TFile* file = new TFile(InputFile);
  TTree* tr = (TTree*)file->Get("ctr");

  if(tr == NULL){
    cout << "could not find tree ctr in file " << endl;
    return 1;
  }

  GretinaCalc* gr = new GretinaCalc;
  S800Calc* s800 = new S800Calc;
  Mode3Calc* m3c = new Mode3Calc;
  tr->SetBranchAddress("gretinacalc",&gr);
  tr->SetBranchAddress("s800calc",&s800);
  tr->SetBranchAddress("mode3calc",&m3c);

  Double_t nentries = tr->GetEntries();

  //Put before-loop actions here.

  for(int i=0; i<nentries;i++){
    if(signal_received){
      break;
    }
    tr->GetEvent(i);

    //Put in-loop actions here.

    if(i%1000 == 0){
      double time_end = get_time();
      cout<<setw(5)<<setiosflags(ios::fixed)<<setprecision(1)<<(100.*i)/nentries<<" % done\t"<<(Float_t)i/(time_end - time_start)<<" events/s " << (nentries-i)*(time_end - time_start)/(Float_t)i<<"s to go \r"<<flush;
    }
  }
  cout << endl;  

  delete tr;

  //Put after-loop actions here.

  double time_end = get_time();
  cout << "Program run time " << time_end - time_start << " s." << endl;
  
  return 0;
}
