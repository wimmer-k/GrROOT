#include "Gretina.hh"
#include "S800.hh"
#include "Trace.hh"
#include "TFile.h"
#include "TTree.h"

void ViewAround(int n){
  TFile* file = TFile::Open("/global/data1x/gretina/e11007_part_I/Run0131/raw.root");
  TTree* tree = (TTree*)file->Get("gtr");

  Gretina* gr = new Gretina;
  S800* s800 = new S800;
  Mode3Event* mode3 = new Mode3Event;

  
  tree->SetBranchAddress("gretina",&gr);
  tree->SetBranchAddress("s800",&s800);
  tree->SetBranchAddress("mode3Event",&mode3);

  for(int i=n-2; i<=n+2; i++){
    tree->GetEvent(i);
    cout << Form("Event Number: %d",i) << endl;
    cout << "S800:   " << s800->GetTS() << ", " << s800->GetInternalTS() << endl;
    for(int cry=0; cry<gr->GetMult(); cry++){
      cout << "Mode2:  " << gr->GetHit(cry)->GetTS() << ", " << gr->GetHit(cry)->GetITS() << endl;
    }
    for(int hit=0; hit<mode3->GetMult(); hit++){
      for(int tr=0; tr<mode3->GetHit(hit)->GetMult(); tr++){
	Trace* trace = mode3->GetHit(hit)->GetTrace(tr);
	cout << "Card29: " << trace->GetTS() << ", " << trace->GetLED() << endl;
      }
    }
    cout << endl;  
  }
}
