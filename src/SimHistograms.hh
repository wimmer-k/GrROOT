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
#ifndef SIM_HISTOGRAMS_HH__
#define SIM_HISTOGRAMS_HH__

#include <iostream>
#include <iomanip>

#include "TObject.h"
#include "TFile.h"
#include "TTree.h"
#include "TList.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TEnv.h"
#include "math.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdexcept>

#include "Settings.hh"

using namespace std;

class GretinaCalc;
class S800Calc;
class Mode3Calc;
class GretinaSim;

class SimHistograms {
public:
  SimHistograms(Settings* set = NULL, int nentries=100000){
    Init(set, nentries);
  }
  ~SimHistograms(){
    delete fhlist;
  }
  void Init(Settings* set, int nentries){
    if (set!=NULL){
      fSett = set;
    } else {
      fSett = new Settings;
    }
    fhlist = new TList;
    fnentries = nentries;
    fentry =0;
  }
  void AddEntry(){
    fentry++;
  }
  TList* GetHList(){return fhlist;}
  void Write(){
    fhlist->Sort();
    fhlist->Write();
  }

  void FillHistograms(GretinaCalc* gr, S800Calc* s800, Mode3Calc* m3c,
		      GretinaSim* gs);

  void Fill(string name,int bins, double low, double high, double value){
    try{
      fhmap.at(name)->Fill(value);
    } catch(out_of_range e) {
      cout << "New 1-d histogram: " << name << endl;
      TH1F* newHist = new TH1F(name.c_str(),name.c_str(),bins,low,high);
      newHist->Fill(value);
      fhlist->Add(newHist);
      fhmap[name] = newHist;
    }
  }
  void Fill(string name,
	    int binsX, double lowX, double highX, double valueX,
	    int binsY, double lowY, double highY, double valueY){
    try{
      fhmap.at(name)->Fill(valueX,valueY);
    } catch(out_of_range e) {
      cout << "New 2-d histogram: " << name << endl;
      TH2F* newHist = new TH2F(name.c_str(),name.c_str(),
			       binsX,lowX,highX,
			       binsY,lowY,highY);
      newHist->Fill(valueX,valueY);
      fhlist->Add(newHist);
      fhmap[name] = newHist;
    }
  }

protected:
  TList* fhlist;
  map<string,TH1*> fhmap;
  int fnentries;
  int fentry;
  Settings* fSett;

};

#endif
