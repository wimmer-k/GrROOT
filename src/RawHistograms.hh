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
////////////                 Kathrin Wimmer  (wimme1k@cmich.edu)
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
#ifndef RAW_HISTOGRAMS_HH__
#define RAW_HISTOGRAMS_HH__

#include <iostream>
#include <iomanip>

#include "TObject.h"
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

class S800;
class Gretina;
class Mode3Event;

class RawHistograms {
public:
  RawHistograms(Settings* set = NULL,int nentries=100000){
    fentry = -1;
    fhlist = new TList;
    fnentries = nentries;
    fSett = set;
    for(int i=0; i<16; i++){
      fcounter[i] = 0;
    }
  }
  ~RawHistograms(){
    delete fhlist;
  }
  
  TList* GetHList(){return fhlist;}
  void Write();
  void FillHistograms(Gretina* gr, S800* s800, Mode3Event* m3r);

  void Fill(string name,int bins, double low, double high, double value){
    try{
      fhmap.at(name)->Fill(value);
    } catch(out_of_range e) {
      //cout << "New 1-d histogram: " << name << endl;
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
      //cout << "New 2-d histogram: " << name << endl;
      TH2F* newHist = new TH2F(name.c_str(),name.c_str(),
			       binsX,lowX,highX,
			       binsY,lowY,highY);
      newHist->Fill(valueX,valueY);
      fhlist->Add(newHist);
      fhmap[name] = newHist;
    }
  }

protected:
  void FillMode2Histograms(Gretina* gr);
  void FillMode3Histograms(Mode3Event* m3r);
  void FillS800Histograms(S800* s800);

  TList* fhlist;
  map<string,TH1*> fhmap;
  Settings* fSett;
  int fnentries;
  int fcounter[16];
  int fentry;
};

#endif
