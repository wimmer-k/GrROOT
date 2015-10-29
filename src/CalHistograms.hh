#ifndef CAL_HISTOGRAMS_HH__
#define CAL_HISTOGRAMS_HH__

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

#define CUTFILE_NAMELENGTH 100

class GretinaCalc;
class S800Calc;
class Mode3Calc;

class CalHistograms {
public:
  CalHistograms(Settings* set = NULL,int tac=0, int cal = 0, char* cutfile=NULL, int nentries=100000){
    Init(set,tac,cal,cutfile,nentries);
  }
  ~CalHistograms(){
    delete fhlist;
  }
  void Init(Settings* set,int tac, int cal, char* cutfile,int nentries){
    if (set!=NULL){
      fSett = set;
    } else {
      fSett = new Settings;
    }
    fCal = cal;
    fhlist = new TList;
    fnentries = nentries;
    if(tac > -1 && tac < 4 ){
      ftac = tac;
    } else {
      cout << "Unknown timing number: " << tac << endl;
      cout << "Please use 0 for TAC timing or 1 for TDC timing if you want to use the OBJ scintillator" << endl;
      cout << "Please use 2 for TAC timing or 3 for TDC timing if you want to use the XFP scintillator" << endl;
      cout << "Assuming that the cuts were made for TAC timing and the OBJ scintillator" << endl;
      ftac = 0;
    }
    if(cutfile!=NULL){
      cout << "using cutfile " << cutfile << endl;
      if(strlen(cutfile)<CUTFILE_NAMELENGTH-1){
	strcpy(fcutfile,cutfile);
	fhasfile = true;
      } else {
	cout << "cutfile filename too long" << endl;
	cout << "please increase CUTFILE_NAMELENGTH in CalHistograms.hh" << endl;
	fhasfile = false;
      }
    } else {
      fhasfile = false;
    }
    fentry =0;

    for(int i=0;i<3;i++){
      fobj_range[i] = fSett->GetRangeOBJ(i);
      fxfp_range[i] = fSett->GetRangeXFP(i);
      ftacobj_range[i] = fSett->GetRangeTACOBJ(i);
      ftacxfp_range[i] = fSett->GetRangeTACXFP(i);
      fIC_range[i] = fSett->GetRangeIC(i);
      fobjC_range[i] = fSett->GetRangeOBJC(i);
      fxfpC_range[i] = fSett->GetRangeXFPC(i);
      ftacobjC_range[i] = fSett->GetRangeTACOBJC(i);
      ftacxfpC_range[i] = fSett->GetRangeTACXFPC(i);
      fPP_range[i] = fSett->GetRangePP(i);
    }

  }
  void AddEntry(){
    fentry++;
  }
  TList* GetHList(){return fhlist;}
  void Write(){
    fhlist->Sort();
    fhlist->Write();
  }

  void FillHistograms(GretinaCalc* gr, S800Calc* s800, Mode3Calc* m3c);

  void FillI(string name,int bins, double low, double high, double value){
    try{
      fhmap1d.at(name)->Fill(value);
    } catch(out_of_range e) {
      TH1I* newHist = new TH1I(name.c_str(),name.c_str(),bins,low,high);
      newHist->Fill(value);
      fhlist->Add(newHist);
      fhmap1d[name] = newHist;
    }
  }
  void FillS(string name,int bins, double low, double high, double value){
    try{
      fhmap1d.at(name)->Fill(value);
    } catch(out_of_range e) {
      TH1S* newHist = new TH1S(name.c_str(),name.c_str(),bins,low,high);
      newHist->Fill(value);
      fhlist->Add(newHist);
      fhmap1d[name] = newHist;
    }
  }
  void Fill(string name,int bins, double low, double high, double value, double weight =1){
    try{
      fhmap1d.at(name)->Fill(value,weight);
    } catch(out_of_range e) {
      TH1F* newHist = new TH1F(name.c_str(),name.c_str(),bins,low,high);
      newHist->Fill(value,weight);
      fhlist->Add(newHist);
      fhmap1d[name] = newHist;
    }
  }
  void Fill(string name,
	    int binsX, double lowX, double highX, double valueX,
	    int binsY, double lowY, double highY, double valueY, 
	    double weight =1){
    try{
      fhmap2d.at(name)->Fill(valueX,valueY,weight);
    } catch(out_of_range e) {
      TH2F* newHist = new TH2F(name.c_str(),name.c_str(),
			       binsX,lowX,highX,
			       binsY,lowY,highY);
      newHist->Fill(valueX,valueY,weight);
      fhlist->Add(newHist);
      fhmap2d[name] = newHist;
    }
  }
  void FillHistogramsNoGate(GretinaCalc* gr, S800Calc* s800, Mode3Calc* m3c);
  void FillHistogramsGateIn(GretinaCalc* gr, S800Calc* s800, Mode3Calc* m3c,const char* inname);
  void FillHistogramsGateOut(GretinaCalc* gr, S800Calc* s800, Mode3Calc* m3c,const char* outname);
  void FillHistogramsGateOnlyOut(GretinaCalc* gr, S800Calc* s800, Mode3Calc* m3c,const char* outname);

protected:
  TList* fhlist;
  map<string,TH1*> fhmap1d;
  map<string,TH2*> fhmap2d;
  int fnentries;
  int fentry;
  char fcutfile[CUTFILE_NAMELENGTH];
  bool fhasfile;
  int ftac;
  int fCal;
  Settings* fSett;

  //histo ranges
  int fobj_range[3];
  int fxfp_range[3];
  int ftacobj_range[3];
  int ftacxfp_range[3];
  int fIC_range[3];
  int fobjC_range[3];
  int fxfpC_range[3];
  int ftacobjC_range[3];
  int ftacxfpC_range[3];
  double fPP_range[3];

};

#undef CUTFILE_NAMELENGTH

#endif
