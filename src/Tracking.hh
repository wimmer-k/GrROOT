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
#ifndef __TRACKING_HH
#define __TRACKING_HH

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "TEnv.h"
#include <algorithm>
#include "GretinaCalc.hh"
#include "GretinaTrack.hh"
#include "TrackSettings.hh"

using namespace std;

#define MAXNIPOINTS 8
#define NPERM MAXNIPOINTS*7*6*5*4*3*2*1

class Tracking {
public:
  Tracking(){
    fset = new TrackSettings();
    GeneratePermutations();
  };
  Tracking(TrackSettings* set){
    fset = set;
    GeneratePermutations();    
  };
  ~Tracking();
  void GeneratePermutations();
  void SetGretina(GretinaCalc *gr){fgr = gr;}
  void SortInClusters();
  vector<double> GetEsum(){return fesum;}
  GretinaEvent* GetEvent(){return fevt;}
  GretinaTrack* GetTrack(int n){return fevt->GetTrack(n);}
  void TrackCluster(GretinaTrack*);
private:
  vector<HitCalc*> SortInCluster(vector<HitCalc*> hits);
  double FigureOfMerit(vector<HitCalc*> hits, int nperm, double esum);
  
  TrackSettings* fset;
  GretinaCalc *fgr;
  GretinaEvent *fevt;
  vector<double> fesum;
  int fperm[MAXNIPOINTS+1][NPERM][MAXNIPOINTS+1];
  int fnperm[MAXNIPOINTS+1];
  int fcurrentsize;
};

#endif
