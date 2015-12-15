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
#ifndef __SCALER_HH
#define __SCALER_HH

#include <iostream>
#include "TObject.h"
#include "Scalerdefs.h"
using namespace std;
class Scaler : public TObject {
public:
  Scaler(){
    Clear();
  }
  void Clear(){
    for(int i=0;i<NSCALER;i++){
      fvalues[i] =-1;
      foverflows[i] =0;
    }
    fstartTime =-1;
    fendTime =-1;
    fts =-1;  
    feventnr = -1;
  }
  void SetOverflows(int i, int n){foverflows[i] = n;}
  void SetValue(int i, ULong64_t val){fvalues[i] = val;}
  void SetStart(ULong64_t t){fstartTime = t;}
  void SetEnd(ULong64_t t){fendTime = t;}
  void SetTS(ULong64_t ts){fts = ts;}
  void SetInternalTS(ULong64_t its){fits = its;}
  void SetEventNr(double evtnr){feventnr = evtnr;}

  ULong64_t GetValue(int i){return fvalues[i];}
  ULong64_t GetStart(){return fstartTime;}
  ULong64_t GetEnd(){return fendTime;}
  ULong64_t GetTS(){return fts;}
  ULong64_t GetInternalTS(){return fits;}
  double GetEventNr(){return feventnr;}
  void Print(){
    cout << "Start Time: " << fstartTime << endl;
    cout << "End Time: " << fendTime << endl;
    cout << "Timestamp: " << fts << endl;
    cout << "Internal Timestamp: " << fits << endl;
    cout << "Event Number: " << feventnr << endl;
    for(int i=0; i<NSCALER; i++){
      cout << "Scaler " << i << ": " << fvalues[i] << endl;
    }
  }

protected:
  ULong64_t fvalues[NSCALER];
  int foverflows[NSCALER];
  ULong64_t fstartTime;                  //  Interval start time in seconds.
  ULong64_t fendTime;                    //  Interval end time in seconds.
  ULong64_t fts;
  ULong64_t fits;
  double feventnr;
  ClassDef(Scaler, 1);
};

#endif
