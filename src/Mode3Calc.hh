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
#ifndef __MODE3CALC_HH
#define __MODE3CALC_HH

#include <iostream>
#include <vector>
#include <math.h>
#include "TObject.h"
#include "TVector3.h"
#include "TMath.h"

using namespace std;
class Mode3Detection : public TObject {
public:
  Mode3Detection(float en, long long int ts, int id){
    //cout << __PRETTY_FUNCTION__ << endl;
    fen = en;
    ftime = ts;
    fid = id;
    fLargestEnergy = en;
  }
  Mode3Detection(Mode3Detection* old){
    //cout << __PRETTY_FUNCTION__ << endl;
    fen = old->GetEn();
    ftime = old->GetTime();
    fid = old->GetID();
    fLargestEnergy = fen;
  }
  Mode3Detection(){
    fen = sqrt(-1.0);
    fLargestEnergy = fen;
    ftime = -1;
    fid = -1;
  }
  void Clear(){
    fen = sqrt(-1.0);
    fLargestEnergy = fen;
    ftime = -1;
    fid = -1;
  }
  
  float GetEn(){
    return fen;
  }
  long long int GetTime(){
    return ftime;
  }
  int GetID(){
    return fid;
  }
  int GetHole(){
    return fid/4;
  }
  int GetCrystal(){
    return fid%4;
  }
  void AddBackHit(Mode3Detection* NewHit){
    float newEnergy = NewHit->GetEn();
    if (newEnergy > fLargestEnergy){
      fid = NewHit->GetID();
      ftime = NewHit->GetTime();
      fLargestEnergy = newEnergy;
    }
    fen += newEnergy;
  }
protected:
  float fen;
  float fLargestEnergy;
  long long int ftime;
  int fid;

  ClassDef(Mode3Detection, 1);
};

class Mode3Calc : public TObject {
public:
  Mode3Calc(){
    fmult =0;
    fhits.clear();
    fmult_ab = 0;
    fhits_ab.clear();
  }
  void Clear(){
    fmult =0;
    fhits.clear();
    fmult_ab = 0;
    fhits_ab.clear();
  }
  void ClearAddBack(){
    fmult_ab = 0;
    fhits_ab.clear();
  }
  void SetMult(int mult){
    fmult = mult;
    fhits.resize(mult);
  }
  void AddHit(int n, Mode3Detection* hit){
    if(n>=fmult || n<0) 
      return;
    fhits[n] = hit;
  }
  void AddHit(Mode3Detection* hit){
    //cout << __PRETTY_FUNCTION__ << " fmult " << fmult << endl;
    fmult++;
    fhits.push_back(hit);
  }

  void SetMultAB(int mult){
    fmult_ab = mult;
    fhits_ab.resize(mult);
  }
  void AddHitAB(int n, Mode3Detection* hit){
    if(n>=fmult_ab || n<0) 
      return;
    fhits_ab[n] = hit;
  }
  void AddHitAB(Mode3Detection* hit){
    fmult_ab++;
    fhits_ab.push_back(hit);
  }

  Short_t GetMult(){
    return fmult;
  }
  vector<Mode3Detection*> GetHits(){
    return fhits;
  }
  Mode3Detection* GetHit(int n){
    if(n>=fmult || n<0) // check everything users are stupid....
      return NULL;
    return fhits[n];
  }
  Short_t GetMultAB(){
    return fmult_ab;
  }
  vector<Mode3Detection*> GetHitsAB(){
    return fhits_ab;
  }
  Mode3Detection* GetHitAB(int n){
    if(n>=fmult_ab || n<0) 
      return NULL;
    return fhits_ab[n];
  }


protected:
  Short_t fmult;
  vector<Mode3Detection*> fhits;
  Short_t fmult_ab;
  vector<Mode3Detection*> fhits_ab;
  
  ClassDef(Mode3Calc, 1);
};
  
  
#endif
