#ifndef __GRETINATRACK_HH
#define __GRETINATRACK_HH

#include <iostream>
#include <vector>
#include <math.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <stdlib.h>

#include "TObject.h"
#include "TVector3.h"
#include "TMath.h"

#include "GretinaCalc.hh"

using namespace std;

//! Class to contain multiple interaction points, used for tracking.
/*!
  Hold an ordered list of HitCalc objects,
    which can be ordered using Tracking::TrackCluster().
  Holds information on the energy of the cluster, and the figure of merit from the tracking.
 */
class GretinaTrack : public TObject {
public:
  GretinaTrack(){
    Clear();
  }
  GretinaTrack(vector<HitCalc*> hits){
    Clear();
    SetHits(hits);
  }
  GretinaTrack(HitCalc* hit){
    Clear();
    AddHit(hit);
  }
  void Clear(){
    fmult = 0;
    fFOM = -2;
    fperm = -1;
    fjumps = 0;
    fesum = -1;
    fDCen = -1;
    for(vector<HitCalc*>::iterator cry=fhits.begin(); cry!=fhits.end(); cry++){
      delete *cry;
    }
    fhits.clear();
  }
  GretinaTrack* Copy(){
    GretinaTrack* output = new GretinaTrack;
    output->fperm = fperm;
    output->fjumps = fjumps;
    output->fesum = fesum;
    output->fDCen = fDCen;
    for(vector<HitCalc*>::iterator hit=fhits.begin(); hit!=fhits.end(); hit++){
      output->AddHit(new HitCalc(*hit));
    }
    return output;
  }
  void SetHits(vector<HitCalc*> hits){
    fmult = hits.size();
    fhits = hits;
  }
  void AddHit(HitCalc* hit){
    fmult++;
    fhits.push_back(hit);
  }
  void DopplerCorrect(double beta, double z=0){
    fDCen = fesum * HitCalc::DopplerCorrectionFactor(fhits[0]->GetPosition(),beta,z);
  }
  void DopplerCorrect(Settings* set, TRACK* track){
    fDCen = fesum * HitCalc::DopplerCorrectionFactor(fhits[0]->GetPosition(),set,track);
  }
  void SetFOM(Double_t fom){ fFOM = fom; }
  void SetPermutation(Int_t perm){ fperm = perm; }
  void SetPermutation(vector<int> perm){ fpermlist = perm; }
  void SetJumps(Double_t jumps){ fjumps = jumps; }
  void SetEsum(Double_t esum){ fesum = esum; }
  void SetDCEnergy(Double_t DCen){ fDCen = DCen; }
  int GetMult(){return fmult;}
  Short_t GetJumps(){return fjumps;}
  Int_t GetPermutation(){return fperm;}
  Double_t GetFOM(){return fFOM;}
  Double_t GetEsum(){return fesum;}
  Double_t GetDCEnergy(){return fDCen;}
  vector<HitCalc*> GetHits(){return fhits;}
  HitCalc* GetHit(int n){return fhits[n];}
  double GetFirstScattering(){
    if (fhits.size()>=2){
      TVector3 pos1 = fhits[0]->GetPosition();
      TVector3 pos2 = fhits[1]->GetPosition();
      return pos1.Angle(pos2-pos1);
    } else {
      return sqrt(-1);
    }
  }
  double GetFirstEnergy(){
    if (fhits.size()>=1){
      return fhits[0]->GetEnergy();
    } else {
      return sqrt(-1);
    }
  }
  double GetFirstEnergyScattering(double initial){
    double en = GetFirstEnergy();
    double costheta = 1 - 511*(1/initial + 1/(initial-en));
    return acos(costheta);
  }
  void WriteIPoints(char* filename){
    FILE* file = fopen(filename,"w");
    for(uint i = 0; i<fhits.size(); i++){
      TVector3 pos =  fhits[i]->GetPosition();
      fprintf(file,
	      "%f\t%f\t%f\t%f\t%d\n",
	      pos.X(),pos.Y(),pos.Z(),
	      fhits[i]->GetEnergy(),
	      fhits[i]->GetIndex());
    }
    fclose(file);
  }
  bool IsTrackedMain(){
    double trackedEn = GetHit(0)->GetEnergy();
    for(int i=1; i<GetMult(); i++){
      if (GetHit(i)->GetEnergy() > trackedEn){
	return false;
      }
    }
    return true;
  }
  void Print(){
    cout << "mult " << fmult << " esum " << fesum << " FOM " << fFOM << endl;
    for(vector<HitCalc*>::iterator hit=fhits.begin(); hit!=fhits.end(); hit++){
      (*hit)->Print();
    }
  }
protected:
  Short_t fmult;
  Short_t fjumps;
  Int_t fperm;
  Double_t fFOM;
  Double_t fesum;
  Double_t fDCen;
  vector<HitCalc*> fhits;
  vector<int> fpermlist;
  ClassDef(GretinaTrack, 1);
};

class GretinaEvent : public TObject {
public:
  GretinaEvent(){
    Clear();
  }
  void Clear(){
    fmult = 0;
    fmult_ab = 0;
    fmult_raw = 0;
    for(vector<GretinaTrack*>::iterator track=ftracks.begin(); track!=ftracks.end(); track++){
      delete *track;
    }
    ftracks.clear();
    for(vector<HitCalc*>::iterator cry=fhits_ab.begin(); cry!=fhits_ab.end(); cry++){
      delete *cry;
    }
    fhits_ab.clear();
    for(vector<HitCalc*>::iterator cry=fhits_raw.begin(); cry!=fhits_raw.end(); cry++){
      delete *cry;
    }
    fhits_raw.clear();
  }
  void AddTrack(GretinaTrack* track){
    ftracks.push_back(track);
    fmult++;
  }
  void SetHitsRaw(vector<HitCalc*> hits){
    fmult_raw = hits.size();
    fhits_raw = hits;
  }
  void SetHitsAB(vector<HitCalc*> hits){
    fmult_ab = hits.size();
    fhits_ab = hits;
  }

  vector<GretinaTrack*> GetTracks(){return ftracks;}
  GretinaTrack* GetTrack(int n){return ftracks[n];}
  Short_t GetMult(){return fmult;}
  Short_t GetMultAB(){return fmult_ab;}
  Short_t GetMultRaw(){return fmult_raw;}
  void PrintEvent(){
    cout << " mult " << fmult << endl;
    for(vector<GretinaTrack*>::iterator track=ftracks.begin(); track!=ftracks.end(); track++){
      (*track)->Print();
    }
  }
  vector<HitCalc*> GetHitsAB(){return fhits_ab;}
  vector<HitCalc*> GetHitsRaw(){return fhits_raw;}
protected:
  Short_t fmult;
  vector<GretinaTrack*> ftracks;
  Short_t fmult_ab;
  Short_t fmult_raw;
  vector<HitCalc*> fhits_ab;
  vector<HitCalc*> fhits_raw;
  ClassDef(GretinaEvent, 1);
};

#endif
