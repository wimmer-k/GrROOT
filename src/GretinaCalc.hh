#ifndef __GRETINACALC_HH
#define __GRETINACALC_HH

#include <iostream>
#include <vector>
#include <math.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <stdexcept>

#include "TObject.h"
#include "TVector3.h"
#include "TMath.h"

#include "Gretinadefs.h"
#include "Gretina.hh"
#include "Trace.hh"

using namespace std;

class TRACK;
class Settings;

//! A class to contain a single calibrated interaction.
/*!
  The HitCalc class is meant to contain a single calibrated interaction.
  This may be the result of multiple addback from multiple sources,
    or may be the result of a single interaction point.
 */
class HitCalc : public TObject {
public:
  //! Required default constructor.
  HitCalc(){
    Clear();
  }
  //! Manually makes a hit.
  HitCalc(Short_t hole, Short_t crystal, Float_t energy, long long int timestamp, TVector3 pos, Float_t t0, Float_t chisq,int index=-1){
    fholenum = hole;
    fcrystal = crystal;
    fen = energy;
    fMaxSingleHit = energy;
    fDCen = sqrt(-1);
    ftimestamp = timestamp;
    fposition = pos;
    ftime = t0;
    fchisq = chisq;
    fHitsAdded = 1;
    fIndex = index;
  }
  //! Creates a hit from the main interaction point of a mode2 crystal.
  HitCalc(Crystal* cry){
    Clear();
    fholenum = cry->GetHole();
    fcrystal = cry->GetCrystal();
    fen = fMaxSingleHit = cry->GetEnergy();
    fDCen = sqrt(-1);
    ftimestamp = cry->GetTS();
    if (cry->GetMaxIPNr()!=-1){
      fposition = cry->GetMaxIP()->GetPosition();
      fTrueFirst = cry->GetIPoint(0)->GetPosition();
      fUsedTrueFirst = (cry->GetMaxIPNr() == 0);
    }
    ftime = cry->GetT0();
    fchisq = cry->GetChiSq();
    fHitsAdded = 1;
  }
  //! Copies an existing hit.
  HitCalc(HitCalc* hit){
    Clear();
    fholenum = hit->GetHole();
    fcrystal = hit->GetCrystal();
    fen = hit->GetEnergy();
    fDCen = hit->GetDCEnergy();
    fMaxSingleHit = hit->GetEnergy();
    fposition = hit->GetPosition();
    ftimestamp = hit->GetTS();
    fchisq = hit->GetChiSq();
    ftime = hit->GetTime();
    fHitsAdded = 1;
    fIndex = hit->GetIndex();
  }
  ~HitCalc(){
    Clear();
  }
  void Clear(){
    fholenum = -1;
    fcrystal = -1;
    ftime = fen = fDCen = fMaxSingleHit = sqrt(-1.0);
    ftimestamp = -1;
    fHitsAdded = 0;
    fchisq = 0;
    fposition.SetXYZ(0,0,0);
    fIndex = -1;
    //Simulation parameters
    fTrueFirst.SetXYZ(0,0,0);
    fUsedTrueFirst = false;
  }
  //! Add another HitCalc to the current HitCalc
  /*!
    Add another HitCalc into this one.
    The timestamp and position are taken from the hit with more energy.
   */
  void AddBackHitCalc(HitCalc* other){
    if (other->GetEnergy() > fMaxSingleHit || isnan(fMaxSingleHit)){
      fholenum = other->GetHole();
      fcrystal = other->GetCrystal();
      ftimestamp = other->GetTS();
      fMaxSingleHit = other->GetEnergy();
      fposition = other->GetPosition();
      ftime = other->GetTime();
      fchisq = other->GetChiSq();
    }
    fen += other->GetEnergy();
    fHitsAdded += other->GetHitsAdded();
  }
  void SetDCEnergy(float dcen){fDCen = dcen;}
  void SetPosition(TVector3 in){fposition = in;}
  void SetChiSq(float chisq){fchisq = chisq;}
  //! An integer identifying the crystal, equal to 4*Hole+(Position-1)
  Short_t GetID(){return fholenum*4 + fcrystal;}
  //! The number of the hole, ranging from 1-30.
  Short_t GetHole(){return fholenum;}
  //! The number of the crystal, ranging from 0-3.
  Short_t GetCrystal(){return fcrystal;}
  //! The energy of the hit (keV).
  Float_t GetEnergy(){return fen;}
  //! The Doppler-corrected energy of the hit.
  Float_t GetDCEnergy(){return fDCen;}
  Float_t GetDCEnergy(float beta, double x=0, double y=0, double z=0){
    TVector3 pos;
    pos.SetXYZ(x,y,z);//in mm
    pos += GetPosition();
    return fen/sqrt(1-beta*beta)*(1-beta*cos(pos.Theta()));
  }
  Float_t GetDCEnergy(float beta, TVector3 gpos){
    return fen/sqrt(1-beta*beta)*(1-beta*cos(gpos.Theta()));
  }
  Float_t GetDCEnergy(float beta, TVector3 gpos, TVector3  s800dir){
    return fen/sqrt(1-beta*beta)*(1-beta*cos(gpos.Angle(s800dir)));
  }
  long long int GetTS(){return ftimestamp;}
  Float_t GetChiSq(){return fchisq;}
  Float_t GetMaxSingleHit(){return fMaxSingleHit;}
  //! The position of the hit in the lab system.
  TVector3 GetPosition(){return fposition;}
  int GetHitsAdded(){return fHitsAdded;}
  Float_t GetTime(){return ftime;}

  void CorrectTime(Trace* tr){
    ftime += ftimestamp - tr->GetLED();
  }
  void DopplerCorrect(double beta, double z = 0){
    fDCen = fen * HitCalc::DopplerCorrectionFactor(GetPosition(),beta,z);
  }
  //! Apply the Doppler correction using the given settings and S800 data.
  /*!
    Apply the Doppler correction using the given settings and S800 data.
    Uses the beta and the target position from the settings file.
    Uses the Phi, Theta, YTA, and DTA from the given TRACK.
   */
  void DopplerCorrect(Settings* set, TRACK* track);
  //! Returns the Doppler-correction factor to correct the energy.
  static double DopplerCorrectionFactor(TVector3 PosToTarget, Settings* set, TRACK* track);
  static double DopplerCorrectionFactor(TVector3 PosToTarget, double beta, double z);
  void Print(){
    cout << "hole " << fholenum << "\tcrystal " << fcrystal << "\ten " << fen << "\tmax hit " << fMaxSingleHit << endl;//"\tipoints " << fipoints.size()<< endl;
    return;
  }

  int GetIndex(){return fIndex;}
  void SetIndex(int i){fIndex = i;}

  //Extra functions for simulated files.
  TVector3 GetTrueFirstPos(){return fTrueFirst;}
  bool GetUsedTrueFirstIP(){return fUsedTrueFirst;}
  double GetThetaFromTrue(){return (fposition.Angle(fTrueFirst));}
  double GetThetaDiff(){return fposition.Theta() - fTrueFirst.Theta();}
  double GetDistanceFromTrue(){return (fposition-fTrueFirst).Mag();}
  double GetPerpDistanceFromTrue(){
    TVector3 diff = fposition - fTrueFirst;
    diff -= diff.Dot(fposition.Unit())*fposition.Unit();
    return diff.Mag();
  }
  double GetDCEnergySimCheat(){return fDCen_simcheat;}
  TVector3 GetTruePos(){return fTrueFirst;}
  TVector3 GetMainPos(){return fposition;}

protected:
  Short_t fholenum;
  Short_t fcrystal;
  Float_t fchisq;
  Float_t fen;
  Float_t fDCen;
  Float_t ftime;
  Float_t fMaxSingleHit; //!
  long long int ftimestamp;
  TVector3 fposition;
  int fHitsAdded;
  int fIndex;
  //Extra variables for simulations
  TVector3 fTrueFirst;
  bool fUsedTrueFirst;
  Float_t fDCen_simcheat;
  ClassDef(HitCalc, 1);
};


//! The calibrated GRETINA object.
/*!
  Holds all calibrated information of GRETINA.
  Holds three lists, keeping the crystal-by-crystal information,
    the add-backed information,
    and the interaction point cluster information.
 */
class GretinaCalc : public TObject {
public:
  GretinaCalc(){
    Clear();
  }
  void Clear(){
    fmult = 0;
    for(vector<HitCalc*>::iterator cry=fhits.begin(); cry!=fhits.end(); cry++){
      delete *cry;
    }
    fhits.clear();
    ClearAddBack();
    ClearCluster();
  }
  void ClearAddBack(){
    fmult_ab = 0;
    for(vector<HitCalc*>::iterator cry=fhits_ab.begin(); cry!=fhits_ab.end(); cry++){
      delete *cry;
    }
    fhits_ab.clear();
  }
  void ClearCluster(){
    for(UShort_t i=0;i<fmult_cl.size();i++){
      for(vector<HitCalc*>::iterator cry=fhits_cl[i].begin(); cry!=fhits_cl[i].end(); cry++){
	delete *cry;
      }
      fhits_cl[i].clear();
    }
    fhits_cl.clear();
    fmult_cl.clear();
    fncluster =0;
  }
  void AddHit(HitCalc* cry){
    fhits.push_back(cry);
    fmult++;
  }
  void AddHitAB(HitCalc* cry){
    fmult_ab++;
    fhits_ab.push_back(cry);
  }
  void AddHitCL(vector<HitCalc*> cry, UShort_t c){
    if(c<fncluster){
      //cout << "invalid " << endl;
      return;
    }
    fncluster++;
    fmult_cl.push_back(cry.size());
    fhits_cl.push_back(cry);
  }
  void AddHitCL(HitCalc* cry, UShort_t c){
    if(c>=fncluster){
      //cout << "invalid " << endl;
      return;
    }
    fmult_cl[c]++;
    fhits_cl[c].push_back(cry);
  }
  void DopplerCorrect(double beta, double z = 0){
    for(vector<HitCalc*>::iterator hit=fhits.begin(); hit!=fhits.end(); hit++){
      (*hit)->DopplerCorrect(beta, z);
    }
    for(vector<HitCalc*>::iterator hit=fhits_ab.begin(); hit!=fhits_ab.end(); hit++){
      (*hit)->DopplerCorrect(beta, z);
    }
  }
  void DopplerCorrect(Settings* set, TRACK* track);
  void CorrectTime(Trace* tr){
    for(vector<HitCalc*>::iterator hit=fhits.begin(); hit!=fhits.end(); hit++){
      (*hit)->CorrectTime(tr);
    }
    for(vector<HitCalc*>::iterator hit=fhits_ab.begin(); hit!=fhits_ab.end(); hit++){
      (*hit)->CorrectTime(tr);
    }
    if(fncluster>0){
       for(UShort_t i=0;i<fncluster;i++){
	for(vector<HitCalc*>::iterator hit=fhits_cl[i].begin(); hit!=fhits_cl[i].end(); hit++){
	  (*hit)->CorrectTime(tr);
	}
      }
    }
  }
  void Print(){
    cout << " singles mult " <<fmult << endl;
    for(vector<HitCalc*>::iterator hit=fhits.begin(); hit!=fhits.end(); hit++){
      (*hit)->Print();
    }
    if(fmult_ab>0){
      cout << "addback mult " <<fmult_ab <<endl;
      for(vector<HitCalc*>::iterator hit=fhits_ab.begin(); hit!=fhits_ab.end(); hit++){
	(*hit)->Print();
      }
    }
    if(fncluster>0){
      cout << fncluster << " cluster created" << endl;
      for(UShort_t i=0;i<fncluster;i++){
	cout << i << " mult " << fmult_cl[i] << endl;
	for(vector<HitCalc*>::iterator hit=fhits_cl[i].begin(); hit!=fhits_cl[i].end(); hit++){
	  (*hit)->Print();
	}
      }
    }
    cout << "-----------------------------------"<<endl;
  }
  void PrintCluster(int i){
    cout << i << " mult " << fmult_cl[i] << endl;
    for(vector<HitCalc*>::iterator hit=fhits_cl[i].begin(); hit!=fhits_cl[i].end(); hit++){
      (*hit)->Print();
    }
  }

  int GetMult(){return fmult;}
  vector<HitCalc*> GetHits(){return fhits;}
  HitCalc* GetHit(int n){return fhits[n];}
  vector<HitCalc*> GetHitsAB(){return fhits_ab;}
  int GetMultAB(){return fmult_ab;}
  HitCalc* GetHitAB(int n){return fhits_ab[n];}

  int GetNCluster(){return fncluster;}
  vector<Short_t> GetMultCL(){return fmult_cl;}
  int GetMultCL(int n){return fmult_cl[n];}
  vector<vector<HitCalc*> > GetClusters(){return fhits_cl;}
  vector<HitCalc*> GetClusters(int n){return fhits_cl[n];}
  HitCalc* GetHitCL(int n, int m){return fhits_cl[n][m];}

  void SetCluster(int n, vector<HitCalc*> hits){
    fhits_cl[n] = hits;
  }

protected:
  Short_t fmult;
  vector<HitCalc*> fhits;
  Short_t fmult_ab;
  vector<HitCalc*> fhits_ab;
  Short_t fncluster;
  vector<Short_t> fmult_cl;
  vector<vector<HitCalc*> > fhits_cl;
  ClassDef(GretinaCalc, 1);
};

//! Allows for comparing of HitCalc objects.
/*!
  Allows for ordering of HitCalc objects by decreasing energy.
  Used as a first guess of an path for tracking.
 */
class HitComparer {
public:
  bool operator() ( HitCalc *lhs, HitCalc *rhs) {
    return (*lhs).GetEnergy() > (*rhs).GetEnergy();
  }
};
#endif

