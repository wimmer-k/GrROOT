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
#ifndef __GRETINASIM_HH
#define __GRETINASIM_HH
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

#include "Simdefs.h"
using namespace std;

typedef struct g4sim_emitted_gamma{
  float e;
  float x, y, z;
  float phi, theta;
  float beta;
} EG;

typedef struct g4sim_abcd1234 {
  int type;          /* defined as abcd1234 */
  int num;           /* # of emitted gammas */
  int full;          /* is full energy */
  EG gammas[MAX_SIM_GAMMAS];
} G4SIM_EGS;

// GEB payload for simulated S800 events
typedef struct S800_physicsdata {
  int type;    /* defined abcd1234 for indicating this version */
  float crdc1_x;   /* Crdc x/y positions in mm */
  float crdc1_y;
  float crdc2_x;
  float crdc2_y;
  float ic_sum;    /* ion chamber energy loss         */
  float tof_xfp;   /* TOF scintillator after A1900    */
  float tof_obj;   /* TOF scintillator in object box  */
  float rf;        /* Cyclotron RF for TOF            */ 
  int trigger; /* Trigger register bit pattern    */
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
  /* from here corrected values extracted from data above */ 
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
  float ic_de;
  /* TOF values with TOF correction applied (from afp/crdc x) */
  float tof_xfpe1;
  float tof_obje1;
  float tof_rfe1;
  /* Trajectory information at target position calculated from 
     a map and afp/bfp/xfp/yfp. New map and you need to re-calc */
  float ata; /* dispersive angle        */
  float bta; /* non-dispersive angle    */
  float dta; /* dT/T T:kinetic energy   */
  float yta; /* non-dispersive position */
} S800_PHYSICSDATA;

class EmittedGamma : public TObject {
public:
  EmittedGamma(){
    Clear();
  }
  EmittedGamma(EG* eg) {
    fenergy = eg->e;
    femissionX = eg->x;
    femissionY = eg->y;
    femissionZ = eg->z;
    femissionPhi = eg->phi;
    femissionTheta = eg->theta;
    femissionBeta = eg->beta;
  }
  EmittedGamma(Float_t e,
	       Float_t x, Float_t y, Float_t z,
	       Float_t ph, Float_t th, Float_t b){
    fenergy = e;
    femissionX = x;
    femissionY = y;
    femissionZ = z;
    femissionPhi = ph;
    femissionTheta = th;
    femissionBeta = b;
  }
  ~EmittedGamma(){
  }
  void Clear(){
    fenergy = 0;
    femissionX = 0;
    femissionY = 0;
    femissionZ = 0;
    femissionPhi = 0;
    femissionTheta = 0;
    femissionBeta = 0;
  }
  Float_t GetEnergy(){ return fenergy; }
  TVector3 GetPos(){
    return TVector3(femissionX, femissionY, femissionZ);
  }
  Float_t GetX(){ return femissionX; }
  Float_t GetY(){ return femissionY; }
  Float_t GetZ(){ return femissionZ; }
  Float_t GetPhi(){ return femissionPhi; }
  Float_t GetTheta(){ return femissionTheta; }
  Float_t GetBeta(){ return femissionBeta; }

  void SetEnergy(Float_t e) { fenergy = e; }
  void SetX(Float_t x) { femissionX = x; }
  void SetY(Float_t y) { femissionY = y; }
  void SetZ(Float_t z) { femissionZ = z; }
  void SetPos(Float_t x, Float_t y, Float_t z) {
    femissionX = x;
    femissionY = y;
    femissionZ = z;
  }
  void SetPhi(Float_t ph) { femissionPhi = ph; }
  void SetTheta(Float_t th) { femissionTheta = th; }
  void SetBeta(Float_t b) { femissionBeta = b; }

protected:
  Float_t fenergy;
  Float_t femissionX;
  Float_t femissionY;
  Float_t femissionZ;
  Float_t femissionPhi;
  Float_t femissionTheta;
  Float_t femissionBeta;

  ClassDef(EmittedGamma, 1);
};

class GretinaSim : public TObject {
public:
  GretinaSim(){
    Clear();
  }
  ~GretinaSim(){
    Clear();
  }
  void Clear(){
    ftimestamp = 0;
    fmult = 0;
    for(vector<EmittedGamma*>::iterator eg=fgammas.begin(); eg!=fgammas.end(); eg++){
      delete *eg;
    }
    fgammas.clear();
  }
  EmittedGamma* GetEmittedGamma(Int_t i){return fgammas[i];}
  Int_t GetMult(){return fmult;}
  long long int GetTimeStamp(){return ftimestamp;}

  void SetTimeStamp(long long int ts){ ftimestamp = ts; }
  void AddEmittedGamma(Float_t e, Float_t x, Float_t y, Float_t z,
		       Float_t ph, Float_t th, Float_t b){
    fmult++;
    fgammas.push_back( new EmittedGamma(e, x, y, z, ph, th, b) );
  }
  void AddEmittedGamma(EG* eg){
    fmult++;
    fgammas.push_back( new EmittedGamma(eg->e,
					eg->x, eg->y, eg->z,
					eg->phi, eg->theta, eg->beta) );
  }

protected:
  long long int ftimestamp;
  Short_t fmult;
  vector<EmittedGamma*> fgammas;
  ClassDef(GretinaSim, 1);
};

#endif
