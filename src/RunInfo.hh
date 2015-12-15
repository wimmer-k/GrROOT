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
#ifndef __RUNINFO_HH
#define __RUNINFO_HH

#include <iostream>
#include <iomanip>
#include "TObject.h"
#include "TMath.h"
#include "Scalerdefs.h"

using namespace std;

class RunInfo : public TObject {
public:
  RunInfo(){
    Clear();
  }
  ~RunInfo(){
    Clear();
  }
  void Clear(){
    frunnr = -1;
    fruntime = -1;
    for(int i=0;i<NSCALER;i++){
      fintegral[i] =0;
    }
    fentries = -1;
    fevents = -1;
    fbuffer = -1;
    fsbuffer = -1;
    fbrho = -1;
    fmass = -1;
    fcharge = -1;
    fich = -1;
    fhodo = -1;
    fgreta = -1;
    fcard29 = -1;
    fobjeff = -1;
    fxfpeff = -1;
    ftofeff = -1;
    fpadeff[0] = -1;
    fpadeff[1] = -1;
    ftrackeff = -1;
    fppaceff[0] = -1;
    fppaceff[1] = -1;
    fiitrackeff = -1;
    fcard29eff = -1;
  }
  void SetRunNr(int i){
    frunnr = i;
  }
  void SetRunTime(double t){
    fruntime = t;
  }
  void Integral(int i, ULong64_t val){
    fintegral[i] = val;
  }
  void SetCounters(Long64_t ic, Long64_t hodo, Long64_t card29, Long64_t greta, Long64_t scint=0){
    fich = ic;
    fscint = scint;
    if(scint<1){//for backwards compatability
      fscint = ic;
    }
    fhodo = hodo;
    fgreta = greta;
    fcard29 = card29;
  }
  void SetEff(Long64_t ic, Long64_t obj, Long64_t xfp, Long64_t tof, Long64_t pad[2], Long64_t track, Long64_t ppac[2], Long64_t iitrack, Long64_t card29, Long64_t scint=0){
    fich = ic;
    fobjeff = (float)obj/ic;
    fxfpeff = (float)xfp/ic;
    ftofeff = (float)tof/ic;
    fppaceff[0] = (float)ppac[0]/ic;
    fppaceff[1] = (float)ppac[1]/ic;
    fiitrackeff = (float)iitrack/ic;
    if(fgreta>0)
      fcard29eff = (float)card29/fgreta;
    fscint = scint;
    if(scint<1){//for backwards compatability
      fscint = ic;
    }
    SetTrackEff(scint, pad, track);
  }
  void SetTrackEff(Long64_t scint, Long64_t pad[2], Long64_t track){
    fpadeff[0] = (float)pad[0]/scint;
    fpadeff[1] = (float)pad[1]/scint;
    ftrackeff = (float)track/scint;
  }
  void SetEvents(Int_t n){
    fevents = n;
  }
  void SetEntries(Int_t n){
    fentries = n;
  }
  void SetBuffer(Int_t n, Int_t nscaler){
    fbuffer = n;
    fsbuffer = nscaler;
  }
  void SetBeam(Float_t brho, Int_t mass, Int_t z){
    fbrho = brho;
    fmass = mass;
    fcharge = z;
  }
  int GetRunNr(){
    return frunnr;
  }
  double GetRunTime(){
    return fruntime;
  }
  int GetBuffers(){
    return fbuffer;
  }
  int GetSBuffers(){
    return fsbuffer;
  }
  int GetEvents(){
    return fevents;
  }
  int GetCalEntries(){
    return fentries;
  }
  Long64_t GetIntegral(int i){
    return fintegral[i];
  }
  Long64_t GetIC(){
    return fich;
  }
  float GetOBJEff(){
    return fobjeff;
  }
  float GetXFPEff(){
    return fxfpeff;
  }
  float GetTOFEff(){
    return ftofeff;
  }
  float GetCRDCEff(int i){
    return fpadeff[i];
  }
  float GetTRACKEff(){
    return ftrackeff;
  }
  float GetPPACEff(int i){
    return fppaceff[i];
  }
  float GetIITRACKEff(){
    return fiitrackeff;
  }
  float GetCard29Eff(){
    return fcard29eff;
  }
  float GetLifeTime(){
    return (double)fintegral[12]/(double)fintegral[11];
  }
  Long64_t GetTotalOBJ(){
    return fintegral[24];
  }
  Long64_t GetTotalXFP(){
    return fintegral[25];
  }
  Float_t GetBrho(){
    return fbrho;
  }
  Int_t GetMass(){
    return fmass; 
  }
  Int_t GetZ(){
    return fcharge; 
  }
  Float_t GetBeta(){
    return GetBetaGamma()/GetGamma();
  }
  Float_t GetBetaGamma(){
    return fbrho/3.107*(float)fcharge/(float)fmass;
  }
  Float_t GetGamma(){
    return sqrt(1+GetBetaGamma()*GetBetaGamma());
  }
  Float_t GetEkin(){
    return fmass*931.5016*(GetGamma()-1);
  }
  Float_t GetEkinpA(){
    return 931.5016*(GetGamma()-1);
  }
  Float_t GetP(){
    return fbrho/3.335*fcharge;
  }
  
  void Print();
protected:
  int frunnr;
  int fbuffer;
  int fsbuffer;
  int fentries;
  int fevents;
  Float_t fbrho;
  Int_t fmass;
  Int_t fcharge;
  double fruntime; 
  Long64_t fintegral[NSCALER];
  Long64_t fich;
  Long64_t fscint;
  Long64_t fhodo;
  Long64_t fgreta;
  Long64_t fcard29;
  float fobjeff;
  float fxfpeff;
  float ftofeff;
  float fpadeff[2];
  float ftrackeff;
  float fppaceff[2];
  float fiitrackeff;
  float fcard29eff;
  ClassDef(RunInfo, 1);
};

#endif
