#ifndef __GRETINA_HH
#define __GRETINA_HH

#include <iostream>
#include <vector>
#include <math.h>
#include "TObject.h"
#include "TVector3.h"
#include "TMath.h"
#include "Gretinadefs.h"

using namespace std;

//interaction points
struct ip{
  float x, y, z;
  float e;
  int seg;
  float seg_ener;
};
//crystal stuff
struct crys_ips_abcd1234{
  int type;
  int crystal_id;
  int num;
  float tot_e;
  long long int timestamp;
  long long trig_time;
  float t0;
  float cfd;
  float chisq;
  float norm_chisq;
  float baseline;
  int pad;
  ip ips[MAX_INTPTS];
};
struct crys_ips_abcd5678 {
  int type;          /* defined as abcd5678 */
  int crystal_id;
  int num;           /* # of int pts from decomp, or # of nets on decomp error */
  float tot_e;       /* dnl corrected */
  int core_e[4];     /* 4 raw core energies from FPGA filter (no shift) */
  long long int timestamp;
  long long trig_time;    /* not yet impl */
  float t0;
  float cfd;
  float chisq;
  float norm_chisq;
  float baseline;
  float prestep;    /* avg trace value before step */
  float poststep;   /* avg trace value following step */
  int pad;          /* non-0 on decomp error, value gives error type */
  ip ips[MAX_INTPTS];
};


class IPoint : public TObject {
public:
  IPoint();
  IPoint(Float_t en, Float_t x, Float_t y, Float_t z, Int_t seg, Float_t segen);
  IPoint(IPoint* old);
  ~IPoint(){Clear();}
  void Clear();
  void SetEnergy(Float_t en){fen = en;}
  void SetPosition(Float_t x, Float_t y, Float_t z){fposition.SetXYZ(x,y,z);}
  void SetPosition(TVector3 pos){fposition = pos;}

  Float_t GetEnergy(){return fen;}
  Float_t GetSegEnergy(){return fseg_en;}
  Int_t GetSeg(){return fseg;}
  TVector3 GetPosition(){return fposition;}
  void PrintEvent();

protected:
  //! The energy of the interaction point.
  Float_t fen;
  //! The position relative to the crystal of the interaction point.
  TVector3 fposition;
  //! The segment number (0-35) of the interaction point.
  Int_t fseg;
  //! The total energy deposited in this segment.  (May contain multiple interaction points.)
  Float_t fseg_en;
  ClassDef(IPoint, 1);
};

class Crystal : public TObject {
public:
  Crystal();
  Crystal(crys_ips_abcd1234 inbuf);
  Crystal(crys_ips_abcd5678 inbuf);
  Crystal(Crystal* old);
  ~Crystal();
  void Clear();
  void AddBackCrystal(Crystal* other);
  void AddIP(IPoint *ip);

  void SetTrigTime(long long int trigtime){ftrig_time = trigtime;}
  void SetT0(Float_t t0){ft0 = t0;}
  void SetCFD(Float_t cfd){fcfd = cfd;}
  void SetChiSq(Float_t chisq){fchisq = chisq;}
  void SetNChiSq(Float_t norm_chisq){fnorm_chisq = norm_chisq;}
  void SetBaseline(Float_t baseline){fbaseline = baseline;}
  void SetError(UShort_t error){ferror = error;}
  Short_t GetID(){return fholenum*4 + fcrystalid;}
  Short_t GetHole(){return fholenum;}
  Short_t GetCrystal(){return fcrystalid;}
  Float_t GetEnergy(){return fen;}
  void SetEnergy(Float_t val){fen = val;}
  Float_t GetIPSum();
  Float_t GetSegmentSum();
  long long int GetTS(){return ftimestamp;}
  long long int GetITS(){return fits;}
  long long int GetTrigTime(){return ftrig_time;}
  Float_t GetT0(){return ft0;}
  Float_t GetCFD(){return fcfd;}
  Float_t GetChiSq(){return fchisq;}
  Float_t GetNChiSq(){return fnorm_chisq;}
  Float_t GetBaseline(){return fbaseline;}
  Short_t GetMult(){return fmult;}
  UShort_t GetError(){return ferror;}
  int GetCoreE(int i){return fcore_e[i];}
  vector<IPoint*> GetIPoints(){return fipoints;}
  IPoint* GetIPoint(int n){return fipoints[n];}
  Short_t GetMaxIPNr(){return fmaxip;}
  IPoint* GetMaxIP(){return fipoints[fmaxip];}
  Float_t GetMaxEn(){return fmaxen;}
  Float_t GetMaxSingleCrystal(){return fMaxSingleCrystal;}
  float GetPreStep(){return fprestep;}
  float GetPostStep(){return fpoststep;}
  void PrintEvent();

protected:
  //! The hole number (1-30) of the hit.
  Short_t fholenum;
  //! The crystal number (0-3) of the hit
  Short_t fcrystalid;
  //! The energy (keV) of the hit.
  Float_t fen;
  //! No longer used.
  Float_t fMaxSingleCrystal;
  //! The number of interaction points in the crystal.
  Short_t fmult;
  //! Empty as of 2012-08-13
  long long int  ftrig_time;
  //! The CFD time, as calculated by the decomposition.
  Float_t ft0;
  //! Empty as of 2012-08-13
  Float_t fcfd;
  //! chi square of decomposition
  Float_t fchisq;
  //! normalized chi square of decomposition
  Float_t fnorm_chisq;
  //! The baseline for the decomposition.
  Float_t fbaseline;
  vector<IPoint*> fipoints;
  //! The position in fipoints of the IP with the most deposited energy
  Short_t fmaxip;
  //! The energy of the segment with the most deposited energy.
  Float_t fmaxen;
  //! The 4 core energies.
  /*!
    The 4 core energies.
    According to Mario, for backwards compatibility reasons, the first is always the selected range.
    The remaining 3 are the other ranges, in order.
   */
  int fcore_e[4];
  //! The GEB header timestamp.
  long long int ftimestamp;
  //! The internal timestamp.  (Should always be equal to ftimestamp)
  long long int fits;
  //! A diagnostic added by Mario.
  float fprestep;
  //! A diagnostic added by Mario.
  float fpoststep;
  //! The error code reported by the decomposition.
  UShort_t ferror;

  ClassDef(Crystal, 1);
};

class Gretina : public TObject {
public:
  Gretina();
  ~Gretina(){Clear();}
  void Clear();
  void AddHit(Crystal* cry);

  int GetHitPattern(){return fhitpattern;}
  int GetMult(){return fmult;}
  vector<Crystal*> GetHits(){return fcrystals;}
  Crystal* GetHit(int n){return fcrystals[n];}
  void PrintEvent();
protected:
  //! An integer whose n-th bit is 1 iff the detector in hole n fired.
  int fhitpattern;
  //! The crystal multiplicity of the event.
  Short_t fmult;
  vector<Crystal*> fcrystals;
  ClassDef(Gretina, 1);
};


#endif
