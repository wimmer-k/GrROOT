#ifndef __S800CALC_HH
#define __S800CALC_HH

#include <iostream>
#include <vector>
#include <math.h>
#include "TObject.h"
#include "TVector3.h"

using namespace std;
//CRDC

class Calibration;
class S800;

class PAD : public TObject {
public:
  PAD(){
    fx = sqrt(-1.0);
    fy = sqrt(-1.0);
    fx_gravity = sqrt(-1.0);
    fcal.clear();
    fchan.clear();
    fid = -1;
    fmaxpad =-1;
    fpadmax.resize(3);
  }
  ~PAD(){
    Clear();
  };

  void Clear(){
    fx = sqrt(-1.0);
    fy = sqrt(-1.0);
    fx_gravity = sqrt(-1.0);
    fcal.clear();
    fchan.clear();
    fid = -1;
    fmaxpad =-1;
    fpadmax.resize(3);
  }
  void SetXY(Float_t x, Float_t y){
    fx = x;
    fy = y;
  }
  //! Sets the pad values of the CRDC.
  /*!
    Sets the pad values of the CRDC.
    Any pads whose value is NaN are thrown out.
   */
  void SetCal(vector<Float_t> cal){
    for(UShort_t i=0;i<cal.size();i++){
      if(!isnan(cal[i])){
	fcal.push_back(cal[i]);
	fchan.push_back(i);
      }
    }
  }
  void SetID(int id){fid = id;}
  void SetXGravity(Float_t x_gravity){fx_gravity = x_gravity;}
  void SetXFit(Float_t x_fit){fx_fit = x_fit;}
  void SetTAC(float tac){ftac = tac;}
  void SetMaxPad(Short_t maxpad){fmaxpad = maxpad;}
  void SetPadMax(Short_t i,Float_t padmax){fpadmax[i] = padmax;}

  Int_t GetID(){return fid;}
  Float_t GetX(){return fx;}
  Float_t GetY(){return fy;}
  vector<Float_t> GetCal(){return fcal;}
  Float_t GetCal(int i){return fcal[i];}
  vector<Int_t> GetChan(){return fchan;}
  Int_t GetChan(int i){return fchan[i];}

#ifdef S800_DETAILEDTREE
  Float_t GetTAC(){return ftac;}
  Short_t GetMaxPad(){return fmaxpad;}
  vector<Float_t> GetPadMax(){return fpadmax;}
  Float_t GetXGravity(){return fx_gravity;}
  Float_t GetXFit(){return fx_fit;}
#endif
protected:
  Float_t fx;
  Float_t fy;
  Int_t fid;
#ifdef S800_DETAILEDTREE
  vector<Float_t> fcal;
  vector<Int_t> fchan;
  Float_t ftac;
  Short_t fmaxpad;
  vector<Float_t> fpadmax;
  Float_t fx_gravity;
  Float_t fx_fit;
#else
  vector<Float_t> fcal; //!
  vector<Int_t> fchan; //!
  Float_t ftac; //!
  Short_t fmaxpad; //!
  vector<Float_t> fpadmax; //!
  Float_t fx_gravity; //!
  Float_t fx_fit;  //!
#endif

  ClassDef(PAD, 1);
};
class PPAC : public TObject {
public:
  PPAC(){
    fid = -1;
    fx = sqrt(-1.0);
    fy = sqrt(-1.0);
    fxstrip.clear();
    fystrip.clear();
    fxmult = 0;
    fymult = 0;
    fxmax = 0;
    fymax = 0;
 }
  ~PPAC(){
    Clear();
  };

  void Clear(){
    fid = -1;
    fx = sqrt(-1.0);
    fy = sqrt(-1.0);
    fxstrip.clear();
    fystrip.clear();
    fxmult = 0;
    fymult = 0;
    fxmax = 0;
    fymax = 0;
  }
  //setting
  void SetXY(Float_t x, Float_t y){
    fx = x;
    fy = y;
  }
  void SetStrips(vector<Short_t> x, vector<Short_t> y){
    fxstrip = x;
    fystrip = y;
  }
  void SetMult(Short_t xm, Short_t ym){
    fxmult = xm;
    fymult = ym;
  }
  void SetMax(Short_t xm, Short_t ym){
    fxmax = xm;
    fymax = ym;
  }
  Int_t GetID(){return fid;}
  Float_t GetX(){return fx;}
  Float_t GetY(){return fy;}
#ifdef S800_DETAILEDTREE
  vector<Short_t> GetXStrips(){return fxstrip;}
  vector<Short_t> GetYStrips(){return fystrip;}
  Short_t GetXMult(){return fxmult;}
  Short_t GetYMult(){return fymult;}
  Short_t GetXMax(){return fxmax;}
  Short_t GetYMax(){return fymax;}
#endif

protected:
  Float_t fx;
  Float_t fy;
  Int_t fid;
#ifdef S800_DETAILEDTREE
  vector<Short_t> fxstrip;
  vector<Short_t> fystrip;
  Short_t fxmult;
  Short_t fymult;
  Short_t fxmax;
  Short_t fymax;
#else
  vector<Short_t> fxstrip; //!
  vector<Short_t> fystrip; //!
  Short_t fxmult; //!
  Short_t fymult; //!
  Short_t fxmax; //!
  Short_t fymax; //!
#endif
  ClassDef(PPAC, 1);
};
//IONCHAMBER
class IC : public TObject {
public:
  IC(){
    fcal.clear();
    fchan.clear();
    fsum = sqrt(-1.0);
    fde = sqrt(-1.0);
  }
  ~IC(){
    Clear();
  };
  void Clear(){
    fcal.clear();
    fchan.clear();
    fsum = sqrt(-1.0);
    fde = sqrt(-1.0);
  }
  //! Sets the calibrated data.
  /*!
    Sets the calibrated data of the IC.
    Any value passed that is NaN will be thrown out and not saved.
   */
  void SetCal(vector<Float_t> cal){
    fcal.clear();
    fchan.clear();
    for(UShort_t i=0;i<cal.size();i++){
      if(!isnan(cal[i])&&cal[i]>0.0){
	fcal.push_back(cal[i]);
	fchan.push_back(i);
      }
    }
  }
  void SetSum(Float_t sum){fsum = sum;}
  void SetDE(Float_t de){fde = de;}
  Float_t GetSum(){return fsum;}
  Float_t GetDE(){return fde;}

  vector<Float_t> GetCal(){return fcal;}
#ifdef S800_DETAILEDTREE
  vector<Int_t> GetChan(){return fchan;}
#endif

protected:
#ifdef S800_DETAILEDTREE
  vector<Float_t> fcal;
  vector<Int_t> fchan;
#else
  vector<Float_t> fcal; //!
  vector<Int_t> fchan; //!
#endif
  Float_t fsum;
  Float_t fde;
  ClassDef(IC, 1);
};

//TIMEOFFLIGHT
class TOF : public TObject {
public:
  TOF(){
    frf = sqrt(-1.0);
    fobj = sqrt(-1.0);
    fxfp = sqrt(-1.0);
    frfc = sqrt(-1.0);
    fobjc = sqrt(-1.0);
    fxfpc = sqrt(-1.0);

    ftac_obj = sqrt(-1.0);
    ftac_xfp = sqrt(-1.0);
    ftac_objc = sqrt(-1.0);
    ftac_xfpc = sqrt(-1.0);

  }
  ~TOF(){
    Clear();
  };
  void Clear(){
    frf = sqrt(-1.0);
    fobj = sqrt(-1.0);
    fxfp = sqrt(-1.0);
    frfc = sqrt(-1.0);
    fobjc = sqrt(-1.0);
    fxfpc = sqrt(-1.0);
    ftac_obj = sqrt(-1.0);
    ftac_xfp = sqrt(-1.0);
    ftac_objc = sqrt(-1.0);
    ftac_xfpc = sqrt(-1.0);
  }
  void Set(Float_t rf, Float_t obj, Float_t xfp){
    frf = rf;
    fobj = obj;
    fxfp = xfp;
  }
  void SetTAC(Float_t obj, Float_t xfp){
    ftac_obj = obj;
    ftac_xfp = xfp;
  }
  void SetCorr(Float_t rf, Float_t obj, Float_t xfp){
    frfc = frf+rf;
    fobjc = fobj+obj;
    fxfpc = fxfp+xfp;
  }
  void SetCorrVals(Float_t rfc, Float_t objc, Float_t xfpc){
    frfc = rfc;
    fobjc = objc;
    fxfpc = xfpc;
  }
  void SetTACCorr(Float_t obj, Float_t xfp){
    ftac_objc = ftac_obj+obj;
    ftac_xfpc = ftac_xfp+xfp;
  }
  void SetTACCorrVals(Float_t objc, Float_t xfpc){
    ftac_objc = objc;
    ftac_xfpc = xfpc;
  }
  Float_t GetRF(){return frf;}
  Float_t GetOBJ(){return fobj;}
  Float_t GetXFP(){return fxfp;}
  Float_t GetTDiff(){return fobj - fxfp;}
  Float_t GetTACOBJ(){return ftac_obj;}
  Float_t GetTACXFP(){return ftac_xfp;}
  Float_t GetTACTDiff(){return ftac_obj - ftac_xfp;}
  //corr
  Float_t GetRFC(){return frfc;}
  Float_t GetOBJC(){return fobjc;}
  Float_t GetXFPC(){return fxfpc;}
  Float_t GetTDiffC(){return fobjc - fxfpc;}
  Float_t GetTACOBJC(){return ftac_objc;}
  Float_t GetTACXFPC(){return ftac_xfpc;}
  Float_t GetTACTDiffC(){return ftac_objc - ftac_xfpc;}


protected:
  //! The rf timing, measured with the TDC.
  Float_t frf;
  //! The obj timing, measured with the TDC.
  Float_t fobj;
  //! The xfp timing, measured with the TDC.
  Float_t fxfp;
  //! The rf timing, measured with the TDC, corrected for angle and position.
  Float_t frfc;
  //! The obj timing, measured with the TDC, corrected for angle and position.
  Float_t fobjc;
  //! The xfp timing, measured with the TDC, corrected for angle and position.
  Float_t fxfpc;
  //! The obj timing, measured with the TAC.
  Float_t ftac_obj;
  //! The xfp timing, measured with the TAC.
  Float_t ftac_xfp;
  //! The obj timing, measured with the TAC, corrected for angle and position.
  Float_t ftac_objc;
  //! The xfp timing, measured with the TAC, corrected for angle and position.
  Float_t ftac_xfpc;
  ClassDef(TOF, 1);
};

//SCINTILLATOR
class SCINT : public TObject {
public:
  SCINT(){
    ftime = 0.0;
    fde = 0.0;
  }
  ~SCINT(){
    Clear();
  };
  void Clear(){
    ftime = 0.0;
    fde = 0.0;
  }
  void Set(Float_t time, Float_t de){
    ftime = time;
    fde = de;
  }
  void SetTime(Float_t tup, Float_t tdown){
    ftimeup = tup;
    ftimedown = tdown;
    ftime = (tup+tdown)/2.;
  }
  void SetDE(Float_t de_up,Float_t de_down){
    fdeup = de_up;
    fdedown = de_down;
    fde = sqrt(1.0*de_up*de_up + 1.0*de_down*de_down);
  }
  Float_t GetTime(){return ftime;}
  Float_t GetDE(){return fde;}

protected:
  Float_t ftime;
  Float_t fde;
#ifdef S800_DETAILEDTREE
  Float_t ftimeup;
  Float_t ftimedown;
  Float_t fdeup;
  Float_t fdedown;
#else
  Float_t ftimeup; //!
  Float_t ftimedown; //!
  Float_t fdeup; //!
  Float_t fdedown; //!
#endif
  ClassDef(SCINT, 1);
};

//TRACK
class TRACK : public TObject {
public:
  TRACK(){
    fxfp = sqrt(-1.0);
    fafp = sqrt(-1.0);
    fyfp = sqrt(-1.0);
    fbfp = sqrt(-1.0);
    fata = sqrt(-1.0);
    fyta = sqrt(-1.0);
    fbta = sqrt(-1.0);
    fdta = sqrt(-1.0);
    fazita = sqrt(-1.0);
    fscatter = sqrt(-1.0);
    fptot = sqrt(-1.0);
    fppar = sqrt(-1.0);
    fptra = sqrt(-1.0);
    fetot = sqrt(-1.0);
  }
  ~TRACK(){
    Clear();
  };
  void Clear(){
    fxfp = sqrt(-1.0);
    fafp = sqrt(-1.0);
    fyfp = sqrt(-1.0);
    fbfp = sqrt(-1.0);
    fata = sqrt(-1.0);
    fyta = sqrt(-1.0);
    fbta = sqrt(-1.0);
    fdta = sqrt(-1.0);
    fazita = sqrt(-1.0);
    fscatter = sqrt(-1.0);
    fptot = sqrt(-1.0);
    fppar = sqrt(-1.0);
    fptra = sqrt(-1.0);
    fetot = sqrt(-1.0);
  }

  Float_t GetXFP(){return fxfp;}
  Float_t GetAFP(){return fafp;}
  Float_t GetYFP(){return fyfp;}
  Float_t GetBFP(){return fbfp;}
  Float_t GetATA(){return fata;}
  Float_t GetYTA(){return fyta;}
  Float_t GetBTA(){return fbta;}
  Float_t GetDTA(){return fdta;}
  Float_t GetPhi(){return fazita;}
  Float_t GetTheta(){return fscatter;}
  Float_t GetPtot(){return fptot;}
  Float_t GetPpar(){return fppar;}
  Float_t GetPtra(){return fptra;}
  Float_t GetEtot(){return fetot;}

  void SetXFP(Float_t value){fxfp = value;}
  void SetAFP(Float_t value){fafp = value;}
  void SetYFP(Float_t value){fyfp = value;}
  void SetBFP(Float_t value){fbfp = value;}
  void SetATA(Float_t value){fata = value;}
  void SetYTA(Float_t value){fyta = value;}
  void SetBTA(Float_t value){fbta = value;}
  void SetDTA(Float_t value){fdta = value;}
  void SetPhi(Float_t value){fazita = value;}
  void SetTheta(Float_t value){fscatter = value;}
  void SetPtot(Float_t value){fptot = value;}
  void SetPpar(Float_t value){fppar = value;}
  void SetPtra(Float_t value){fptra = value;}
  void SetEtot(Float_t value){fetot = value;}

protected:
  //! The x location at CRDC1
  Float_t fxfp;
  //! The x angle at CRDC1
  Float_t fafp;
  //! The y location at CRDC1
  Float_t fyfp;
  //! The y angle at CRDC1
  Float_t fbfp;
  //! The x angle at the target. (mrad)
  Float_t fata;
  //! The y position at the target. (mm)
  Float_t fyta;
  //! The y angle at the target. (mrad)
  Float_t fbta;
  //! The relative energy difference from the central brho (% dE/E).
  Float_t fdta;
  //! The azimuthal angle of the projectile. (deg)
  Float_t fazita;
  //! The scattering angle of the projectile. (deg)
  Float_t fscatter;
  Float_t fptot;
  Float_t fppar;
  Float_t fptra;
  Float_t fetot;
  ClassDef(TRACK, 1);
};

//IITRACK
class IITRACK : public TObject {
public:
  IITRACK(){
    fxii = sqrt(-1.0);
    faii = sqrt(-1.0);
    fyii = sqrt(-1.0);
    fbii = sqrt(-1.0);
    fazita = sqrt(-1.0);
    fscatter = sqrt(-1.0);
  }
  ~IITRACK(){
    Clear();
  };
  void Clear(){
    fxii = sqrt(-1.0);
    faii = sqrt(-1.0);
    fyii = sqrt(-1.0);
    fbii = sqrt(-1.0);
    fazita = sqrt(-1.0);
    fscatter = sqrt(-1.0);
  }

  Float_t GetXII(){return fxii;}
  Float_t GetAII(){return faii;}
  Float_t GetYII(){return fyii;}
  Float_t GetBII(){return fbii;}
  Float_t GetPhi(){return fazita;}
  Float_t GetTheta(){return fscatter;}

  void SetXII(Float_t value){fxii = value;}
  void SetAII(Float_t value){faii = value;}
  void SetYII(Float_t value){fyii = value;}
  void SetBII(Float_t value){fbii = value;}
  void SetPhi(Float_t value){fazita = value;}
  void SetTheta(Float_t value){fscatter = value;}

protected:
  Float_t fxii;
  Float_t faii;
  Float_t fyii;
  Float_t fbii;
  Float_t fazita;
  Float_t fscatter;
  ClassDef(IITRACK, 1);
};

//HODO
class HODO : public TObject {
public:
  HODO(){
    Clear();
  }
  ~HODO(){
    Clear();
  };
  void Clear(){
    ftime =-1;
    fmult = 0;
    fmultab = 0;
    fen.clear();
    fch.clear();
    fenab.clear();
    fchab.clear();
    fmaxhit = 0;
  }
  void SetTime(int time){ftime = time;}
  void Set(Double_t data, Short_t ch){
    fen.push_back(data);
    fch.push_back(ch);
    fmult++;
  }
  void SetAB(Double_t data, Short_t ch){
    fenab.push_back(data);
    fchab.push_back(ch);
    fmultab++;
  }
  void AddAB(Double_t en, Short_t ch){
    fenab.back() = fenab.back()+en;
    if(en > fmaxhit)
      fchab.back() = ch;
  }
  void AddFirstAB(Double_t en, Short_t ch){
    fenab.push_back(en);
    fchab.push_back(ch);
    fmultab++;
    fmaxhit = en;
  }

  Int_t GetTime(){return ftime;}
  Short_t GetMult(){return fmult;}
  Short_t GetMultAB(){return fmultab;}
  vector<Double_t>* GetEnergy(){return &fen;}
  vector<Short_t>* GetChannel(){return &fch;}
  vector<Double_t>* GetEnergyAB(){return &fenab;}
  vector<Short_t>* GetChannelAB(){return &fchab;}

protected:
  Int_t ftime;
  Short_t fmult;
  vector<Double_t> fen;
  vector<Short_t> fch;
  Short_t fmultab;
  vector<Double_t> fenab;
  vector<Short_t> fchab;
  Double_t fmaxhit;

  ClassDef(HODO, 1);
};

//S800 CALCULATED
class S800Calc : public TObject {
public:
  S800Calc(){}
  void Clear(){
    ftimes800 =0;
    fregistr =0;
    fts =0;
    fPAD[0].Clear();
    fPAD[1].Clear();
    fPPAC[0].Clear();
    fPPAC[1].Clear();
    fIC.Clear();
    fTOF.Clear();
    fSCINT[0].Clear();
    fSCINT[1].Clear();
    fSCINT[2].Clear();
    fTRACK.Clear();
    fIITRACK.Clear();
    fHODO.Clear();
  }

  void SetTimeS800(Float_t time){ftimes800 = time;}
  Float_t GetTimeS800(){return ftimes800;};
  void SetRegistr(Short_t registr){fregistr = registr;}
  Short_t GetRegistr(){return fregistr;};
  void SetTimeStamp(long long int time){fts = time;}
  long long int GetTimeStamp(){return fts;};
  bool IsCoinc(){
    return ( fregistr & (1<<1) );
  }
  bool IsSingle(){
    return ( fregistr & (1<<0) );
  }
  bool IsHodo(){
    return ( fregistr & (1<<3) );
  }

  void SetPAD(PAD pad, int id){fPAD[id] = pad;}
  void SetPPAC(PPAC ppac, int id){fPPAC[id] = ppac;}
  void SetIC(IC ic){fIC = ic;}
  void SetTOF(TOF tof){fTOF = tof;}
  void SetSCINT(SCINT scint, int id){fSCINT[id] = scint;}
  void SetTRACK(TRACK track){fTRACK = track;}
  void SetIITRACK(IITRACK track){fIITRACK = track;}
  void SetHODO(HODO hodo){fHODO = hodo;}
  PAD* GetPAD(int id){return &fPAD[id];}
  PPAC* GetPPAC(int id){return &fPPAC[id];}
  IC* GetIC(){return &fIC;}
  TOF* GetTOF(){return &fTOF;}
  SCINT* GetSCINT(int id){return &fSCINT[id];}
  TRACK* GetTRACK(){return &fTRACK;}
  IITRACK* GetIITRACK(){return &fIITRACK;}
  HODO* GetHODO(){return &fHODO;}

protected:
  PAD fPAD[2];
  PPAC fPPAC[2];
  IC fIC;
  TOF fTOF;
  SCINT fSCINT[3];
  TRACK fTRACK;
  IITRACK fIITRACK;
  HODO fHODO;
  Float_t ftimes800;
  Short_t fregistr;
  long long int fts;

  ClassDef(S800Calc, 1);
};


#endif
