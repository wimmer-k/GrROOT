#ifndef __S800_HH
#define __S800_HH

#include <iostream>
#include <vector>
#include "TObject.h"
#include "TString.h"
#include "S800defs.h"

using namespace std;

// TIME OF FLIGHT
class GTimeOfFlight : public TObject {
public:
  GTimeOfFlight(){
    frf = 0;
    fobj = 0;
    fxfp = 0;
    ftar = 0;
    ftac_obj = 0;
    ftac_xfp = 0;
  }
  ~GTimeOfFlight(){
    Clear();
  };
  void Clear(){
    frf = 0;
    fobj = 0;
    fxfp = 0;
    ftar = 0;
    ftac_obj = 0;
    ftac_xfp = 0;
  }
  void Set(Short_t rf, Short_t obj, Short_t xfp, Short_t tar){
    frf = rf;
    fobj = obj;
    fxfp = xfp;
    ftar = tar;
  }
  void SetTAC(Short_t tac_obj, Short_t tac_xfp){
    ftac_obj = tac_obj;
    ftac_xfp = tac_xfp;
  }
  Short_t GetRF(){return frf;}
  Short_t GetOBJ(){return fobj;}
  Short_t GetXFP(){return fxfp;}
  Short_t GetTAR(){return ftar;}
  Short_t GetTACOBJ(){return ftac_obj;}
  Short_t GetTACXFP(){return ftac_xfp;}

protected:
  //! Timing from rf of cyclotron
  Short_t frf;
  //! Timing of object scintillator
  Short_t fobj;
  //! Timing of A1900 extended focal plane
  Short_t fxfp;
  Short_t ftar;
  //! TAC timing of obj
  Short_t ftac_obj;
  //! TAC timing of xfp
  Short_t ftac_xfp;


  ClassDef(GTimeOfFlight, 1);
};

// TRIGGER
class GTrigger : public TObject {
public:
  GTrigger(){
    fregistr = 0;
    fs800 = 0;
    fexternal1 = 0;
    fexternal2 = 0;
    fsecondary = 0;
  }
  ~GTrigger(){
    Clear();
  };
  void Clear(){
    fregistr = 0;
    fs800 = 0;
    fexternal1 = 0;
    fexternal2 = 0;
    fsecondary = 0;
  }
  void Set(int registr, int s800, int external1, int external2, int secondary){
    fregistr = registr;
    fs800 = s800;
    fexternal1 = external1;
    fexternal2 = external2;
    fsecondary = secondary;
  }
  Short_t GetRegistr(){return fregistr;}
  Short_t GetS800(){return fs800;}
  Short_t GetExternal1(){return fexternal1;}
  Short_t GetExternal2(){return fexternal2;}
  Short_t GetSecondary(){return fsecondary;}

protected:
  //! Describes trigger setting
  Short_t fregistr;
  //! trigger on s800
  Short_t fs800;
  Short_t fexternal1;
  Short_t fexternal2;
  Short_t fsecondary;


  ClassDef(GTrigger, 1);
};

//SCINTILLATOR
class GScintillator : public TObject {
public:
  GScintillator(){
    fID =-1;
    fde_up = -1.0;
    fde_down = -1.0;
    ftime_up = -1.0;
    ftime_down = -1.0;
  }
  ~GScintillator(){
    Clear();
  };
  void Clear(){
    fID =-1;
    fde_up = -1.0;
    fde_down = -1.0;
    ftime_up = -1.0;
    ftime_down = -1.0;
  }
  void SetID(int id){
    fID = id;
  }
  void Set(int de_up, int time_up, int de_down, int time_down){
    if(de_up>-1)
      fde_up = de_up;
    if(time_up>-1)
      ftime_up = time_up;
    if(de_down>-1)
      fde_down = de_down;
    if(time_down>-1)
      ftime_down = time_down;
  }
  Int_t GetID(){return fID;}
  Float_t GetDE_up(){return fde_up;}
  Float_t GetTime_up(){return ftime_up;}
  Float_t GetDE_down(){return fde_down;}
  Float_t GetTime_down(){return ftime_down;}

protected:
  Int_t fID;
  //! Energy loss for top PMT
  Float_t fde_up;
  //! Energy loss for bottom PMT
  Float_t fde_down;
  Float_t ftime_up;
  Float_t ftime_down;

  ClassDef(GScintillator, 1);
};

class GLaBr : public TObject {
public:
  GLaBr(){
    Clear();
  }
  ~GLaBr(){
    Clear();
  };

  void Clear(){
    for(int i=0; i<S800_LABR_DETS; i++){
      fADC[i] = 0;
      fTDC[i] = 0;
    }
  }
  void SetADC(int ch, short val){
    if (ch<0 || ch>=S800_LABR_DETS){
      cout << "Improper LaBr detector number: " << ch << endl;
    } else {
      fADC[ch] = val;
    }
  }
  void SetTDC(int ch, short val){
    if (ch<0 || ch>=S800_LABR_DETS){
      cout << "Improper LaBr detector number: " << ch << endl;
    } else {
      fTDC[ch] = val;
    }
  }
  short GetADC(int ch){return fADC[ch];}
  short GetTDC(int ch){return fTDC[ch];}

protected:
  short fADC[S800_LABR_DETS];
  short fTDC[S800_LABR_DETS];

  ClassDef(GLaBr, 1);
};

//IONCHAMBER
class GIonChamber : public TObject {
public:
  GIonChamber(){
    fchannels.clear();
    fdata.clear();
  }
  ~GIonChamber(){
    Clear();
  };
  void Clear(){
    fchannels.clear();
    fdata.clear();
  }
  void Set(int ch, int data){
    fchannels.push_back(ch);
    fdata.push_back(data);
  }

  vector<int>* GetChannels(){return &fchannels;}
  vector<float>* GetData(){return &fdata;}

protected:
  //! A vector of which pads fired.
  vector<int> fchannels;
  //! The energies recorded by the pad corresponding to each entry in fchannels.
  vector<float> fdata;

  ClassDef(GIonChamber, 1);
};
//CRDC
class GCrdc : public TObject {
public:
  GCrdc(){
    fID =-1;
    fanode = -1;
    ftac = -1;
    fdata.clear();
    fsample.clear();
    fchannels.clear();
  }
  ~GCrdc(){
    Clear();
  };
  void Clear(){
    fID =-1;
    fanode = -1;
    ftac = -1;
    fdata.clear();
    fsample.clear();
    fchannels.clear();
  }
  void SetID(int id){fID = id;}
  void SetAnodeTAC(int anode, int tac){
    fanode = anode;
    ftac = tac;
  }
  void SetSampleWidth(int width){fwidth = width;}
  void Set(Short_t data, Short_t sample, Short_t ch){
    fdata.push_back(data);
    fsample.push_back(sample);
    fchannels.push_back(ch);
  }
  Float_t GetAnode(){return fanode;}
  Float_t GetTAC(){return ftac;}
  Int_t GetID(){return fID;}
  vector<Short_t>* GetData(){return &fdata;}
  vector<Short_t>* GetSample(){return &fsample;}
  vector<Short_t>* GetChannels(){return &fchannels;}

protected:
  //! ID of the CRDC.  (0 for CRDC1, 1 for CRDC2)
  Int_t fID;
  //! Total charge in CRDC
  Float_t fanode;
  //! Drift time of the electrons
  Float_t ftac;
  Short_t fwidth;
  //! Which pads fired.
  vector<Short_t> fchannels;
  //! The energy read for each pad.  In the same order as fchannels.
  vector<Short_t> fdata;
  //! How many samples were taken with flash ADC.  Not implemented, and is read from the settings file.
  vector<Short_t> fsample;

  ClassDef(GCrdc, 1);
};

//T PPAC
class GTppac : public TObject {
public:
  GTppac(){
    fID =-1;
    fdata.clear();
    fsample.clear();
    fchannels.clear();
  }
  ~GTppac(){
    Clear();
  };
  void Clear(){
    fID =-1;
    fdata.clear();
    fsample.clear();
    fchannels.clear();
  }
  void SetID(int id){fID = id;}
  void SetSampleWidth(int width){fwidth = width;}
  void Set(Short_t data, Short_t sample, Short_t ch){
    fdata.push_back(data);
    fsample.push_back(sample);
    fchannels.push_back(ch);
  }
  vector<Short_t>* GetData(){return &fdata;}
  vector<Short_t>* GetSample(){return &fsample;}
  vector<Short_t>* GetChannels(){return &fchannels;}

protected:
  Int_t fID;
  Short_t fwidth;
  vector<Short_t> fdata;
  vector<Short_t> fsample;
  vector<Short_t> fchannels;

  ClassDef(GTppac, 1);
};

//HODOSCOPE
class GHodoscope : public TObject {
public:
  GHodoscope(){
    ftime =-1;
    fhitpattern[0] = 0;
    fhitpattern[1] = 0;
    fdata.clear();
    fchannels.clear();
  }
  ~GHodoscope(){
    Clear();
  };
  void Clear(){
    ftime =-1;
    fhitpattern[0] = 0;
    fhitpattern[1] = 0;
    fdata.clear();
    fchannels.clear();
  }
  void SetTime(int time){ftime = time;}
  void SetHitpattern(int hp0, int hp1){
    fhitpattern[0] = hp0;
    fhitpattern[1] = hp1;
  };
  void Set(Short_t data, Short_t ch){
    fdata.push_back(data);
    fchannels.push_back(ch);
  }
  vector<Short_t>* GetData(){return &fdata;}
  vector<Short_t>* GetChannels(){return &fchannels;}
  Int_t GetTime(){return ftime;}
  Int_t GetHitPattern(int i){
    if(i>-1&&i<2)
      return fhitpattern[i];
    return 0;
  }

protected:
  Int_t ftime;
  Int_t fhitpattern[2];
  vector<Short_t> fdata;
  vector<Short_t> fchannels;

  ClassDef(GHodoscope, 1);
};

//GALOTTE
class GGalotte : public TObject {
public:
  GGalotte(){
    for(int i=0;i<4;i++)
      fdata[i] = -1;
  }
  ~GGalotte(){
    Clear();
  };
  void Clear(){
    for(int i=0;i<4;i++)
      fdata[i] = -1;
  }
  void Set(Short_t data, Short_t ch){
    if(ch>-1 && ch<4)
      fdata[ch] = data;
  }
  int GetData(int i){return fdata[i];}

protected:
  Int_t fdata[4];
  ClassDef(GGalotte, 1);
};

class S800 : public TObject {
public:
  S800(){
    Clear();
  }
  ~S800(){
    Clear();
  }
  void Clear(){
    fTof.Clear();
    fTrigger.Clear();
    for(int i=0;i<3;i++)
      fScintillator[i].Clear();
    fIonChamber.Clear();
    for(int i=0;i<2;i++)
      fCrdc[i].Clear();
    for(int i=0;i<2;i++)
      fTCrdc[i].Clear();
    fTPpac.Clear();
    fHodoscope.Clear();
    fts = -1;
  }
  void SetTS(long long int ts){fts = ts;}
  void SetInternalTS(long long int ts){fits = ts;}
  void SetEvtNr(long long int nr){fevtnr = nr;}
  void SetObPin(unsigned short i){fObPin = i;}
  GTimeOfFlight* GetTimeOfFlight(){return &fTof;}
  GTrigger* GetTrigger(){return &fTrigger;}
  GScintillator* GetScintillator(int id){return &fScintillator[id];}
  GIonChamber* GetIonChamber(){return &fIonChamber;}
  GCrdc* GetCrdc(int id){return &fCrdc[id];}
  GCrdc* GetTCrdc(int id){return &fTCrdc[id];}
  GTppac* GetTPpac(){return &fTPpac;}
  GLaBr* GetLaBr(){return &fLaBr;}
  GHodoscope* GetHodoscope(){return &fHodoscope;}
  GGalotte* GetGalotte(){return &fGalotte;}
  unsigned short GetObPin(){return fObPin;}
  long long int GetTS(){return fts;}
  long long int GetInternalTS(){return fits;}
  long long int GetEvtNr(){return fevtnr;}


protected:

  GTimeOfFlight fTof;
  GLaBr fLaBr;
  GTrigger fTrigger;
  GScintillator fScintillator[3];
  GIonChamber fIonChamber;
  GCrdc fCrdc[2];
  GCrdc fTCrdc[2];
  GTppac fTPpac;
  GHodoscope fHodoscope;
  GGalotte fGalotte;
  unsigned short fObPin;
  //! The timestamp as read from the GEB header.
  long long int fts; //global
  //! The timestamp as read from the internal data structure.
  /*!
    The timestamp as read from the internal data structure.
    fts = 8*fits, since the S800 has a 12.5 MHz clock as opposed to gretina's 100 MHz clock.
   */
  long long int fits;
  //! The event number of the S800.
  /*!
    The event number of the S800.
    This is counted internally by the S800 each time it fires.
   */
  long long int fevtnr;
  ClassDef(S800, 1);
};

#endif
