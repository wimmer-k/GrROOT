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
#ifndef __SETTINGS_HH
#define __SETTINGS_HH

#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <vector>

#include "TSystem.h"
#include "TEnv.h"
#include "TVector3.h"

using namespace std;

class Settings : public TObject {
public:
  Settings();//default ctor
  Settings(const char* filename);
  Settings(vector<char*> files);
  ~Settings();

  void ReadSettings(TEnv* set);
  void PrintSettings();


  const string InputFile(){return fInputFiles[0];}

  int VLevel(){return fVerboseLevel;}
  int EventTimeDiff(){return fEventTimeDiff;}
  int LastEvent(){return fLastEvent;}
  int CrdcChannels(){return fCrdcChannels;}
  int PpacChannels(){return fPpacChannels;}
  int PpacWidth(){return fPpacWidth;}
  float PpacPitch(){return fPpacPitch;}
  float PpacGap(){return fPpacGap;}
  float PpacZ(){return fPpacZ;}
  int SampleWidth(){return fSampleWidth;}
  int Saturation(int ch){return fSaturation[ch];}
  int MSaturation(int ch){return fMSaturation[ch];}
  int IonChamberChannels(){return fIonChamberChannels;}
  int GravityWidth(int ch){return fGravityWidth[ch];}
  int FitWidth(int ch){return fFitWidth[ch];}
  int Method(int ch){return fMethod[ch];}
  Float_t XOffset(int ch){return fxOffset[ch];}
  Float_t XSlope(int ch){return fxSlope[ch];}
  Float_t YOffset(int ch){return fyOffset[ch];}
  Float_t YCoincHeadstart(int ch){return fyCoincHeadstart[ch];}
  Float_t YSlope(int ch){return fySlope[ch];}
  Float_t YCorOffset(int ch){return fyCorOffset[ch];}
  Float_t YCorSlope(int ch){return fyCorSlope[ch];}
  const char* CalFile(){return fCalFile.c_str();}
  const char* PedestalFile(){return fPedestalFile.c_str();}
  const char* BadFile(){return fBadFile.c_str();}
  const char* MapFile(){return fMapFile.c_str();}
  const char* CalFileIC(){return fCalFileIC.c_str();}
  int TrackParams(){return fTrackParams;}
  int TrackCoefs(){return fTrackCoefs;}
  Float_t Gap(){return fGap;}
  Float_t FPShift(){return fFPShift;}
  Float_t ShiftY(){return fShiftY;}
  Float_t AngleA(){return fAngleA;}
  Float_t AngleB(){return fAngleB;}
  Float_t GretinaAngleA(){return fGretAngleA;}
  Float_t GretinaAngleB(){return fGretAngleB;}
  Float_t RFe1C(){return frfE1C;}
  Float_t RFC(){return frfC;}
  Float_t OBJe1C(){return fobjE1C;}
  Float_t OBJC(){return fobjC;}
  Float_t XFPe1C(){return fxfpE1C;}
  Float_t XFPC(){return fxfpC;}
  Float_t TACOBJe1C(){return ftac_objE1C;}
  Float_t TACOBJC(){return ftac_objC;}
  Float_t TACXFPe1C(){return ftac_xfpE1C;}
  Float_t TACXFPC(){return ftac_xfpC;}
  Float_t dE_ytilt(){return fdE_ytilt;}
  Float_t dE_x0tilt(){return fdE_x0tilt;}
  Float_t dE_xtilt(){return fdE_xtilt;}
  Float_t ICThresh(){return fICthresh;}
  const char* TimeCorFile(){return fTimeCorFile.c_str();}
  const char* EvtNrFile(){return fEvtNrFile.c_str();}
  Float_t TargetX(){return ftargetX;}
  Float_t TargetY(){return ftargetY;}
  Float_t TargetZ(){return ftargetZ;}
  TVector3 TargetPos(){return TVector3(ftargetX,ftargetY,ftargetZ);}
  Float_t TargetBeta(){return ftargetBeta;}
  const char* AccepFile(){return fAccepFile.c_str();}
  float AccepCutOff(int i){return fAccepCutOff[i];}
  float AccepCutOff(){return fAccepCutOff[1];}
  float AccepExtrapol(){return fAccepCutOff[0];}
  float DTACorr(){return fDTACorr;}
  const char* MatrixFile(){return fMatrixFile.c_str();}
  const char* NeighborFile(){return fNeighborFile.c_str();}
  const char* GretinaCalFile(){return fGretinaCalFile.c_str();}
  const char* GretinaRecalFile(){return fGretinaRecalFile.c_str();}
  int AddBackType(){return fAddBackType;}
  double ClusterAngle(){return fClusterAngle;}
  int StoreAllIPoints(){return fStoreAllIPoints;}
  void SetAddBackType(int val){fAddBackType = val;}
  double OverflowThreshold(){return fOverflowThreshold;}

  int Hole2Det(int hole);
  int Det2Hole(int det);

  //hodoscope
  const char* HodoCalFile(){return fHodoCalFile.c_str();}
  Float_t HodoThresh(){return fHodoThresh;}

  //scintillator
  double GetScintThresh(){return fScintThresh;}

  void SetTracking(bool tracking){fTracking = tracking;}
  bool GetTracking(){return fTracking;}

  //histo ranges
  int GetRangeOBJ(int i){return fobj_range[i];}
  int GetRangeXFP(int i){return fxfp_range[i];}
  int GetRangeTACOBJ(int i){return ftacobj_range[i];}
  int GetRangeTACXFP(int i){return ftacxfp_range[i];}
  int GetRangeIC(int i){return fIC_range[i];}
  int GetRangeOBJC(int i){return fobjC_range[i];}
  int GetRangeXFPC(int i){return fxfpC_range[i];}
  int GetRangeTACOBJC(int i){return ftacobjC_range[i];}
  int GetRangeTACXFPC(int i){return ftacxfpC_range[i];}
  double GetRangePP(int i){return fPP_range[i];}

protected:
  int fEventTimeDiff;

  int fVerboseLevel;
  int fLastEvent;
  int fCrdcChannels;
  int fSampleWidth;
  int fPpacChannels;
  int fPpacWidth;
  float fPpacPitch;
  float fPpacGap;
  float fPpacZ;
  int fIonChamberChannels;
  vector<string> fInputFiles;
  string fCalFile;
  string fPedestalFile;
  string fBadFile;
  string fMapFile;
  string fCalFileIC;

  int fSaturation[2];
  int fMSaturation[2];
  Int_t fGravityWidth[2];
  Int_t fFitWidth[2];
  Int_t fMethod[2];
  Float_t fxOffset[2];
  Float_t fxSlope[2];
  Float_t fyOffset[2];
  Float_t fyCoincHeadstart[2];
  Float_t fySlope[2];
  Float_t fyCorOffset[2];
  Float_t fyCorSlope[2];
  int fTrackParams;
  int fTrackCoefs;
  Float_t fGap;
  Float_t fFPShift;
  Float_t fShiftY;
  Float_t fAngleA;
  Float_t fAngleB;
  Float_t fGretAngleA;
  Float_t fGretAngleB;
  Float_t frfE1C;
  Float_t frfC;
  Float_t fobjE1C;
  Float_t fobjC;
  Float_t fxfpE1C;
  Float_t fxfpC;
  Float_t ftac_objE1C;
  Float_t ftac_objC;
  Float_t ftac_xfpE1C;
  Float_t ftac_xfpC;
  Float_t fdE_ytilt;
  Float_t fdE_xtilt;
  Float_t fdE_x0tilt;

  Float_t fICthresh;

  string fTimeCorFile;
  string fEvtNrFile;

  Float_t ftargetX;
  Float_t ftargetY;
  Float_t ftargetZ;
  Float_t ftargetBeta;

  string fAccepFile;
  float fAccepCutOff[2];
  float fDTACorr;

  //gretina
  string fMatrixFile;
  string fNeighborFile;
  string fGretinaCalFile;
  string fGretinaRecalFile;

  map<int,int> fdet2hole;
  map<int,int> fhole2det;

  int fAddBackType;
  double fClusterAngle;
  int fStoreAllIPoints;
  double fOverflowThreshold;

  //hodo
  string fHodoCalFile;
  Float_t fHodoThresh;

  //E1 scintillator
  double fScintThresh;

  bool fTracking;

  //histo ranges
  int fobj_range[3];
  int fxfp_range[3];
  int ftacobj_range[3];
  int ftacxfp_range[3];
  int fIC_range[3];
  int fobjC_range[3];
  int fxfpC_range[3];
  int ftacobjC_range[3];
  int ftacxfpC_range[3];
  double fPP_range[3];
  
  ClassDef(Settings, 1)
};

#endif
