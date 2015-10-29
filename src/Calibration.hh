#ifndef __CALIBRATION_HH
#define __CALIBRATION_HH

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <map>

#include "TSystem.h"
#include "TEnv.h"
#include "TRandom.h"
#include "TMath.h"
#include "TFile.h"
#include "TH1F.h"

#include "Settings.hh"
#include "S800Calc.hh"
#include "Mode3Calc.hh"
#include "Gretinadefs.h"
#include "Fit.hh"
#include "Tracking.hh"

using namespace std;

class GTimeOfFlight;
class GIonChamber;
class GHodoscope;
class GScintillator;
class SCINT;
class IITRACK;
class TRACK;
class Gretina;
class Crystal;
class IPoint;
class Mode3Event;
class Trace;
class GretinaCalc;
class HitCalc;
class GCrdc;
struct recal;

static const unsigned short tppacmapx[64]=
{ 30, 31, 28, 29,   26, 27, 24, 25,    22, 23, 20, 21,   18, 19, 16, 17,  
	14, 15, 12, 13,   10, 11,  8,  9,     6,  7,  4,  5,    2,  3,  0,  1,
	33, 32, 35, 34,   37, 36, 39, 38,    41, 40, 43, 42,   45, 44, 47, 46,
	49, 48, 51, 50,   53, 52, 55, 54,    57, 56, 59, 58,   61, 60, 63, 62};
static const unsigned short tppacmapy[64]=
{  0,  1,  2,  3,    4,  5,  6,  7,     8,  9, 10, 11,   12, 13, 14, 15,
	16, 17, 18, 19,   20, 21, 22, 23,    24, 25, 26, 27,   28, 29, 30, 31,
  63, 62, 61, 60,   59, 58, 57, 56,    55, 54, 53, 52,   51, 50, 49, 48,
  47, 46, 45, 44,   43, 42, 41, 40,    39, 38, 37, 36,   35, 34, 33, 32};

class Calibration {
public:
  Calibration();
  
//! Constructs the Calibration object based on the given settings file.
/*!
  Constructs the Calibration objects, reading in the various settings.
  @param setting The Settings object from which to read calibration data.
  @param event The event number from which to start counting, used for time-dependent corrections.
 */
  Calibration(Settings*, int event =0);
  ~Calibration();
  //! Read the CRDC calibrations as given from the calfile and pedestalfile.
  /*!
    In the following lines, %d refers to the CRDC number, either 0 or 1,
      and %03d refers to the pad number, 0-127, padded with leading zeros to three digits.
    For each pad, the following line is expected in pedestalfile to give the pedestal values.
      - Crdc.%d.Ped.%03d: 0.0
    
    For each pad, the following lines are expected in filename to give the pad calibration.
      - Crdc.%d.Slope.%03d:  1.0
      - Crdc.%d.Offset.%03d: 0.0
   */
  void ReadCalibration(const char* filename, const char* pedestalfile);
  //! Read the GretinaRecalibration file.
  /*!
    For each file to be recalibrated, there should be a line as follows.
      - Detector.%d.Crystal.%d.Recal: 1

    Next, have two lines for the slope and offset to be used in recalibration.
      - Detector.%d.Crystal.%d.Slope:
      - Detector.%d.Crystal.%d.Offset:
  */
  void ReadGretinaRecalibration(const char* filename);
  void ReadBadPads(const char* filename);
  //! Reads the calibration of the ion chamber.
  /*!
    For each ion chamber pad, numbered 0-15, there should be the following lines.
      - IonChamber.Slope.%d:  1.0
      - IonChamber.Offset.%d: 0.0

    In addition, the total energy can be calibrated with the following parameters.
      - IonChamber.Slope.DE:  1.0
      - IonChamber.Offset.DE: 0.0
   */
  void ReadICCalibration(const char *filename);
  //! Read the time dependent corrections
  void ReadTimeCorrections(const char *filename);
  //! Read the calibration for the hodoscope
  /*!
    For each hodoscope segment, numbered 0-31, the following parameters give the calibration.
      - Hodoscope.Slope.%d:  1.0
      - Hodoscope.Offset.%d: 0.0
  */
  void ReadHodoCalibration(const char *filename);
  //! Read the matrix by which crystal coordinates are converted to lab coordinates.
  void ReadMatrix(const char* filename);
  void PrintCalibration();
  //! Read the S800 map
  /*!
    Read an S800 inverse map for use later.
    The map is expected to be in the ASCII output format as given by http://maps.nscl.msu.edu/~s800maps/
   */
  void ReadMap(const char* filename);
  //! From a given input, use the inverse map for the S800.
  Float_t MapCalc(int calcorder, int parameter, Float_t *input);


  //! Construct all calibrated objects.
  /*!
    First three parameters are the three raw data structures, S800, Gretina, and Mode3Event.
    Last three parameters are the three calibrated data structures to write to.
    Internally, calls BuildS800Calc(), BuildGretinaCalc(), and BuildMode3Calc().
    Additionally, performs calibrations that depend on multiple systems, such as the doppler correction.
    @param inS800 The raw S800 object as input.
    @param inGret The raw Gretina object as input
    @param inMode3 The raw Mode3Event object as input
    @param outS800 A pointer to the S800Calc object to be built.
    @param outGret A pointer to the GretinaCalc object to be built.
    @param outMode3 A pointer to the Mode3Calc object to be built.
   */
  void BuildAllCalc(S800* inS800, Gretina* inGret, Mode3Event* inMode3,
		    S800Calc* outS800, GretinaCalc* outGret, Mode3Calc* outMode3);

  //! Construct calibrated S800 object.
  void BuildS800Calc(S800* in, S800Calc* out);
  //! Construct the calibrated TRACK object.
  /*!
    First four parameters are the calibrated x and y positions in the CRDCs.
    Last parameter is the TRACK object to be built.
    Uses the inverse map of the S800 to calculate both focal plane and target parameters.
   */
  void BuildTrackCalc(Float_t x0, Float_t y0, Float_t x1, Float_t y1,TRACK* out);
  //! Contruct fcrdccal containing each calibrated pad.
  void CrdcCal(vector<Short_t>* channel, vector<Short_t>* data, int nr);
  //! Construct the calibrated PPAC objects, fppac[2].  Currently not used.
  void PpacSort(vector<Short_t>* channel, vector<Short_t>* data, PPAC* out1, PPAC* out2);
  void BuildIITrackCalc(Float_t x0, Float_t y0, Float_t x1, Float_t y1,IITRACK* out);
  bool BuildIonChamber(GIonChamber* in, TRACK* track, IC* out);
  //! Construct a vector containing the calibrated IC pad data.
  vector<Float_t> ICCal(GIonChamber* in);
  //! Return the average of the given IC pad data.
  Float_t ICSum(vector<Float_t> cal);
  Float_t ICDE(Float_t sum, TRACK* track);
  void BuildHodoscope(GHodoscope* in, HODO* out);
  bool BuildScint(GScintillator* in, SCINT* out);
  //! Subtract the two times.  (Why do I even have this function?)
  Float_t TimeOffset(Float_t time1, Float_t time2);
  //! Return the difference of the given time and the S800 time.
  Float_t TimeOffset(Float_t time);
  Float_t RFCorrection(TRACK* tr);
  Float_t OBJCorrection(TRACK* tr);
  Float_t XFPCorrection(TRACK* tr);
  Float_t TACOBJCorrection(TRACK* tr);
  Float_t TACXFPCorrection(TRACK* tr);

  //! Perform a gravity fit of the CRDC data.
  Float_t CalcGravityClassic();
  //! Perform a gravity fit of the CRDC data.
  Float_t CalcGravity();
  //! Fit the CRDC data, return the centroid.
  Float_t CalcFit();

  //! Build the GretinaCalc object, given a raw Gretina object.
  void BuildGretinaCalc(Gretina* in, GretinaCalc* out);
  //! Read the neighbor file, defining which Gretina crystals are adjacent.
  void ReadNeighbors(const char* filename);
  //! Read the energy calibration of Gretina.  (Used only for mode 3 data.)
  /*!
    Read the energy calibration from the file specified.
    Store the calibration in fgslope[7][4] and fgoffset[7][4]
    This is only used when reading mode 3 data, as mode 2 data already has the energies calibrated.
   */
  void ReadGretinaCal(const char* filename);

  //! Build the Mode3Calc object, given a raw Mode3Event object.
  /*
    Construct the calibrated Mode3Calc object.
    For each Trace present in the Mode3Event, construct a Mode3Detection in the Mode3Calc.
    Traces with hole-number equal to 31 are discarded, as these are card29 events.
   */
  bool BuildMode3Calc(Mode3Event* in, Mode3Calc* out);
  //! Calibrate a single Trace, adding it to the specified Mode3Calc.
  bool Mode3Calibration(Trace* trace,Mode3Calc* out);
  void AddBackMode3(Mode3Calc* mode3);
  void AddBackMode3Hole(Mode3Calc* mode3);
  void AddBackMode3Neighbors(Mode3Calc* mode3);
  void AddBackMode3Everything(Mode3Calc* mode3);
  float CalibratedEnergy(Trace* trace);

  //! Modify the input Gretina object.
  /*!
    The Gretina object is largely calibrated after the decomposition.
    However, it may be necessary to modify values.
    For example, the calibration used in a crystal may be incorrect.
    In this case, the crystal should be recalibrated as part of the overall calibration process.
   */
  bool GretinaRecalibrate(Gretina* gr);
  vector<HitCalc*> FilterOverflows(GretinaCalc* gr);
  void AddBackGretinaCrystal(GretinaCalc* gr);
  void AddBackGretinaHole(GretinaCalc* gr);
  void AddBackGretinaNeighbors(GretinaCalc* gr);
  void AddBackGretinaNeighborsNoSearch(GretinaCalc* gr);
  void AddBackGretinaEverything(GretinaCalc* gr);
  void ClusterGretina(GretinaCalc* gr, Gretina *in);
  vector<HitCalc*> ExtractAllHits(Gretina* in);
  void AllGretinaHits(GretinaCalc* gr, Gretina *in);

  //! Construct tracked gamma events
  void GammaTrack(GretinaCalc* gr, GretinaEvent* gt);

  //! Return whether a specified CRDC channel is bad.
  /*!
    The parameter is the channel in the CRDC to be checked.
    @param crdcnum The number of the CRDC to be checked.  (0 for CRDC1, 1 for CRDC2)
    @param ch The channel of the CRDC to be checked.
   */
  bool IsBad(int crdcnum,int ch){
    for(UShort_t b=0;b<fbad[crdcnum].size();b++){
      if(ch==fbad[crdcnum][b])
	return true;
    }
    return false;
  };
  //! Builds a calibrated PAD object in the class variable fpad.
  void SetPad(GCrdc* in);
  void SetTS800(Short_t ts800){
    fts800 = ts800;
    fs800valid = -1;
  }
  PAD GetPad(){return fpad;}
  void BuildTof(GTimeOfFlight* in, TRACK* track, TOF* out);

  int S800Valid(){return fs800valid;}
  Float_t GetBrho(){return fbrho;}
  Int_t GetMass(){return fmass; }
  Int_t GetZ(){return fcharge; }
  Mode3Calc* GetMode3Calc(){return fmode3;}

  void SetHasMap(bool hasmap){fhasmap = hasmap;}

  Long64_t GetICHCtr(){return fichctr;}
  Long64_t GetSCINTCtr(){return fscintctr;}
  Long64_t GetOBJCtr(){return fobjctr;}
  Long64_t GetXFPCtr(){return fxfpctr;}
  Long64_t GetTOFCtr(){return ftofctr;}
  Long64_t GetPADCtr(int pad){return fpadctr[pad];}
  Long64_t* GetPADCtr(){return fpadctr;}
  Long64_t GetTRACKCtr(){return ftrackctr;}
  Long64_t GetPPACCtr(int ppac){return fppacctr[ppac];}
  Long64_t* GetPPACCtr(){return fppacctr;}
  Long64_t GetIITRACKCtr(){return fiitrackctr;}
  Long64_t GetCARD29Ctr(){return fcard29ctr;}
  Long64_t GetGRETACtr(){return fgretactr;}
  Long64_t GetHodoCtr(){return fhodoctr;}
  
  void PrintCtrs();

  
private:
  void ResetCtrs();
  bool IsNeighbor(int d1,int c1, int d2, int c2);
  bool IsNeighbor(int ID1, int ID2);
  int NumNeighbors(int det,int cry);
  void TransformCoordinates(HitCalc* hit);
  TVector3 TransformCoordinates(int hole, int cry, TVector3 local);
  void CalibrateIPoints(Crystal* cry);
  bool HodoNeighbor(Short_t ch1, Short_t ch2);

  Settings* fSett;
  int fverbose;
  TRandom* fRand;

  int fevent;

  PAD fpad;
  vector<vector<Float_t> > fped;  
  vector<vector<Float_t> > fslope;
  vector<vector<Float_t> > foffset;
  vector<Float_t> fICoffset;
  vector<Float_t> fICslope;
  Float_t fde_slope;
  Float_t fde_offset;
  vector<vector<Int_t> > fbad;
  vector<Float_t> fcrdccal;
  vector<Short_t> fmaxcoefficient;
  vector<vector<Short_t> > forder;
  vector<vector<vector<Short_t> > > fexponent;
  vector<vector<Float_t> > fcoefficient;
  Short_t fts800;  
  Short_t fmaxorder;  
  Float_t fbrho; 
  Int_t fmass; 
  Int_t fcharge; 
  int fs800valid;

  //! time dependent corrections
  TH1F* fobj_cor;
  TH1F* fxfp_cor;
  TH1F* fobjtac_cor;
  TH1F* fxfptac_cor;
  TH1F* fy_cor[2];
  TH1F* fIC_cor[2];

  //! The matrix used to transform from crystal coordinates to lab coordinates.
  float fcrmat[MAXDETPOS][MAXCRYSTALNO][4][4];

  //! Specifies which crystals should be recalibrated and the slope and offset to be used.
  vector<recal> fGretRecal;

  //gretina

  Mode3Calc* fmode3;

  int fnrofneigh[7][4];
  int fneighbors[7][4][6];
  double fgslope[7][4];
  double fgoffset[7][4];

  
  double fbeta;
  int fAddBackType;
  bool fhasmap;

  //hodo
  double fhodogain[32];
  double fhodooffset[32];

  //counters for efficiencies
  Long64_t fichctr;
  Long64_t fscintctr;
  Long64_t fobjctr;
  Long64_t fxfpctr;
  Long64_t ftofctr;
  Long64_t fpadctr[2];
  Long64_t ftrackctr;
  Long64_t fppacctr[2];
  Long64_t fiitrackctr;
  Long64_t fcard29ctr;
  //you have to look into the future! maybe sometime there will be Greta
  Long64_t fgretactr;
  Long64_t fhodoctr;

  Tracking* ftracking;

};

#endif
