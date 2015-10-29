#ifndef UNPACKED_EVENT_HH__
#define UNPACKED_EVENT_HH__

#include <iostream>
#include <iomanip>

#include "TObject.h"
#include "TTree.h"
#include "TList.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TEnv.h"
#include "math.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include "Calibration.hh"

#include "Gretinadefs.h"
#include "Gretina.hh"

#include "Trace.hh"

#include "S800defs.h"
#include "S800.hh"

#include "Scalerdefs.h"
#include "Scaler.hh"

#include "S800Calc.hh"
#include "GretinaCalc.hh"
#include "Mode3Calc.hh"

#include "Settings.hh"
#include "RawHistograms.hh"
#include "CalHistograms.hh"

class UnpackedEvent {
public:
  UnpackedEvent(){};
  UnpackedEvent(Settings* settings);
  ~UnpackedEvent(){
    delete fGretina;
    delete fS800;
    delete fMode3Event;
    delete fS800Calc;
    delete fMode3Calc;
    delete fGretinaCalc;
    delete fScaler;
    delete frhist;
    delete fchist;
  }
  
  void SetVL(int vl){
    fvl = vl;
  }
  void SetCalibration(Calibration* cal){
    fcal = cal;
  }
  void SetWrite(bool wtree, bool whist, bool wctree, bool wchist){
    fwtree = wtree;
    fwhist = whist;
    fwcaltree = wctree;
    fwcalhist = wchist;
  }
  void SetRecalibrate(bool recalibrate){
    frecalibrate = recalibrate;
  }
  void SetWrite(bool wtree, bool whist){
    fwtree = wtree;
    fwhist = whist;
  }
  void Init();
  //! Passed a gretina Crystal, make a new crystal in the Gretina object.
  int DecodeGretina(Crystal* cryst, long long int gts);
  //! Read the specified buffer and make the S800 event.
  int DecodeS800(unsigned short *pevent, long long int ts, unsigned short twords);
  //! Read the specified buffer, make a Scaler, and write it to the scaler tree.
  int DecodeScaler(unsigned short *pevent, long long int ts);
  //! Read the specified buffer and make the Mode3Event.
  int DecodeMode3(char* cBuf, int len, long long int ts, bool card29 = false);

  //! Write the last event to file.
  /*!
    Write the last event to file.
    This event would not be written otherwise, as there are no following buffers.
    Usually, events are written whenever the following entry is past some time window.
   */
  void WriteLastEvent();
  void PrintHit(struct crys_ips_abcd1234 inbuf);
  void PrintHit(struct crys_ips_abcd5678 inbuf);

  int NrOfEvents(){return fnentries;}
  int NrOfCalEvents(){return fncalentries;}
  int NrOfHits(){return fhits;}
  int NrOfStrangeHits(){return fstrangehits;}
  TTree* GetTree(){return ftr;}
  TTree* GetCalTree(){return fcaltr;}
  TTree* GetTraceTree(){return ftr;}
  Scaler* GetScaler(){return fScaler;};
  bool HasScaler(){return fhasscaler;}
  double GetRunTime();

protected:
  //! Create a single trace
  Trace DecodeTrace(unsigned short** wBuf_p, int length, long long int gts);
  unsigned short* DecodeS800Crdc(unsigned short *pevent, int id, bool fp);
  unsigned short* DecodeS800CrdcRaw(unsigned short *pevent, int id, bool fp);
  unsigned short* DecodeS800Ppac(unsigned short *pevent);
  unsigned short* DecodeS800PpacRaw(unsigned short *pevent);
  unsigned short* DecodeS800IonChamber(unsigned short *pevent);
  unsigned short* DecodeS800Scintillator(unsigned short *pevent, unsigned short updown, int id);
  unsigned short* DecodeS800LaBr(unsigned short* pevent);
  unsigned short* DecodeS800TimeOfFlight(unsigned short *pevent);
  unsigned short* DecodeS800Trigger(unsigned short *pevent);
  unsigned short* DecodeS800Hodoscope(unsigned short *pevent, unsigned short length);
  unsigned short* DecodeS800Galotte(unsigned short *pevent, unsigned short length);
  unsigned short* DecodeS800ObPin(unsigned short *pevent);

  //! Performs end of event actions.
  /*!
    Performs end of event actions.
    The actions performed depend on the flags given.
    - rt, write out the raw tree.
    - rh, fill raw histograms, to be written later.
    - ct, write out the calibrated tree.
    - ch, fill calibrated histograms, to be written later.

    Note that calibrations are only performed in GrROOT if either "ct" or "ch" are given as flags.
   */
  void CloseEvent();
  //! Clears memory of current event.
  void ClearEvent();

  TTree *ftr;
  TTree *fcaltr;

  Gretina *fGretina;
  S800 *fS800;
  Scaler *fScaler;
  Mode3Event *fMode3Event;

  S800Calc *fS800Calc;
  GretinaCalc* fGretinaCalc;
  Mode3Calc *fMode3Calc;

  bool fhasdata;
  bool fhasscaler;
  long long int fcurrent_ts;
  long long int ffirst_ts;
  int fvl;
  int fnentries;
  int fncalentries;
  int fhits;
  int fctr;
  bool fwtree;
  bool fwhist;
  bool fwcaltree;
  bool fwcalhist;
  int fstrangehits;
  int frecalibrate;
  int fEventTimeDiff;

  Calibration* fcal;

  RawHistograms* frhist;
  CalHistograms* fchist;

  double fslope[7][4];
  double foffset[7][4];

  Settings* fSett;
  map<int,bool> holesFound;

  ULong64_t fprevscaler[NSCALER];
  int foverflows[NSCALER];


};

#endif
