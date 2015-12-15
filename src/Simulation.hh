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
#ifndef SIMULATION_HH__
#define SIMULATION_HH__

#include <iostream>
#include <iomanip>

#include "TObject.h"
#include "TTree.h"
#include "TEnv.h"
#include "TRandom.h"

#include "S800Calc.hh"
#include "GretinaCalc.hh"
#include "GretinaSim.hh"
#include "SimSettings.hh"
#include "UnpackedEvent.hh"

struct simresolution{
  // resolution = A*sqrt(1 + B*Energy)
  int Detector;
  int Crystal;
  double A;
  double B;
  double C;
};

struct simthreshold{
  // threshold = 0.5*( 1 + tanh( (Energy - E) / dE ) )
  int Detector;
  int Crystal;
  double E;
  double dE;
};

class Simulation : public UnpackedEvent{
public:
  Simulation(){};
  Simulation(SimSettings* setting);
  ~Simulation(){}
  //! Initialize, read parameters and setup the tree
  void Init();
  //! Read the resolutions to be applied to the simulated data
  void ReadSimResolution(const char* filename);
  //! Read the energy thresholds to be applied to the simulated data
  void ReadSimThresholds(const char* filename);
  //! Read the specified buffer and make the GretinaG4Sim event.
  int DecodeGretinaG4Sim(G4SIM_EGS*, long long int);
  //! Read the specified buffer and make the S800PhysicsData event.
  int DecodeS800PhysicsData(S800_PHYSICSDATA*, long long int);
  //! The number of simulated events
  int NrOfSimEvents(){return fnsimentries;}
  //! The tree containing the simulated data
  TTree* GetSimTree(){return fsimtr;}
  void SetWrite(bool wtree, bool whist, bool wctree, bool wchist, bool wstree = false){
    UnpackedEvent::SetWrite(wtree, whist, wctree, wchist);
    fwsimtree = wstree;
  }
  bool SimResolution(Gretina* gr);
  bool SimThresholds(Gretina* gr);
  virtual void CloseEvent();
protected:
  TRandom* fRand;
  TTree *fsimtr;
  GretinaSim* fGretinaSim;
  int fnsimentries;
  bool fwsimtree;
  SimSettings* fSett;
  vector<simresolution> fSimResolutions;
  vector<simthreshold>  fSimThresholds;

};

#endif
