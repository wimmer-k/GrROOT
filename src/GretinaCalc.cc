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
#include "GretinaCalc.hh"

#include "S800Calc.hh"
#include "Settings.hh"

void HitCalc::DopplerCorrect(Settings* set, TRACK* track){
  fDCen = fen*HitCalc::DopplerCorrectionFactor(GetPosition(),set,track);
  fDCen_simcheat = fen*HitCalc::DopplerCorrectionFactor(fTrueFirst,set,track);
}

double HitCalc::DopplerCorrectionFactor(TVector3 PosToTarget, Settings* set, TRACK* track){

#ifdef SPECTCL_MODE
  //Find the angle between the particle and the gamma.
  PosToTarget.SetY(PosToTarget.Y() + track->GetYTA());
  double coor_theta = PosToTarget.Theta();
  double coor_phi = PosToTarget.Phi();
  
  double s800_theta = track->GetTheta(); //In rad, already has ata,bta correction applied.
  double s800_phi = track->GetPhi(); //In rad, already has ata,bta correction applied.
  double CosDop = 
    sin(coor_theta)*sin(s800_theta) *
    (sin(coor_phi)*sin(s800_phi)+cos(coor_phi)*cos(s800_phi)) +
    cos(coor_theta)*cos(s800_theta);
  
  //Adjust the beta according to the dta from the map.
  double beta = set->TargetBeta();
  double gamma = 1/sqrt(1.0 - beta*beta);
  double dp_p = gamma/(1+gamma)*(track->GetDTA()/100); //DTA is in % dE/E when stored in TRACK
  beta *= (1+dp_p/(gamma*gamma));
  gamma = 1/sqrt(1-beta*beta);
  
  return gamma*(1-beta*CosDop);

#else
  PosToTarget.SetY(PosToTarget.Y() - track->GetYTA());
  TVector3 BeamDir = TVector3(0,0,1);
  //theta and phi do have the y -> -y rotation already in ata and bta not
  // void Calibration::BuildTrackCalc(Float_t x0, Float_t y0, Float_t x1, Float_t y1, TRACK* out)
  // lines 557 and following
  BeamDir.RotateY( (track->GetATA() - set->GretinaAngleA() )/1000 ); 
  BeamDir.Rotate(  (track->GetBTA() - set->GretinaAngleB() )/1000, -BeamDir.Cross(TVector3(0,1,0)) ); 
  
  //old not Gretina ATA and BTA
  //BeamDir.SetTheta(track->GetTheta());
  //BeamDir.SetPhi(track->GetPhi());
  double CosDop = cos(PosToTarget.Angle(BeamDir));

  double beta = set->TargetBeta();
  double gamma0 = 1/sqrt(1.0 - beta*beta);
  double gamma = gamma0 + (gamma0-1)*track->GetDTA()/100.; //logbook page 42 for derivation
  beta = sqrt( 1 - 1./(gamma*gamma) );
  return gamma*(1-beta*CosDop);

#endif
}

double HitCalc::DopplerCorrectionFactor(TVector3 PosToTarget, double beta, double z=0){
  TVector3 zpos;
  zpos.SetXYZ(0,0,z);
  zpos += PosToTarget;
  return (1-beta*cos(zpos.Theta()))/sqrt(1-beta*beta);
}

void GretinaCalc::DopplerCorrect(Settings* set, TRACK* track){
  for(vector<HitCalc*>::iterator hit=fhits.begin(); hit!=fhits.end(); hit++){
    (*hit)->DopplerCorrect(set,track);
  }
  for(vector<HitCalc*>::iterator hit=fhits_ab.begin(); hit!=fhits_ab.end(); hit++){
    (*hit)->DopplerCorrect(set,track);
  }
}
  
