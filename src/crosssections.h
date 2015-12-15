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
#include <iostream>
#include <stdlib.h>

using namespace std;

double co(Double_t E){
  //klein-nishina formula
  double e=E/511.; // E/m_e
  double op2e = 1+2*e; 
  double r = 2.8179e-15; // classical electron radius
  double z_ge = 32; // Z for germanium
  return 2*TMath::Pi()*r*r*z_ge* ( (1+e)/e/e * ( (2+2*e)/(op2e) - log(op2e)/e ) + log(op2e)/2/e - (op2e+e)/op2e/op2e) *1e28;
  //siegbahn page 51
}
double ph(Double_t E){
  double e=511./E;
  double a = 7.2973525698e-3; //alpha finestructure constant
  double r = 2.8179e-15; // classical electron radius
  double z_ge = 32; // Z for germanium
  // the 9./8 takes into a account the KL shell
  return 9./8*8./3.*TMath::Pi()*r*r *pow(z_ge,5)*pow(a,4)*pow(2,5./2.)*pow(e,7./2.)*1e28 *1./e;
  //siegbahn page 43
}
double pa(Double_t E){
  if(E<1022)
    return 0;
  double par[4];
  par[0] = 0.948261;
  par[1] = -0.1332e-3;
  par[2] = 0.15567e-6;
  par[3] = 0.792189;
  
  //fit to data E. Storm, H.I. Israel, Nucl Data Tables A7, 565 (1970)
  return log(par[0] + E*par[1] + E*E*par[2])*par[3];

}
