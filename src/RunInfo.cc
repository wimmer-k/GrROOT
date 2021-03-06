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
#include "RunInfo.hh"
void RunInfo::Print(){
  cout << "-------------------------------------------" << endl;
  cout << "Run # " << frunnr << endl;
  cout << "run time " << GetRunTime() << endl;
  cout << "Events " << fevents << endl;
  cout << "Entries: "<< endl;
  cout << "Total\t" << fentries << endl;
  cout << "IC\t" << fich << endl;
  cout << "SCINT\t" << fscint << endl;
  cout << "Hodo\t" << fhodo << endl;
  cout << "Greta\t" << fgreta << endl;
  cout << "Card29\t" << fcard29 << endl;
  cout << "Efficiencies " <<setprecision(5)<< endl;
  cout << "OBJ\t"    << fobjeff*100. <<" %" << endl;
  cout << "XFP\t"    << fxfpeff*100. <<" %" << endl;
  cout << "TOF\t"    << ftofeff*100. <<" %" << endl;
  cout << "CRDC 0\t" << fpadeff[0]*100. <<" %" <<  endl;
  cout << "CRDC 1\t" << fpadeff[1]*100. <<" %" << endl;
  cout << "TRACK\t"  << ftrackeff*100. <<" %" << endl;
  cout << "PPAC 0\t" << fppaceff[0]*100. <<" %" <<  endl;
  cout << "PPAC 1\t" << fppaceff[1]*100. <<" %" << endl;
  cout << "IITRACK\t"  << fiitrackeff*100. <<" %" << endl;
  cout << "Card 29\t"  << fcard29eff*100. <<" %" << endl;
  cout << "LifeTime " << GetLifeTime() <<endl;
  cout << "Rates " << endl;
  cout << "XFP\t" << (double)fintegral[25]/GetRunTime() << " /s " << endl;
  cout << "OBJ\t" << (double)fintegral[24]/GetRunTime() << " /s " << endl;
  cout << "S800\t" << (double)fintegral[4]/GetRunTime() << " /s " << endl;
  cout << "Raw\t" << (double)fintegral[9]/GetRunTime() << " /s " << endl;
  cout << "Live\t" << (double)fintegral[10]/GetRunTime() << " /s " << endl;
  cout << "Coinc\t" << (double)fintegral[5]/GetRunTime() << " /s " << endl;
  cout << "Second\t" << (double)fintegral[8]/GetRunTime() << " /s " << endl;
  cout << "-------------------------------------------" << endl;

}
