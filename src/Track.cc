#include <iostream>
#include <iomanip>
#include <string.h>
#include <sys/time.h>
#include <signal.h>


#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH2S.h"
#include "TH1S.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TMath.h"
#include "TCutG.h"
#include "TEnv.h"
#include "TKey.h"
#include "TDirectory.h"

#include "Tracking.hh"
#include "TrackSettings.hh"
#include "Calibration.hh"
#include "CommandLineInterface.hh"
#include "S800Calc.hh"
#include "GretinaCalc.hh"
#include "GretinaTrack.hh"

#define TRACKING true

using namespace TMath;
using namespace std;

bool signal_received = false;
void signalhandler(int sig){
  if (sig == SIGINT){
    signal_received = true;
  }
}

double get_time(){
  struct timeval t;
  gettimeofday(&t, NULL);
  double d = t.tv_sec + (double) t.tv_usec/1000000;
  return d;
}
int main(int argc, char* argv[]){
  signal(SIGINT,signalhandler);
  double time_start = get_time();
  vector<char*> InputFiles;
  char* SettingFile = NULL;
  char* OutputFile = NULL;
  char* TreeName = (char*)"ctr";
  int nmax =0;
  int vl =0;

  CommandLineInterface* interface = new CommandLineInterface();

  interface->Add("-i", "inputfiles", &InputFiles);
  interface->Add("-o", "outputfile", &OutputFile);
  interface->Add("-tn", "name of the tree", &TreeName);
  interface->Add("-s", "settingsfile", &SettingFile);
  interface->Add("-n", "nmax", &nmax);
  interface->Add("-v", "verbose", &vl);
  interface->CheckFlags(argc, argv);

  if(InputFiles.size() == 0 || OutputFile == NULL){
    cerr<<"You have to provide at least one input file and the output file!"<<endl;
    exit(1);
  }
  cout<<"input file(s):"<<endl;
  for(unsigned int i=0; i<InputFiles.size(); i++){
    cout<<InputFiles[i]<<endl;
  }
  cout<<"output file: "<<OutputFile<< endl;

  TChain* tr = new TChain(TreeName);
  for(unsigned int i=0; i<InputFiles.size(); i++){
    tr->Add(InputFiles[i]);
  }
  if(tr == NULL){
    cout << "could not find tree ctr in file " << endl;
    for(unsigned int i=0; i<InputFiles.size(); i++){
      cout<<InputFiles[i]<<endl;
    }
    return 3;
  }
  if(SettingFile == NULL){
    cout << "No settings file given " << endl;
    return 3;
  }
  cout << "creating outputfile " << endl;
  TFile* outfile = new TFile(OutputFile,"recreate");
  TrackSettings* set = new TrackSettings(SettingFile);
  Calibration *cal = new Calibration(set,0);
  TEnv *cutsset = new TEnv(SettingFile); 

  GretinaCalc* gr = new GretinaCalc;
  S800Calc* s800 = new S800Calc;
  tr->SetBranchAddress("gretinacalc",&gr);
  tr->SetBranchAddress("s800calc",&s800);


  Tracking *track = new Tracking(set);
  TTree* ttr = new TTree("ttr","Gretina tracked data");
  GretinaEvent *ge = new GretinaEvent;
  ttr->Branch("ge",&ge,320000);
  ttr->Branch("s800calc",&s800,320000);
  ttr->BranchRef();

  Double_t nentries = tr->GetEntries();
  cout << nentries << " entries in tree " << endl;

  if(outfile->IsZombie()){
    return 4;
  }
  double pprange[3];
  double gamrange[3];
  pprange[0] = cutsset->GetValue("PP.NBins",100);
  pprange[1] = cutsset->GetValue("PP.Low",13);
  pprange[2] = cutsset->GetValue("PP.High",14);
  gamrange[0] = cutsset->GetValue("Gam.NBins",4000);
  gamrange[1] = cutsset->GetValue("Gam.Low",0);
  gamrange[2] = cutsset->GetValue("Gam.High",4000);

  TH1F* hbeta = new TH1F("hbeta","hbeta",400,0.25,0.45);
  TH1F* hdta = new TH1F("hdta","hdta",200,-10,10);
  TH1F* hindex = new TH1F("hindex","hindex",20,-1.5,18.5);

  TH1F* hegam = new TH1F("hegam","hegam",gamrange[0],gamrange[1],gamrange[2]);
  TH1F* hegamdc = new TH1F("hegamdc","hegamdc",gamrange[0],gamrange[1],gamrange[2]);
  TH2F* hegamdc_dta = new TH2F("hegamdc_dta","hegamdc_dta",gamrange[0],gamrange[1],gamrange[2],100,-10,10);
  TH2F* hegamdc_theta = new TH2F("hegamdc_theta","hegamdc_theta",gamrange[0],gamrange[1],gamrange[2],100,0,3);
  TH2F* hegamdc_scatt = new TH2F("hegamdc_scatt","hegamdc_scatt",gamrange[0],gamrange[1],gamrange[2],100,0,3);
  TH2F* hegamdc_phi = new TH2F("hegamdc_phi","hegamdc_phi",gamrange[0],gamrange[1],gamrange[2],100,-7,0);
  TH2F* hegamdc_dphi = new TH2F("hegamdc_dphi","hegamdc_dphi",gamrange[0],gamrange[1],gamrange[2],100,0,7);
  TH2F* hegamegamdc = new TH2F("hegamegamdc","hegamegamdc",gamrange[0],gamrange[1],gamrange[2],gamrange[0],gamrange[1],gamrange[2]);

  TH1F* hegamAB = new TH1F("hegamAB","hegamAB",gamrange[0],gamrange[1],gamrange[2]);
  TH1F* hegamABdc = new TH1F("hegamABdc","hegamABdc",gamrange[0],gamrange[1],gamrange[2]);
  TH2F* hegamABdc_dta = new TH2F("hegamABdc_dta","hegamABdc_dta",gamrange[0],gamrange[1],gamrange[2],100,-10,10);
  TH2F* hegamABdc_theta = new TH2F("hegamABdc_theta","hegamABdc_theta",gamrange[0],gamrange[1],gamrange[2],100,0,3);
  TH2F* hegamABdc_scatt = new TH2F("hegamABdc_scatt","hegamABdc_scatt",gamrange[0],gamrange[1],gamrange[2],100,0,3);
  TH2F* hegamABdc_phi = new TH2F("hegamABdc_phi","hegamABdc_phi",gamrange[0],gamrange[1],gamrange[2],100,-7,0);
  TH2F* hegamABdc_dphi = new TH2F("hegamABdc_dphi","hegamABdc_dphi",gamrange[0],gamrange[1],gamrange[2],100,0,7);
  TH2F* hegamegamABdc = new TH2F("hegamegamABdc","hegamegamABdc",gamrange[0],gamrange[1],gamrange[2],gamrange[0],gamrange[1],gamrange[2]);

  TH1F* hegamTR = new TH1F("hegamTR","hegamTR",gamrange[0],gamrange[1],gamrange[2]);
  TH1F* hegamTRdc = new TH1F("hegamTRdc","hegamTRdc",gamrange[0],gamrange[1],gamrange[2]);
  TH1F* hegamTRdccut = new TH1F("hegamTRdccut","hegamTRdccut",gamrange[0],gamrange[1],gamrange[2]);
  TH2F* hegamTRdc_dta = new TH2F("hegamTRdc_dta","hegamTRdc_dta",gamrange[0],gamrange[1],gamrange[2],100,-10,10);
  TH2F* hegamTRdc_theta = new TH2F("hegamTRdc_theta","hegamTRdc_theta",gamrange[0],gamrange[1],gamrange[2],100,0,3);
  TH2F* hegamTRdc_scatt = new TH2F("hegamTRdc_scatt","hegamTRdc_scatt",gamrange[0],gamrange[1],gamrange[2],100,0,3);
  TH2F* hegamTRdc_phi = new TH2F("hegamTRdc_phi","hegamTRdc_phi",gamrange[0],gamrange[1],gamrange[2],100,-7,0);
  TH2F* hegamTRdc_dphi = new TH2F("hegamTRdc_dphi","hegamTRdc_dphi",gamrange[0],gamrange[1],gamrange[2],100,0,7);
  TH2F* hegamTR_FOM = new TH2F("hegamTR_FOM","hegamTR_FOM",2000,-1,1,gamrange[0],gamrange[1],gamrange[2]);
  TH2F* hegamTRdc_FOM = new TH2F("hegamTRdc_FOM","hegamTRdc_FOM",2000,-1,1,gamrange[0],gamrange[1],gamrange[2]);
  TH2F* hegamTRdc_FOMf = new TH2F("hegamTRdc_FOMf","hegamTRdc_FOMf",1000,0.9,1,gamrange[0],gamrange[1],gamrange[2]);
  TH2F* hegamTRdc_pp = new TH2F("hegamTRdc_pp","hegamTRdc_pp",gamrange[0],gamrange[1],gamrange[2],pprange[0],pprange[1],pprange[2]);
  TH2F* hegamTRdc_ppcut = new TH2F("hegamTRdc_ppcut","hegamTRdc_ppcut",gamrange[0],gamrange[1],gamrange[2],pprange[0],pprange[1],pprange[2]);
  TH2F* hegamegamTRdc = new TH2F("hegamegamTRdc","hegamegamTRdc",gamrange[0],gamrange[1],gamrange[2],gamrange[0],gamrange[1],gamrange[2]);
  TH2F* hegamegamTRdccut = new TH2F("hegamegamTRdccut","hegamegamTRdccut",gamrange[0],gamrange[1],gamrange[2],gamrange[0],gamrange[1],gamrange[2]);



  Int_t nbytes = 0;
  Int_t status;

  double beta0 = set->TargetBeta();
  double gamma0 = 1/sqrt(1.0 - beta0*beta0);
  double gamma, beta;

  TVector3 PosToTarget,BeamDir;
  HitCalc* hit;
  vector<HitCalc*> hitsAB;
  vector<HitCalc*> hitsRaw;

  int writetree = 1;
  writetree = cutsset->GetValue("WriteTree",1);
  if(writetree>0)
    cout << "writing tree " << endl;
  else
    cout << "no tree " << endl;
  int simulation = 0;
  simulation = cutsset->GetValue("Simulation",0);
  if(simulation>0)
    cout << "analyzing simulation " << endl;
  else
    cout << "analyzing data " << endl;

  cout << nentries << " entries in tree " << endl;

  if(nmax>0)
    nentries = nmax;
  for(int i=0; i<nentries;i++){
    if(signal_received){
      break;
    }
    if(vl>1)
      cout << "--------------------------"<< "getting entry " << i << "--------------------------"<< endl;
    status = tr->GetEvent(i);

    if(status == -1){
      cerr<<"Error occured, couldn't read entry "<<i<<" from tree "<<tr->GetName()<<" in file "<<tr->GetFile()->GetName()<<endl;
      return 5;
    }
    else if(status == 0){
      cerr<<"Error occured, entry "<<i<<" in tree "<<tr->GetName()<<" in file "<<tr->GetFile()->GetName()<<" doesn't exist"<<endl;
      return 6;
    }
    nbytes += status;

    TRACK *s800track = s800->GetTRACK();
    if(set->RedoMap())
      cal->BuildTrackCalc(s800->GetPAD(0)->GetX(),s800->GetPAD(0)->GetY(),s800->GetPAD(1)->GetX(),s800->GetPAD(1)->GetY(),
			  s800track);

    BeamDir = TVector3(0,0,1);
    //theta and phi do have the y -> -y rotation already in, ata and bta not
    // void Calibration::BuildTrackCalc(Float_t x0, Float_t y0, Float_t x1, Float_t y1, TRACK* out)
    // lines 557 and following
    BeamDir.RotateY( (s800->GetTRACK()->GetATA() - set->GretinaAngleA() )/1000 );
    BeamDir.Rotate(  (s800->GetTRACK()->GetBTA() - set->GretinaAngleB() )/1000, -BeamDir.Cross(TVector3(0,1,0)) );

    gamma = gamma0 + (gamma0-1)*s800->GetTRACK()->GetDTA()/100.; //logbook page 42 for derivation
    beta = sqrt( 1 - 1./(gamma*gamma) );
    hbeta->Fill(beta);
    hdta->Fill(s800->GetTRACK()->GetDTA());
    //loop over raw hits


    hitsAB.clear();
    hitsRaw.clear();

    int highestg=-1;
    double highesten =0;

    for(UShort_t g=0;g<gr->GetMult();g++){ // looping over gamma events
      hit = gr->GetHit(g);
      //PosToTarget = hit->GetPosition() - set->TargetPos();
      PosToTarget = hit->GetPosition();
      PosToTarget.SetZ(PosToTarget.Z() + set->TargetPos().Z()); 
      PosToTarget.SetX(PosToTarget.X() + set->TargetPos().X());       
      PosToTarget.SetY(PosToTarget.Y() - s800->GetTRACK()->GetYTA());
      double CosDop = cos(PosToTarget.Angle(BeamDir));
      hit->SetDCEnergy(hit->GetEnergy()*gamma*(1-beta*CosDop));
      hegam->Fill(hit->GetEnergy());
      hegamdc->Fill(hit->GetDCEnergy());
      hegamdc_dta->Fill(hit->GetDCEnergy(),s800->GetTRACK()->GetDTA());
      hegamdc_theta->Fill(hit->GetDCEnergy(),hit->GetPosition().Theta());
      hegamdc_scatt->Fill(hit->GetDCEnergy(),PosToTarget.Angle(BeamDir));
      hegamdc_phi->Fill(hit->GetDCEnergy(),s800->GetTRACK()->GetPhi());
      double dphi = s800->GetTRACK()->GetPhi() - PosToTarget.Phi();
      while(dphi<0)
	dphi+=2*TMath::Pi();
      while(dphi>2*TMath::Pi())
	dphi-=2*TMath::Pi();
      hegamdc_dphi->Fill(hit->GetDCEnergy(),dphi);
      hitsRaw.push_back(hit);
      if(hit->GetDCEnergy()>highesten){
	highesten = hit->GetDCEnergy();
	highestg = g;
      }	
    }
    if(gr->GetMult()>=2){
      for(UShort_t g=0;g<gr->GetMult();g++){ // looping over gamma events
	if(highestg>-1 && g!=highestg){
	  hegamegamdc->Fill(highesten,gr->GetHit(g)->GetDCEnergy());
	}
      }
    }

    highestg=-1;
    highesten =0;
    //loop over addback hits
    for(UShort_t g=0;g<gr->GetMultAB();g++){ // looping over gamma events
      hit = gr->GetHitAB(g);
      //PosToTarget = hit->GetPosition() - set->TargetPos();
      PosToTarget = hit->GetPosition();
      PosToTarget.SetZ(PosToTarget.Z() + set->TargetPos().Z()); 
      PosToTarget.SetX(PosToTarget.X() + set->TargetPos().X()); 
      PosToTarget.SetY(PosToTarget.Y() - s800->GetTRACK()->GetYTA());
      double CosDop = cos(PosToTarget.Angle(BeamDir));
      hit->SetDCEnergy(hit->GetEnergy()*gamma*(1-beta*CosDop));
      hegamAB->Fill(hit->GetEnergy());
      hegamABdc->Fill(hit->GetDCEnergy());
      hegamABdc_dta->Fill(hit->GetDCEnergy(),s800->GetTRACK()->GetDTA());
      hegamABdc_theta->Fill(hit->GetDCEnergy(),hit->GetPosition().Theta());
      hegamABdc_scatt->Fill(hit->GetDCEnergy(),PosToTarget.Angle(BeamDir));
      hegamABdc_phi->Fill(hit->GetDCEnergy(),s800->GetTRACK()->GetPhi());
      double dphi = s800->GetTRACK()->GetPhi() - PosToTarget.Phi();
      while(dphi<0)
	dphi+=TMath::Pi();
      while(dphi>2*TMath::Pi())
	dphi-=TMath::Pi();
      hegamABdc_dphi->Fill(hit->GetDCEnergy(),dphi);
      hitsAB.push_back(hit);
      if(hit->GetDCEnergy()>highesten){
	highesten = hit->GetDCEnergy();
	highestg = g;
      }	
    }
    if(gr->GetMult()>=2){
      for(UShort_t g=0;g<gr->GetMultAB();g++){ // looping over gamma events
	if(highestg>-1 && g!=highestg){
	  hegamegamABdc->Fill(highesten,gr->GetHitAB(g)->GetDCEnergy());
	}
      }
    }

    if(TRACKING){
      track->SetGretina(gr);
      track->SortInClusters();
      ge = track->GetEvent();
      highestg=-1;
      highesten =0;
      int highestgcut=-1;
      double highestencut =0;
      for(UShort_t g=0;g<ge->GetMult();g++){ // looping over gamma events
	if(ge->GetTrack(g)->GetMult()==0)
	  continue;
	hit = ge->GetTrack(g)->GetHits()[0];
	//PosToTarget = hit->GetPosition() - set->TargetPos();
	PosToTarget = hit->GetPosition();
	PosToTarget.SetZ(PosToTarget.Z() + set->TargetPos().Z()); 
	PosToTarget.SetX(PosToTarget.X() + set->TargetPos().X()); 
	PosToTarget.SetY(PosToTarget.Y() - s800->GetTRACK()->GetYTA());
	double CosDop = cos(PosToTarget.Angle(BeamDir));
	double DCen = ge->GetTrack(g)->GetEsum()*gamma*(1-beta*CosDop);
	hegamTR->Fill(ge->GetTrack(g)->GetEsum());
	hegamTR_FOM->Fill(ge->GetTrack(g)->GetFOM(),ge->GetTrack(g)->GetEsum());
	ge->GetTrack(g)->SetDCEnergy(DCen);
	hegamTRdc->Fill(DCen);
	hindex->Fill(hit->GetIndex());
	hegamTRdc_FOM->Fill(ge->GetTrack(g)->GetFOM(),DCen);
	hegamTRdc_FOMf->Fill(ge->GetTrack(g)->GetFOM(),DCen);
	hegamTRdc_dta->Fill(DCen,s800->GetTRACK()->GetDTA());
	hegamTRdc_theta->Fill(DCen,hit->GetPosition().Theta());
	hegamTRdc_scatt->Fill(DCen,PosToTarget.Angle(BeamDir));
	hegamTRdc_phi->Fill(DCen,s800->GetTRACK()->GetPhi());
	hegamTRdc_pp->Fill(DCen,s800->GetTRACK()->GetPpar());
	if(ge->GetTrack(g)->GetFOM()>0.98 && ge->GetTrack(g)->GetFOM()<1.1){
	  hegamTRdc_ppcut->Fill(DCen,s800->GetTRACK()->GetPpar());
	  hegamTRdccut->Fill(DCen);
	}
	double dphi = s800->GetTRACK()->GetPhi() - PosToTarget.Phi();
	while(dphi<0)
	  dphi+=TMath::Pi();
	while(dphi>2*TMath::Pi())
	  dphi-=TMath::Pi();
	hegamTRdc_dphi->Fill(DCen,dphi);
	if(DCen>highesten){
	  highesten = DCen;
	  highestg = g;
	}	
	if(ge->GetTrack(g)->GetFOM()>0.98 && ge->GetTrack(g)->GetFOM()<1.1 && DCen>highestencut){
	  highestencut = DCen;
	  highestgcut = g;
	}	
      }//gammas
      if(ge->GetMult()>=2){
	for(UShort_t g=0;g<ge->GetMult();g++){ // looping over gamma events
	  if(highestg>-1 && g!=highestg){
	    hegamegamTRdc->Fill(highesten,ge->GetTrack(g)->GetDCEnergy());
	  }
	  if(highestgcut>-1 && g!=highestgcut && ge->GetTrack(g)->GetFOM()>0.98 && ge->GetTrack(g)->GetFOM()<1.1){
	    hegamegamTRdccut->Fill(highestencut,ge->GetTrack(g)->GetDCEnergy());
	  }
	}
      }
    }//tracking ON
    ge->SetHitsRaw(hitsRaw);
    ge->SetHitsAB(hitsAB);
    
    if(writetree>0)
      ttr->Fill();

    
    if(i%1000 == 0){
      double time_end = get_time();
      cout << setw(5) << setiosflags(ios::fixed) << setprecision(1) << (100.*i)/nentries <<
	" % done\t" << (Float_t)i/(time_end - time_start) << " events/s " <<
	(nentries-i)*(time_end - time_start)/(Float_t)i << "s to go \r" << flush;
    }
  }
  cout << endl;
  outfile->cd();
  if(writetree>0){
    cout << "writing tree " << endl;
    ttr->Write("",TObject::kOverwrite);
  }

  hbeta->Write("",TObject::kOverwrite);
  hdta->Write("",TObject::kOverwrite);
  hegam->Write("",TObject::kOverwrite);
  hegamdc->Write("",TObject::kOverwrite);
  hegamdc_dta->Write("",TObject::kOverwrite);
  hegamdc_theta->Write("",TObject::kOverwrite);
  hegamdc_scatt->Write("",TObject::kOverwrite);
  hegamdc_phi->Write("",TObject::kOverwrite);
  hegamdc_dphi->Write("",TObject::kOverwrite);
  hegamegamdc->Write("",TObject::kOverwrite);

  hindex->Write("",TObject::kOverwrite);

  hegamAB->Write("",TObject::kOverwrite);
  hegamABdc->Write("",TObject::kOverwrite);
  hegamABdc_dta->Write("",TObject::kOverwrite);
  hegamABdc_theta->Write("",TObject::kOverwrite);
  hegamABdc_scatt->Write("",TObject::kOverwrite);
  hegamABdc_phi->Write("",TObject::kOverwrite);
  hegamABdc_dphi->Write("",TObject::kOverwrite);
  hegamegamABdc->Write("",TObject::kOverwrite);

  hegamTR->Write("",TObject::kOverwrite);
  hegamTR_FOM->Write("",TObject::kOverwrite);
  hegamTRdc->Write("",TObject::kOverwrite);
  hegamTRdccut->Write("",TObject::kOverwrite);
  hegamTRdc_FOM->Write("",TObject::kOverwrite);
  hegamTRdc_FOMf->Write("",TObject::kOverwrite);
  hegamTRdc_dta->Write("",TObject::kOverwrite);
  hegamTRdc_theta->Write("",TObject::kOverwrite);
  hegamTRdc_scatt->Write("",TObject::kOverwrite);
  hegamTRdc_phi->Write("",TObject::kOverwrite);
  hegamTRdc_dphi->Write("",TObject::kOverwrite);
  hegamTRdc_pp->Write("",TObject::kOverwrite);
  hegamTRdc_ppcut->Write("",TObject::kOverwrite);
  hegamegamTRdc->Write("",TObject::kOverwrite);
  hegamegamTRdccut->Write("",TObject::kOverwrite);

  outfile->Close();
  delete tr;

  double time_end = get_time();
  cout << "Run time " << time_end - time_start << " s." << endl;

  return 0;
}
