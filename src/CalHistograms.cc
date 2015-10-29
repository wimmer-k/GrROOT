#include <iostream>
#include <iomanip>
#include <string.h>
#include <sys/time.h>

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

#include "CommandLineInterface.hh"
#include "S800Calc.hh"
#include "GretinaCalc.hh"
#include "Mode3Calc.hh"
#include "Scaler.hh"
#include "CalHistograms.hh"
using namespace TMath;
using namespace std;

//not nice but works at least
static TCutG* TimeCut;
static bool foundTCut = false;

void CalHistograms::FillHistograms(GretinaCalc* gr, S800Calc* s800, Mode3Calc* m3c){

  //bool hass800 = s800->GetTimeS800()!=0;
  bool hasgret = gr->GetMult()>0;

  //Declare the histograms here.
  static vector<TCutG*> InPartCut;
  static vector<vector<TCutG*> > OutPartCut;
 
  static PAD* pad[2];
  static PPAC* ppac[2];
  static IC* ich;
  static SCINT* scint;
  static TOF* tof;
  static TRACK* track;
  static IITRACK* iitrack;
  static HODO* hodo;

  static bool foundCuts = false;
  static bool NoInCuts = false;
  if (!foundCuts){
    //Read in the cuts file for incoming and outgoing particle ID
    if(fhasfile){
      cout << "Cuts were created for ";
      if (ftac&(1<<0))
	cout << "TDC";
      else
	cout << "TAC";
      cout << " data and are applied to the ";
      if (ftac&(1<<1))
	cout << "XFP";
      else
	cout << "OBJ";
      cout << " scintillator" << endl;
      

      //Remember the current directory so that we can revert to it later.
      TDirectory* outfile = gDirectory;

      char* Name = NULL;
      char* Name2 = NULL;
      TFile *cFile = new TFile(fcutfile);
      TIter next(cFile->GetListOfKeys());
      TKey *key;
      while((key=(TKey*)next())){
	if(strcmp(key->GetClassName(),"TCutG") == 0){
	  Name = (char*)key->GetName();
	  if(strstr(Name,"in") && !strstr(Name,"out")){
	    cout << "incoming cut found "<<Name << endl;
	    InPartCut.push_back((TCutG*)cFile->Get(Name));
	  }
	  if(strstr(Name,"tgam")){
	    cout << "timing cut found "<<Name << endl;
	    TimeCut = (TCutG*)cFile->Get(Name);
	    foundTCut = true;
	  }
	}
      }
      TIter next2(cFile->GetListOfKeys());
      if(InPartCut.size()>0){
	OutPartCut.resize(InPartCut.size());
	while((key=(TKey*)next2())){
	  if(strcmp(key->GetClassName(),"TCutG") == 0){
	    Name = (char*)key->GetName();
	    if(strstr(Name,"in") && strstr(Name,"out")){
	      for(unsigned short i=0;i<InPartCut.size();i++){
		Name2 = (char*)InPartCut[i]->GetName();
		if(strstr(Name,strstr(Name2,Name2))){
		  OutPartCut[i].push_back((TCutG*)cFile->Get(Name));
		  cout << "outgoing cut found "<<Name << endl;
		}
	      }
	    }
	  }
	}
      }
      else{
	OutPartCut.resize(1);
	while((key=(TKey*)next2())){
	  if(strcmp(key->GetClassName(),"TCutG") == 0){
	    Name = (char*)key->GetName();
	    if(!strstr(Name,"in") && strstr(Name,"out")){
	      OutPartCut[0].push_back((TCutG*)cFile->Get(Name));
	      cout << "outgoing cut found (no incoming!) "<<Name << endl;
	      NoInCuts = true;
	    }
	  }
	}
      }
      cFile->Close();
      outfile->cd();
 
      foundCuts = true;
      cout << "Cuts found" << endl;
    }
  }

  //-------------------------------------------------------------------------
  //*************************************************************************
  //Fill the histograms here.
  //*************************************************************************
  //-------------------------------------------------------------------------
 
  int hp = s800->GetRegistr();
  FillI("registr",16,0,16,hp);
  for(int j=0;j<16;j++){
    if(hp & (1<<j))
      FillI("trigbit",16,0,16,j);
  }

  ich = s800->GetIC();
  tof = s800->GetTOF();
  track = s800->GetTRACK();
  iitrack = s800->GetIITRACK();
  hodo = s800->GetHODO();
  for(UShort_t p=0; p<2;p++){
    pad[p] = s800->GetPAD(p);
    ppac[p] = s800->GetPPAC(p);
  }

  FillHistogramsNoGate(gr,s800,m3c);
  FillHistogramsGateIn(gr,s800,m3c,"all");
  FillHistogramsGateOut(gr,s800,m3c,"all");
  for(UShort_t in=0;in<InPartCut.size();in++){ // loop over incoming cuts
    if( (ftac%2==1 && InPartCut[in]->IsInside(tof->GetOBJ(),tof->GetXFP()))
	|| (ftac%2==0 && InPartCut[in]->IsInside(tof->GetTACOBJ(),tof->GetTACXFP())) ){
      const char* inname = InPartCut[in]->GetName();

      FillHistogramsGateIn(gr,s800,m3c,inname);
      if(hasgret){
	FillHistogramsGateIn(gr,s800,m3c,Form("%s_coinc",inname));
      }
      FillHistogramsGateOut(gr,s800,m3c,inname);

      for(UShort_t ou=0;ou<OutPartCut[in].size();ou++){ // loop over PID cuts
	if((ftac == 1 && OutPartCut[in][ou]->IsInside(tof->GetOBJC(),ich->GetDE()))
	   ||(ftac == 0 && OutPartCut[in][ou]->IsInside(tof->GetTACOBJC(),ich->GetDE()))
	   ||(ftac == 3 && OutPartCut[in][ou]->IsInside(tof->GetXFPC(),ich->GetDE()))
	   ||(ftac == 2 && OutPartCut[in][ou]->IsInside(tof->GetTACXFPC(),ich->GetDE()))){
	  const char* outname = OutPartCut[in][ou]->GetName();

	  FillHistogramsGateOut(gr,s800,m3c,outname);
	  if (hasgret){
	    FillHistogramsGateOut(gr,s800,m3c,Form("%s_coinc",outname));
	  }

	}

      }
    }
  }
  if(NoInCuts){
    // for pure beams, i.e. no incoming cut required or possible
    for(UShort_t ou=0;ou<OutPartCut[0].size();ou++){ // loop over PID cuts
      if((ftac == 1 && OutPartCut[0][ou]->IsInside(tof->GetOBJC(),ich->GetDE()))
	 ||(ftac == 0 && OutPartCut[0][ou]->IsInside(tof->GetTACOBJC(),ich->GetDE()))
	 ||(ftac == 3 && OutPartCut[0][ou]->IsInside(tof->GetXFPC(),ich->GetDE()))
	 ||(ftac == 2 && OutPartCut[0][ou]->IsInside(tof->GetTACXFPC(),ich->GetDE()))){
	const char* outname = OutPartCut[0][ou]->GetName();
	
	FillHistogramsGateOnlyOut(gr,s800,m3c,outname);
	
      }
    }
  }
  
}
void CalHistograms::FillHistogramsNoGate(GretinaCalc* gr, S800Calc* s800, Mode3Calc* m3c){

  PAD* pad[2];
  PPAC* ppac[2];
  for(UShort_t p=0; p<2;p++){
    pad[p] = s800->GetPAD(p);
    ppac[p] = s800->GetPPAC(p);
  }
  IC* ich = s800->GetIC();
  TOF* tof = s800->GetTOF();
  TRACK* track = s800->GetTRACK();
  IITRACK* iitrack = s800->GetIITRACK();
  HODO* hodo = s800->GetHODO();
  bool hasgret = gr->GetMult()>0;

  // Double_t sum =0;
  // for(UShort_t c=0; c<hodo->GetEnergy()->size();c++){
  //   Short_t ch = hodo->GetChannel()->at(c);
  //   Fill(Form("hodo_%d",ch),1000,0,4000,hodo->GetEnergy()->at(c));
  //   Fill("hodo_vs_ch",32,0,32,ch,400,0,4000,hodo->GetEnergy()->at(c));
  //   Fill("hodo_all",1000,0,4000,hodo->GetEnergy()->at(c));
  //   Fill(Form("hodotime_%d",ch),4000,0,4000,hodo->GetTime());
  //   Fill("hodotime_vs_ch",32,0,32,ch,1000,0,4000,hodo->GetTime());
  //   Fill("hodotime_all",4000,0,4000,hodo->GetTime());
  //   if(hodo->GetEnergy()->at(c)>50)
  //     sum+=hodo->GetEnergy()->at(c);
  // }// hodo energy
  // Fill("hodo_sum",1000,0,4000,sum);

  if(fCal==0){
    Fill("obj_vs_time",
	 fnentries/10000+1,0,(fnentries/10000+1)*10000,fentry,
	 fobj_range[0]/2,fobj_range[1],fobj_range[2],tof->GetOBJ());
    Fill("xfp_vs_time",
	 fnentries/10000+1,0,(fnentries/10000+1)*10000,fentry,
	 fxfp_range[0]/2,fxfp_range[1],fxfp_range[2],tof->GetXFP());
    Fill("objtac_vs_time",
	 fnentries/10000+1,0,(fnentries/10000+1)*10000,fentry,
	 ftacobj_range[0]/2,ftacobj_range[1],ftacobj_range[2],tof->GetTACOBJ());
    Fill("xfptac_vs_time",
	 fnentries/10000+1,0,(fnentries/10000+1)*10000,fentry,
	 ftacxfp_range[0]/2,ftacxfp_range[1],ftacxfp_range[2],tof->GetTACXFP());
    Fill("IC_vs_time",
	 fnentries/10000+1,0,(fnentries/10000+1)*10000,fentry,
	 fIC_range[0]/2,fIC_range[1],fIC_range[2],ich->GetDE());
  }



  for(UShort_t p=0; p<2;p++){
    Fill(Form("ICde_vs_x_%d",p),
	 600,-300,300,pad[p]->GetX(),
	 fIC_range[0],fIC_range[1],fIC_range[2],ich->GetDE());
    Fill(Form("ICde_vs_y_%d",p),
	 200,-100,100,pad[p]->GetY(),
	 fIC_range[0],fIC_range[1],fIC_range[2],ich->GetDE());
    Fill(Form("ICsum_vs_x_%d",p),
	 600,-300,300,pad[p]->GetX(),
	 fIC_range[0],fIC_range[1],fIC_range[2],ich->GetSum());
    Fill(Form("ICsum_vs_y_%d",p),
	 200,-100,100,pad[p]->GetY(),
	 fIC_range[0],fIC_range[1],fIC_range[2],ich->GetSum());
  }

  if(hasgret){
    Fill("xfp_vs_obj_coinc",
	 fobj_range[0],fobj_range[1],fobj_range[2],tof->GetOBJ(),
	 fxfp_range[0],fxfp_range[1],fxfp_range[2],tof->GetXFP());
  }


  Fill("xfp_vs_obj",
       fobj_range[0],fobj_range[1],fobj_range[2],tof->GetOBJ(),
       fxfp_range[0],fxfp_range[1],fxfp_range[2],tof->GetXFP());
  Fill("x_vs_obj",
       fobj_range[0],fobj_range[1],fobj_range[2],tof->GetOBJ(),
       600,-300,300,pad[0]->GetX());
  Fill("afp_vs_obj",
       fobj_range[0],fobj_range[1],fobj_range[2],tof->GetOBJ(),
       100,-100,100,track->GetAFP());
  Fill("x_vs_objC",
       fobjC_range[0],fobjC_range[1],fobjC_range[2],tof->GetOBJC(),
       600,-300,300,pad[0]->GetX());
  Fill("afp_vs_objC",
       fobjC_range[0],fobjC_range[1],fobjC_range[2],tof->GetOBJC(),
       100,-50,50,track->GetAFP());

  Fill("x_vs_xfp",
       fxfp_range[0],fxfp_range[1],fxfp_range[2],tof->GetXFP(),
       600,-300,300,pad[0]->GetX());
  Fill("afp_vs_xfp",
       fxfp_range[0],fxfp_range[1],fxfp_range[2],tof->GetXFP(),
       100,-50,50,track->GetAFP());
  Fill("x_vs_xfpC",
       fxfpC_range[0],fxfpC_range[1],fxfpC_range[2],tof->GetXFPC(),
       600,-300,300,pad[0]->GetX());
  Fill("afp_vs_xfpC",
       fxfpC_range[0],fxfpC_range[1],fxfpC_range[2],tof->GetXFPC(),
       100,-50,50,track->GetAFP());

  Fill("xfp_vs_objtac",
       ftacobj_range[0],ftacobj_range[1],ftacobj_range[2],tof->GetTACOBJ(),
       ftacxfp_range[0],ftacxfp_range[1],ftacxfp_range[2],tof->GetTACXFP());
  Fill("x_vs_objtac",
       ftacobj_range[0],ftacobj_range[1],ftacobj_range[2],tof->GetTACOBJ(),
       600,-300,300,pad[0]->GetX());
  Fill("afp_vs_objtac",
       ftacobj_range[0],ftacobj_range[1],ftacobj_range[2],tof->GetTACOBJ(),
       100,-50,50,track->GetAFP());
  Fill("x_vs_objtacC",
       ftacobjC_range[0],ftacobjC_range[1],ftacobjC_range[2],tof->GetTACOBJC(),
       600,-300,300,pad[0]->GetX());
  Fill("afp_vs_objtacC",
       ftacobjC_range[0],ftacobjC_range[1],ftacobjC_range[2],tof->GetTACOBJC(),
       100,-50,50,track->GetAFP());

  Fill("x_vs_xfptac",
       ftacxfp_range[0],ftacxfp_range[1],ftacxfp_range[2],tof->GetTACXFP(),
       600,-300,300,pad[0]->GetX());
  Fill("afp_vs_xfptac",
       ftacxfp_range[0],ftacxfp_range[1],ftacxfp_range[2],tof->GetTACXFP(),
       100,-50,50,track->GetAFP());
  Fill("x_vs_xfptacC",
       ftacxfpC_range[0],ftacxfpC_range[1],ftacxfpC_range[2],tof->GetTACXFPC(),
       600,-300,300,pad[0]->GetX());
  Fill("afp_vs_xfptacC",
       ftacxfpC_range[0],ftacxfpC_range[1],ftacxfpC_range[2],tof->GetTACXFPC(),
       100,-50,50,track->GetAFP());
}


void CalHistograms::FillHistogramsGateIn(GretinaCalc* gr, S800Calc* s800, Mode3Calc* m3c, const char* inname){
  PAD* pad[2];
  PPAC* ppac[2];
  for(UShort_t p=0; p<2;p++){
    pad[p] = s800->GetPAD(p);
    ppac[p] = s800->GetPPAC(p);
  }
  IC* ich = s800->GetIC();
  TOF* tof = s800->GetTOF();
  TRACK* track = s800->GetTRACK();
  IITRACK* iitrack = s800->GetIITRACK();
  HODO* hodo = s800->GetHODO();
  bool hasgret = gr->GetMult()>0;

  Fill("hcrdcrawrange0xy",
       1200,0,600, s800->GetPAD(0)->GetX(),
       4096,0,4096,s800->GetPAD(0)->GetY());
  Fill("hcrdcrawrange1xy",
       1200,0,600, s800->GetPAD(1)->GetX(),
       4096,0,4096,s800->GetPAD(1)->GetY());
  Fill("hcrdc0xy",
       600,-300,300,s800->GetPAD(0)->GetX(),
       600,-300,300,s800->GetPAD(0)->GetY());
  Fill("hcrdc1xy",
       600,-300,300,s800->GetPAD(1)->GetX(),
       600,-300,300,s800->GetPAD(1)->GetY());


  Fill(Form("ICde_vs_obj_%s",inname),
       fobj_range[0],fobj_range[1],fobj_range[2],tof->GetOBJ(),
       fIC_range[0],fIC_range[1],fIC_range[2],ich->GetDE());
  Fill(Form("ICde_vs_objc_%s",inname),
       fobjC_range[0],fobjC_range[1],fobjC_range[2],tof->GetOBJC(),
       fIC_range[0],fIC_range[1],fIC_range[2],ich->GetDE());
  if (hasgret){
    Fill(Form("ICde_vs_objc_coinc_%s",inname),
	 fobjC_range[0],fobjC_range[1],fobjC_range[2],tof->GetOBJC(),
	 fIC_range[0],fIC_range[1],fIC_range[2],ich->GetDE());
  }
  Fill(Form("ICde_vs_objtacc_%s",inname),
       ftacobjC_range[0],ftacobjC_range[1],ftacobjC_range[2],tof->GetTACOBJC(),
       fIC_range[0],fIC_range[1],fIC_range[2],ich->GetDE());
  Fill(Form("ICde_vs_xfpc_%s",inname),
       fxfpC_range[0],fxfpC_range[1],fxfpC_range[2],tof->GetXFPC(),
       fIC_range[0],fIC_range[1],fIC_range[2],ich->GetDE());
  Fill(Form("ICde_vs_xfptacc_%s",inname),
       ftacxfpC_range[0],ftacxfpC_range[1],ftacxfpC_range[2],tof->GetTACXFPC(),
       fIC_range[0],fIC_range[1],fIC_range[2],ich->GetDE());

  Fill(Form("ICsum_vs_objc_%s",inname),
       fobjC_range[0],fobjC_range[1],fobjC_range[2],tof->GetOBJC(),
       fIC_range[0],fIC_range[1],fIC_range[2],ich->GetSum());
  Fill(Form("ICsum_vs_objtacc_%s",inname),
       ftacobjC_range[0],ftacobjC_range[1],ftacobjC_range[2],tof->GetTACOBJC(),
       fIC_range[0],fIC_range[1],fIC_range[2],ich->GetSum());
  Fill(Form("ICsum_vs_xfpc_%s",inname),
       fxfpC_range[0],fxfpC_range[1],fxfpC_range[2],tof->GetXFPC(),
       fIC_range[0],fIC_range[1],fIC_range[2],ich->GetSum());
  Fill(Form("ICsum_vs_xfptacc_%s",inname),
       ftacxfpC_range[0],ftacxfpC_range[1],ftacxfpC_range[2],tof->GetTACXFPC(),
       fIC_range[0],fIC_range[1],fIC_range[2],ich->GetSum());

  for(UShort_t p=0; p<2;p++){
    Fill(Form("ICde_vs_x%d_%s",p,inname),
	 600,-300,300,pad[p]->GetX(),
	 fIC_range[0],fIC_range[1],fIC_range[2],ich->GetDE());
    Fill(Form("ICde_vs_y%d_%s",p,inname),
	 200,-100,100,pad[p]->GetY(),
	 fIC_range[0],fIC_range[1],fIC_range[2],ich->GetDE());
    Fill(Form("ICsum_vs_x%d_%s",p,inname),
	 600,-300,300,pad[p]->GetX(),
	 fIC_range[0],fIC_range[1],fIC_range[2],ich->GetSum());
    Fill(Form("ICsum_vs_y%d_%s",p,inname),
	 200,-100,100,pad[p]->GetY(),
	 fIC_range[0],fIC_range[1],fIC_range[2],ich->GetSum());
  }
}


void CalHistograms::FillHistogramsGateOut(GretinaCalc* gr, S800Calc* s800, Mode3Calc* m3c, const char* outname){
  TString oname = outname;

  PAD* pad[2];
  PPAC* ppac[2];
  for(UShort_t p=0; p<2;p++){
    pad[p] = s800->GetPAD(p);
    ppac[p] = s800->GetPPAC(p);
  }
  IC* ich = s800->GetIC();
  TOF* tof = s800->GetTOF();
  TRACK* track = s800->GetTRACK();
  IITRACK* iitrack = s800->GetIITRACK();
  HODO* hodo = s800->GetHODO();
  bool hasgret = gr->GetMult()>0;

  for (UShort_t g=0; g<gr->GetMult(); g++){
    HitCalc* hit = gr->GetHit(g);
    float energy = hit->GetEnergy();
    float energy_dc = hit->GetDCEnergy();
    Fill(Form("egam_%s",oname.Data()),
	 8000,0,8000,energy);
    Fill(Form("egam_scatter_%s",oname.Data()),
	 500,0,0.5,track->GetTheta(),
	 8000,0,8000,energy);
    Fill(Form("egam_tgam_%s",oname.Data()),
	 8000,0,8000,energy,
	 400,-200,200,hit->GetTime());
    Fill(Form("tgam_egam_%s",oname.Data()),
	 400,-200,200,hit->GetTime(),
	 8000,0,8000,energy);
    Fill(Form("egam_summary_%s",oname.Data()),
	 28,-0.5,27.5,4*fSett->Hole2Det(hit->GetHole())+hit->GetCrystal(),
	 8000,0,8000,energy);
    Fill(Form("egamdc_%s",oname.Data()),
	 8000,0,8000,energy_dc);
    Fill(Form("egamdc_scatter_%s",oname.Data()),
	 500,0,0.5,track->GetTheta(),
	 8000,0,8000,energy_dc);
    Fill(Form("egamdc_grettheta_%s",oname.Data()),
	 8000,0,8000,energy_dc,
	 180,0,3.14159,hit->GetPosition().Theta());

    Fill(Form("egamdc_ppar_%s",oname.Data()),
	 8000,0,8000,energy_dc,
	 500,8,20,track->GetPpar());
    if(foundTCut && TimeCut->IsInside(hit->GetTime(),energy)){
      Fill(Form("egam_%s_tcut",oname.Data()),
	   8000,0,8000,energy);
      Fill(Form("egam_scatter_%s_tcut",oname.Data()),
	   500,0,0.5,track->GetTheta(),
	   8000,0,8000,energy);
      Fill(Form("egam_summary_%s_tcut",oname.Data()),
	   28,-0.5,27.5,4*fSett->Hole2Det(hit->GetHole())+hit->GetCrystal(),
	   8000,0,8000,energy);
      Fill(Form("egamdc_%s_tcut",oname.Data()),
	   8000,0,8000,energy_dc);
      Fill(Form("egamdc_scatter_%s_tcut",oname.Data()),
	   500,0,0.5,track->GetTheta(),
	   8000,0,8000,energy_dc);
      Fill(Form("egamdc_grettheta_%s_tcut",oname.Data()),
	   8000,0,8000,energy_dc,
	   180,0,3.14159,hit->GetPosition().Theta());
      Fill(Form("egamdc_ppar_%s_tcut",oname.Data()),
	   8000,0,8000,energy_dc,
	   500,8,20,track->GetPpar());
    }
    int q = fSett->Hole2Det(hit->GetHole())+1;
    if((q==1 || q==7) && hit->GetCrystal()==1){
      Fill(Form("egamdc_BW_%s",oname.Data()),
	   2000,0,4000,hit->GetDCEnergy(0.3));
    }
    else if((q==2 || q==6 || q==4 || q==5 ) && hit->GetCrystal()==2){
      Fill(Form("egamdc_FW_%s",oname.Data()),
	   2000,0,4000,hit->GetDCEnergy(0.3));
    }
    Fill(Form("egamdc_s800phi_%s",oname.Data()),
	 8000,0,8000,energy_dc,
	 180,-2*3.14159,0,s800->GetTRACK()->GetPhi());
    Fill(Form("egamdc_gretphi_%s",oname.Data()),
	 8000,0,8000,energy_dc,
	 180,-3.14159,3.14159,hit->GetPosition().Phi());
  }

  for(UShort_t g=0;g<gr->GetMultAB();g++){
    HitCalc* hit = gr->GetHitAB(g);
    float energy = hit->GetEnergy();
    float energy_dc = hit->GetDCEnergy();
    Fill(Form("egamAB_%s",oname.Data()),
	 8000,0,8000,energy);
    Fill(Form("egamABdc_%s",oname.Data()),
	 8000,0,8000,energy_dc);
    Fill(Form("egamABdc_ppar_%s",oname.Data()),
	 8000,0,8000,energy_dc,
	 500,8,20,track->GetPpar());
    if(foundTCut && TimeCut->IsInside(hit->GetTime(),energy)){
      Fill(Form("egamAB_%s_tcut",oname.Data()),
	   8000,0,8000,energy);
      Fill(Form("egamABdc_%s_tcut",oname.Data()),
	   8000,0,8000,energy_dc);
      Fill(Form("egamABdc_ppar_%s_tcut",oname.Data()),
	   8000,0,8000,energy_dc,
	   500,8,20,track->GetPpar());
    }
  }

  Fill(Form("ata_%s",oname.Data()),
       200,-100,100,track->GetATA());
  Fill(Form("bta_%s",oname.Data()),
       200,-100,100,track->GetBTA());
  Fill(Form("ata_bta_%s",oname.Data()),
       200,-100,100,track->GetATA(),
       200,-100,100,track->GetBTA());
  Fill(Form("scatter_%s",oname.Data()),
       500,0,0.5,track->GetTheta());

  // Various gamma-gamma spectra with no additional gates.
  {
    int highestg=-1;
    double highesten =0;
    for(UShort_t g=0;g<gr->GetMult();g++){
      if(gr->GetHit(g)->GetEnergy()<1 || gr->GetHit(g)->GetEnergy()>6000)
	continue;
      if(gr->GetHit(g)->GetDCEnergy()>highesten){
	highesten = gr->GetHit(g)->GetDCEnergy();
	highestg =g;
      }
    }
    for(UShort_t g=0;g<gr->GetMult();g++){
      if(gr->GetHit(g)->GetEnergy()<1 || gr->GetHit(g)->GetEnergy()>6000)
	continue;
      if(highestg>-1 && g!=highestg){
	Fill(Form("egamegamdc_fold_%s",oname.Data()),
	     8000,0,8000,highesten,
	     8000,0,8000,gr->GetHit(g)->GetDCEnergy());
      }
    }
    highestg=-1;
    highesten =0;
    for(UShort_t g=0;g<gr->GetMultAB();g++){
      if(gr->GetHitAB(g)->GetEnergy()<1 || gr->GetHitAB(g)->GetEnergy()>6000)
	continue;
      if(gr->GetHitAB(g)->GetDCEnergy()>highesten){
	highesten = gr->GetHitAB(g)->GetDCEnergy();
	highestg =g;
      }
    }
    for(UShort_t g=0;g<gr->GetMultAB();g++){
      if(gr->GetHitAB(g)->GetEnergy()<1 || gr->GetHitAB(g)->GetEnergy()>6000)
	continue;
      if(highestg>-1 && g!=highestg){
	Fill(Form("egamegamABdc_fold_%s",oname.Data()),
	     8000,0,8000,highesten,
	     8000,0,8000,gr->GetHitAB(g)->GetDCEnergy());
      }
    }
  }
  for(UShort_t g1=0;g1<gr->GetMult();g1++){
    for(UShort_t g2=g1+1;g2<gr->GetMult();g2++){
      Fill(Form("egamegamdc_%s",oname.Data()),
	   3500,0,3500,gr->GetHit(g1)->GetDCEnergy(),
	   3500,0,3500,gr->GetHit(g2)->GetDCEnergy());
      Fill(Form("egamegamdc_sym_%s",oname.Data()),
	   3500,0,3500,gr->GetHit(g1)->GetDCEnergy(),
	   3500,0,3500,gr->GetHit(g2)->GetDCEnergy());
      Fill(Form("egamegamdc_sym_%s",oname.Data()),
	   3500,0,3500,gr->GetHit(g2)->GetDCEnergy(),
	   3500,0,3500,gr->GetHit(g1)->GetDCEnergy());
    }
  }
  for(UShort_t g1=0;g1<gr->GetMultAB();g1++){
    for(UShort_t g2=g1+1;g2<gr->GetMultAB();g2++){
      Fill(Form("egamegamABdc_%s",oname.Data()),
	   3500,0,3500,gr->GetHitAB(g1)->GetDCEnergy(),
	   3500,0,3500,gr->GetHitAB(g2)->GetDCEnergy());
      Fill(Form("egamegamABdc_sym_%s",oname.Data()),
	   3500,0,3500,gr->GetHitAB(g1)->GetDCEnergy(),
	   3500,0,3500,gr->GetHitAB(g2)->GetDCEnergy());
      Fill(Form("egamegamABdc_sym_%s",oname.Data()),
	   3500,0,3500,gr->GetHitAB(g2)->GetDCEnergy(),
	   3500,0,3500,gr->GetHitAB(g1)->GetDCEnergy());
    }
  }


  // Various gamma-gamma spectra with an additional time gate.
  {
    int highestg=-1;
    double highesten =0;
    for(UShort_t g=0;g<gr->GetMult();g++){
      HitCalc* hit = gr->GetHit(g);
      if(hit->GetEnergy()<1 || hit->GetEnergy()>6000)
	continue;
      if(foundTCut && TimeCut->IsInside(hit->GetTime(),hit->GetEnergy())){
	if(hit->GetDCEnergy()>highesten){
	  highesten = hit->GetDCEnergy();
	  highestg =g;
	}
      }
    }
    for(UShort_t g=0;g<gr->GetMult();g++){
      HitCalc* hit = gr->GetHit(g);
      if(hit->GetEnergy()<1 || hit->GetEnergy()>6000)
	continue;
      if(foundTCut && TimeCut->IsInside(hit->GetTime(),hit->GetEnergy())){
	if(highestg>-1 && g!=highestg){
	  Fill(Form("egamegamdc_fold_%s_tcut",oname.Data()),
	       8000,0,8000,highesten,
	       8000,0,8000,hit->GetDCEnergy());
	}
      }
    }
    highestg=-1;
    highesten =0;
    for(UShort_t g=0;g<gr->GetMultAB();g++){
      HitCalc* hit = gr->GetHitAB(g);
      if(hit->GetEnergy()<1 || hit->GetEnergy()>6000)
	continue;
      if(foundTCut && TimeCut->IsInside(hit->GetTime(),hit->GetEnergy())){
	if(hit->GetDCEnergy()>highesten){
	  highesten = hit->GetDCEnergy();
	  highestg =g;
	}
      }
    }
    for(UShort_t g=0;g<gr->GetMultAB();g++){
      HitCalc* hit = gr->GetHitAB(g);
      if(hit->GetEnergy()<1 || hit->GetEnergy()>6000)
	continue;
      if(foundTCut && TimeCut->IsInside(hit->GetTime(),hit->GetEnergy())){
	if(highestg>-1 && g!=highestg){
	  Fill(Form("egamegamABdc_fold_%s_tcut",oname.Data()),
	       8000,0,8000,highesten,
	       8000,0,8000,hit->GetDCEnergy());
	}
      }
    }
  }
  for(UShort_t g1=0;g1<gr->GetMult();g1++){
    for(UShort_t g2=g1+1;g2<gr->GetMult();g2++){
      if(foundTCut && TimeCut->IsInside(gr->GetHit(g1)->GetTime(),gr->GetHit(g1)->GetEnergy()) &&
	 TimeCut->IsInside(gr->GetHit(g2)->GetTime(),gr->GetHit(g2)->GetEnergy())){
	Fill(Form("egamegamdc_%s_tcut",oname.Data()),
	     3500,0,3500,gr->GetHit(g1)->GetDCEnergy(),
	     3500,0,3500,gr->GetHit(g2)->GetDCEnergy());
	Fill(Form("egamegamdc_sym_%s_tcut",oname.Data()),
	     3500,0,3500,gr->GetHit(g1)->GetDCEnergy(),
	     3500,0,3500,gr->GetHit(g2)->GetDCEnergy());
	Fill(Form("egamegamdc_sym_%s_tcut",oname.Data()),
	     3500,0,3500,gr->GetHit(g2)->GetDCEnergy(),
	     3500,0,3500,gr->GetHit(g1)->GetDCEnergy());
      }
    }
  }
  for(UShort_t g1=0;g1<gr->GetMultAB();g1++){
    for(UShort_t g2=g1+1;g2<gr->GetMultAB();g2++){
      if(foundTCut && TimeCut->IsInside(gr->GetHitAB(g1)->GetTime(),gr->GetHitAB(g1)->GetEnergy()) &&
	 TimeCut->IsInside(gr->GetHitAB(g2)->GetTime(),gr->GetHitAB(g2)->GetEnergy())){
	Fill(Form("egamegamABdc_%s_tcut",oname.Data()),
	     3500,0,3500,gr->GetHitAB(g1)->GetDCEnergy(),
	     3500,0,3500,gr->GetHitAB(g2)->GetDCEnergy());
	Fill(Form("egamegamABdc_sym_%s_tcut",oname.Data()),
	     3500,0,3500,gr->GetHitAB(g1)->GetDCEnergy(),
	     3500,0,3500,gr->GetHitAB(g2)->GetDCEnergy());
	Fill(Form("egamegamABdc_sym_%s_tcut",oname.Data()),
	     3500,0,3500,gr->GetHitAB(g2)->GetDCEnergy(),
	     3500,0,3500,gr->GetHitAB(g1)->GetDCEnergy());
      }
    }
  }

  int hp = s800->GetRegistr();
  Fill(Form("registr_%s",oname.Data()),
       16,-0.5,15.5,hp);
  for(int j=0;j<16;j++){
    if(hp & (1<<j)){
      Fill(Form("trigbit_%s",oname.Data()),
	   16,-0.5,15.5,j);
    }
  }

  if(fCal&(1<<2)){
    for(UShort_t p=0; p<2;p++){
      Fill(Form("ICde_vs_x%d_%s",p,oname.Data()),
	   600,-300,300,pad[p]->GetX(),
	   fIC_range[0],fIC_range[1],fIC_range[2],ich->GetDE());
      Fill(Form("ICde_vs_y%d_%s",p,oname.Data()),
	   600,-100,100,pad[p]->GetY(),
	   fIC_range[0],fIC_range[1],fIC_range[2],ich->GetDE());
      Fill(Form("ICsum_vs_x%d_%s",p,oname.Data()),
	   600,-300,300,pad[p]->GetX(),
	   fIC_range[0],fIC_range[1],fIC_range[2],ich->GetSum());
      Fill(Form("ICsum_vs_y%d_%s",p,oname.Data()),
	   600,-100,100,pad[p]->GetY(),
	   fIC_range[0],fIC_range[1],fIC_range[2],ich->GetSum());
    }
    Fill(Form("x_vs_obj_%s",oname.Data()),
	 fobj_range[0],fobj_range[1],fobj_range[2],tof->GetOBJ(),
	 600/2,-300,300,pad[0]->GetX());
    Fill(Form("afp_vs_obj_%s",oname.Data()),
	 fobj_range[0],fobj_range[1],fobj_range[2],tof->GetOBJ(),
	 100/2,-100,100,track->GetAFP());
    Fill(Form("x_vs_objC_%s",oname.Data()),
	 fobjC_range[0],fobjC_range[1],fobjC_range[2],tof->GetOBJC(),
	 600/2,-300,300,pad[0]->GetX());
    Fill(Form("afp_vs_objC_%s",oname.Data()),
	 fobjC_range[0],fobjC_range[1],fobjC_range[2],tof->GetOBJC(),
	 100/2,-100,100,track->GetAFP());
    Fill(Form("x_vs_xfp_%s",oname.Data()),
	 fxfp_range[0],fxfp_range[1],fxfp_range[2],tof->GetXFP(),
	 600/2,-300,300,pad[0]->GetX());
    Fill(Form("afp_vs_xfp_%s",oname.Data()),
	 fxfp_range[0],fxfp_range[1],fxfp_range[2],tof->GetXFP(),
	 100/2,-100,100,track->GetAFP());
    Fill(Form("x_vs_xfpC_%s",oname.Data()),
	 fxfpC_range[0],fxfpC_range[1],fxfpC_range[2],tof->GetXFPC(),
	 600/2,-300,300,pad[0]->GetX());
    Fill(Form("afp_vs_xfpC_%s",oname.Data()),
	 fxfpC_range[0],fxfpC_range[1],fxfpC_range[2],tof->GetXFPC(),
	 100/2,-100,100,track->GetAFP());

    Fill(Form("x_vs_objtac_%s",oname.Data()),
	 ftacobj_range[0],ftacobj_range[1],ftacobj_range[2],tof->GetTACOBJ(),
	 600/2,-300,300,pad[0]->GetX());
    Fill(Form("afp_vs_objtac_%s",oname.Data()),
	 ftacobj_range[0],ftacobj_range[1],ftacobj_range[2],tof->GetTACOBJ(),
	 100/2,-100,100,track->GetAFP());
    Fill(Form("x_vs_objtacC_%s",oname.Data()),
	 ftacobjC_range[0],ftacobjC_range[1],ftacobjC_range[2],tof->GetTACOBJC(),
	 600/2,-300,300,pad[0]->GetX());
    Fill(Form("afp_vs_objtacC_%s",oname.Data()),
	 ftacobjC_range[0],ftacobjC_range[1],ftacobjC_range[2],tof->GetTACOBJC(),
	 100/2,-100,100,track->GetAFP());
    Fill(Form("x_vs_xfptac_%s",oname.Data()),
	 ftacxfp_range[0],ftacxfp_range[1],ftacxfp_range[2],tof->GetTACXFP(),
	 600/2,-300,300,pad[0]->GetX());
    Fill(Form("afp_vs_xfptac_%s",oname.Data()),
	 ftacxfp_range[0],ftacxfp_range[1],ftacxfp_range[2],tof->GetTACXFP(),
	 100/2,-100,100,track->GetAFP());
    Fill(Form("x_vs_xfptacC_%s",oname.Data()),
	 ftacxfpC_range[0],ftacxfpC_range[1],ftacxfpC_range[2],tof->GetTACXFPC(),
	 600/2,-300,300,pad[0]->GetX());
    Fill(Form("afp_vs_xfptacC_%s",oname.Data()),
	 ftacxfpC_range[0],ftacxfpC_range[1],ftacxfpC_range[2],tof->GetTACXFPC(),
	 100/2,-100,100,track->GetAFP());
  }

  Fill(Form("xfp_%s",oname.Data()),
       600,-300,300,track->GetXFP());
  Fill(Form("xfpf_%s",oname.Data()),
       600,-0.3,0.3,track->GetXFP()/1000.);
  Fill(Form("xfpafp_%s",oname.Data()),
       120,-60,60,track->GetAFP(),
       600,-300,300,track->GetXFP());
  Fill(Form("xfpazita_%s",oname.Data()),
       360,0,360,track->GetPhi(),
       600,-300,300,track->GetXFP());
  Fill(Form("afp_%s",oname.Data()),
       400,-200,200,track->GetAFP());
  Fill(Form("yfp_%s",oname.Data()),
       200,-100,100,track->GetYFP());
  Fill(Form("bfp_%s",oname.Data()),
       100,-50,50,track->GetBFP());
  Fill(Form("yta_%s",oname.Data()),
       100,-50,50,track->GetYTA());
  Fill(Form("dta_%s",oname.Data()),
       200,-10,10,track->GetDTA());
  Fill(Form("azita_%s",oname.Data()),
       360,0,360,track->GetPhi());
  Fill(Form("ptot_%s",oname.Data()),
       500,8,20,track->GetPtot());
  Fill(Form("ppar_%s",oname.Data()),
       500,8,20,track->GetPpar());
  Fill(Form("ptra_%s",oname.Data()),
       200,0,2,track->GetPtra());
  Fill(Form("etot_%s",oname.Data()),
       1000,2000,3000,track->GetEtot());
  Fill(Form("dta_vs_ata_%s",oname.Data()),
       120,-60,60,track->GetATA(),
       400,-10,10,track->GetDTA());


#ifdef S800_DETAILEDTREE
  if(fCal&(1<<1)){
    for(UShort_t c=0;c<ich->GetChan().size();c++){
      Fill(Form("ic_ch%d_%s",ich->GetChan()[c],oname.Data()),
	   2000,0,2000,ich->GetCal()[c]);
      Fill(Form("ic_%s",oname.Data()),
	   16,0,16,ich->GetChan()[c],
	   2500,0,2500,ich->GetCal()[c]);
    }
  }
  if(fCal&(1<<0)) {
    for(UShort_t p=0; p<2;p++){
      if(pad[p]->GetMaxPad()>0&&pad[p]->GetMaxPad()<224){
	Fill(Form("crdcmaxpad_%d_%s",p,oname.Data()),
	     224,0,224,pad[p]->GetMaxPad(),
	     1024,0,1024,pad[p]->GetPadMax()[1]);

	Fill(Form("crdcpad_%d_%d_%s",p,pad[p]->GetMaxPad(),oname.Data()),
	     1024,0,1024,pad[p]->GetPadMax()[1]);

	//max and neighbours

	Fill(Form("crdcmaxpad0_%d_%s",p,oname.Data()),
	     224,0,224,pad[p]->GetMaxPad(),
	     1024,0,1024,pad[p]->GetPadMax()[1]);
	Fill(Form("crdcmaxpad0_%d_%s",p,oname.Data()),
	     224,0,224,pad[p]->GetMaxPad()-1,
	     1024,0,1024,pad[p]->GetPadMax()[0]);
	Fill(Form("crdcmaxpad0_%d_%s",p,oname.Data()),
	     224,0,224,pad[p]->GetMaxPad()+1,
	     1024,0,1024,pad[p]->GetPadMax()[2]);

	Fill(Form("crdcpad0_%d_%d_%s",p,pad[p]->GetMaxPad(),oname.Data()),
	     1024,0,1024,pad[p]->GetPadMax()[1]);
	Fill(Form("crdcpad0_%d_%d_%s",p,pad[p]->GetMaxPad()-1,oname.Data()),
	     1024,0,1024,pad[p]->GetPadMax()[0]);
	Fill(Form("crdcpad0_%d_%d_%s",p,pad[p]->GetMaxPad()+1,oname.Data()),
	     1024,0,1024,pad[p]->GetPadMax()[2]);


	Fill(Form("xpad_%d_%s",p,oname.Data()),
	     200,-2,2,pad[p]->GetXGravity()-pad[p]->GetXFit(),
	     224,0,224,pad[p]->GetMaxPad());

	Fill(Form("xpad0_%d_%s",p,oname.Data()),
	     224,0,224,pad[p]->GetMaxPad(),
	     1024,0,1024,pad[p]->GetPadMax()[0]);
	Fill(Form("xpad0_%d_%s",p,oname.Data()),
	     224,0,224,pad[p]->GetMaxPad(),
	     1024,0,1024,pad[p]->GetPadMax()[2]);

	Fill(Form("xpadl_%d_%s",p,oname.Data()),
	     224,0,224,pad[p]->GetMaxPad(),
	     1000,-500,500,pad[p]->GetPadMax()[1]-pad[p]->GetPadMax()[0]);
	Fill(Form("xpadr_%d_%s",p,oname.Data()),
	     224,0,224,pad[p]->GetMaxPad(),
	     1000,-500,500,pad[p]->GetPadMax()[1]-pad[p]->GetPadMax()[2]);

	Fill(Form("xgra_%d_%s",p,oname.Data()),
	     200,-2,2,pad[p]->GetXGravity()-pad[p]->GetMaxPad(),
	     224,0,224,pad[p]->GetMaxPad());

	Fill(Form("xfit_%d_%s",p,oname.Data()),
	     200,-2,2,pad[p]->GetXFit()-pad[p]->GetMaxPad(),
	     224,0,224,pad[p]->GetMaxPad());
      }//pads
    }//crdc

  }
#endif

}
void CalHistograms::FillHistogramsGateOnlyOut(GretinaCalc* gr, S800Calc* s800, Mode3Calc* m3c, const char* outname){
  //in case you have pure beam, fill here
}
