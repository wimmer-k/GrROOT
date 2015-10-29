// Fitting routines for corrections to the calibration for S800
// K. Wimmer (NSCL)
// root -l corrections.C
// or on the root commandline
// .x corrections.C
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <cmath>

#include "TGraph.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TEnv.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraph.h"
#include "TFile.h"
#include "TF1.h"
#include "TMath.h"
#define isnan(x) ((x) != (x))
using namespace std;
Float_t OBJe1 = -1;
Float_t OBJ = -1;
Float_t XFPe1 = -1;
Float_t XFP = -1;
Float_t TACOBJe1 = -1;
Float_t TACOBJ = -1;
Float_t TACXFPe1 = -1;
Float_t TACXFP = -1;

Int_t steps = -1;
Float_t start = -1;
Float_t inc = -1;

Int_t fiter =0;

Float_t ICXCorr = -1;
Float_t ICX0Corr = -1;

void timecor(char* filename, int afp, int iter = 0);
void afpcor(char* filename, float ustart, float ustop, int usteps, int iter = 0);
void momcor(int type, float ustart, float ustop, int usteps);
void iicor(int type);
void iicor_pro(int type);
void xcor(char* filename, float ustart, float ustop, int usteps, int iter = 0);
void reset();
void saveset(char* filename);
Double_t fscheitel(Double_t *x, Double_t *par);
void corrections(){
  cout << "================================================================================" << endl;
  cout << "Corrections for "<<endl;
  cout << "- time and IC depending on the position and angle in the focal planeof the S800" << endl;
  cout << "- the momentum dependence on the incoming angle at the intermediate image" << endl;
  cout << "     K. Wimmer (NSCL) " << endl;
  cout << "--------------------------------------------------------------------------------" << endl;  
  cout << "- for time corrections run:" << endl;
  cout << "  \"afpcor(\"FILENAME\", float \"FROM\", float \"TO\", int \"STEPS\", int \"ITER\")\"" << endl;
  cout << "  (you should run this before the x correction and\n  apply the new settings after it.)" << endl;
  cout << "  \"xcor(\"FILENAME\", float \"FROM\", float \"TO\", int \"STEPS\", int \"ITER\")\"" << endl;
  cout << "- for IC corrections run:" << endl;
  cout << "  \"iccor(\"FILENAME\")\"" << endl;
    cout << "- for II tracking corrections run:" << endl;
  cout << "  \"momcor(int \"TYPE\", float \"FROM\", float \"TO\", int \"STEPS\")\"" << endl;
  cout << "--------------------------------------------------------------------------------" << endl;  
  cout << "- to save your new corrections type:" << endl;
  cout << "  \"saveset(\"FILENAME\")\"" << endl;
  cout << "- to reset:" << endl;
  cout << "  \"reset()\"" << endl;  
  cout << "================================================================================" << endl;
  reset();
}
void reset(){
  OBJe1 = sqrt(-1.0);
  OBJ = sqrt(-1.0);
  XFPe1 = sqrt(-1.0);
  XFP = sqrt(-1.0);
  TACOBJe1 = sqrt(-1.0);
  TACOBJ = sqrt(-1.0);
  TACXFPe1 = sqrt(-1.0);
  TACXFP = sqrt(-1.0);

  ICXCorr = sqrt(-1.0);
  ICX0Corr = sqrt(-1.0);


}
void saveset(char* filename){
  TEnv *set = new TEnv(filename);
  if(fiter == 0){
    if(!isnan(OBJe1))
      set->SetValue("OBJe1.Corr",OBJe1);
    if(!isnan(OBJ))
      set->SetValue("OBJ.Corr",OBJ);
    if(!isnan(TACOBJe1))
      set->SetValue("TACOBJe1.Corr",TACOBJe1);
    if(!isnan(TACOBJ))
      set->SetValue("TACOBJ.Corr",TACOBJ);
    if(!isnan(XFPe1))
      set->SetValue("XFPe1.Corr",XFPe1);
    if(!isnan(XFP))
      set->SetValue("XFP.Corr",XFP);
    if(!isnan(TACXFPe1))
      set->SetValue("TACXFPe1.Corr",TACXFPe1);
    if(!isnan(TACXFP))
      set->SetValue("TACXFP.Corr",TACXFP);
    if(!isnan(ICXCorr))
      set->SetValue("IC.X.Corr",ICXCorr);
    if(!isnan(ICX0Corr))
      set->SetValue("IC.X0.Corr",ICX0Corr);
  }
  else{
    if(!isnan(OBJe1))
      set->SetValue("OBJe1.Corr",set->GetValue("OBJe1.Corr",0.0)+OBJe1);
    if(!isnan(OBJ))
      set->SetValue("OBJ.Corr",set->GetValue("OBJ.Corr",0.0)+OBJ);
    if(!isnan(TACOBJe1))
      set->SetValue("TACOBJe1.Corr",set->GetValue("TACOBJe1.Corr",TACOBJe1)+TACOBJe1);
    if(!isnan(TACOBJ))
      set->SetValue("TACOBJ.Corr",set->GetValue("TACOBJ.Corr",TACOBJ)+TACOBJ);
    if(!isnan(XFPe1))
      set->SetValue("XFPe1.Corr",set->GetValue("XFPe1.Corr",0.0)+XFPe1);
    if(!isnan(XFP))
      set->SetValue("XFP.Corr",set->GetValue("XFP.Corr",0.0)+XFP);
    if(!isnan(TACXFPe1))
      set->SetValue("TACXFPe1.Corr",set->GetValue("TACXFPe1.Corr",TACXFPe1)+TACXFPe1);
    if(!isnan(TACXFP))
      set->SetValue("TACXFP.Corr",set->GetValue("TACXFP.Corr",TACXFP)+TACXFP);
    if(!isnan(ICXCorr))
      set->SetValue("IC.X.Corr",ICXCorr);
    if(!isnan(ICX0Corr))
      set->SetValue("IC.X0.Corr",ICX0Corr);
  }
  set->SaveLevel(kEnvLocal);
}
void momcor(int type, float ustart, float ustop, int usteps){
  steps = usteps;
  start = ustart;
  inc = (ustop-ustart)/usteps;
  iicor(type);
}
void afpcor(char* filename, float ustart, float ustop, int usteps, int iter){
  steps = usteps;
  start = ustart;
  inc = (ustop-ustart)/usteps;
  timecor(filename, 1, iter);
  fiter = iter;
}
void xcor(char* filename, float ustart, float ustop, int usteps, int iter){
  steps = usteps;
  start = ustart;
  inc = (ustop-ustart)/usteps;
  timecor(filename, 0, iter);
  fiter = iter;
}
void timecor(char* filename, int afp, int iter){
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1); 
  TFile *f = new TFile(filename);

  TH2F *h[4];
  TCanvas *c1 = new TCanvas("c","c",0,0,600,400);
  c1->Divide(4,2);
  c1->cd(1);
  
  //adjust here!
  double fitlow[4];
  double fithigh[4];
  fitlow[0] = 1900;
  fitlow[1] = 50;
  fitlow[2] = 1750;
  fitlow[3] = 1950;
  fithigh[0] = 2000;
  fithigh[1] = 250;
  fithigh[2] = 1850;
  fithigh[3] = 2050;
  
  if(afp>0){
    if(iter == 0){
      h[0] = (TH2F*)f->Get("afp_vs_obj");
      h[1] = (TH2F*)f->Get("afp_vs_xfp");
      h[2] = (TH2F*)f->Get("afp_vs_objtac");
      h[3] = (TH2F*)f->Get("afp_vs_xfptac");
    }
    else{
      h[0] = (TH2F*)f->Get("afp_vs_objC");
      h[1] = (TH2F*)f->Get("afp_vs_xfpC");
      h[2] = (TH2F*)f->Get("afp_vs_objtacC");
      h[3] = (TH2F*)f->Get("afp_vs_xfptacC");
    }
  }
  else{
    if(iter == 0){
      h[0] = (TH2F*)f->Get("x_vs_obj");
      h[1] = (TH2F*)f->Get("x_vs_xfp");
      h[2] = (TH2F*)f->Get("x_vs_objtac");
      h[3] = (TH2F*)f->Get("x_vs_xfptac");
    }
    else{
      h[0] = (TH2F*)f->Get("x_vs_objC");
      h[1] = (TH2F*)f->Get("x_vs_xfpC");
      h[2] = (TH2F*)f->Get("x_vs_objtacC");
      h[3] = (TH2F*)f->Get("x_vs_xfptacC");
    }
  }
  int from =0;
  int to = 430;
  TH2F *hcor;
  TH2F* hres[4];
  TH1D *hpro[4][200];
  double cor[4][200];
  double width[4][200];
  for(int j=0;j<4;j++){
    int nx = h[j]->GetNbinsX();
    int lx = h[j]->GetXaxis()->GetBinLowEdge(1);
    int ux = h[j]->GetXaxis()->GetBinUpEdge(h[j]->GetNbinsX());
    int ny = h[j]->GetNbinsY();
    int ly = h[j]->GetYaxis()->GetBinLowEdge(1);
    int uy = h[j]->GetYaxis()->GetBinUpEdge(h[j]->GetNbinsY());
    for(int i=0;i<steps;i++){
      if(hcor)
	delete hcor;
      hcor = new TH2F(Form("corrected%d%d",j,i),Form("corrected%d%d",j,i),nx,lx,ux,ny,ly,uy);
      TF1 *wi = new TF1("w","gaus",fitlow[j],fithigh[j]);
      wi->SetLineColor(3);
      wi->SetLineWidth(1);
     cor[j][i] = start+i*inc;
      for(int x=1;x<nx;x++){//obj
	double xp = h[j]->GetXaxis()->GetBinCenter(x);
	for(int y=1;y<ny;y++){//afp
	  double yp = h[j]->GetYaxis()->GetBinCenter(y);
	  if(h[j]->GetBinContent(x,y)>0)
	    hcor->Fill(xp+cor[j][i]*yp,yp,h[j]->GetBinContent(x,y));
	}
	
      }  
      if(afp>0)
	hpro[j][i] = hcor->ProjectionX(Form("hpro%d%d",j,i));
      else
	hpro[j][i] = hcor->ProjectionX(Form("hpro%d%d",j,i),from,to);
      hpro[j][i]->Fit(wi,"Rq");
      width[j][i] = wi->GetParameter(2);
    }
  }
  double result[4];
  for(int j=0;j<4;j++){
    int nx = h[j]->GetNbinsX();
    int lx = h[j]->GetXaxis()->GetBinLowEdge(1);
    int ux = h[j]->GetXaxis()->GetBinUpEdge(h[j]->GetNbinsX());
    int ny = h[j]->GetNbinsY();
    int ly = h[j]->GetYaxis()->GetBinLowEdge(1);
    int uy = h[j]->GetYaxis()->GetBinUpEdge(h[j]->GetNbinsY());
   c1->cd(1+j);
    TGraph *g = new TGraph(steps,&cor[j][0],&width[j][0]);
    g->SetMarkerColor(3);
    g->SetMarkerStyle(20);
    g->Draw("AP");
    TF1 *res = new TF1("res","pol2",0.5,1.5);
    g->Fit(res,"q");
    
    result[j] = -res->GetParameter(1)/2./res->GetParameter(2);
    cout <<h[j]->GetName() <<"\t"<< result[j] << endl;
    hres[j] = new TH2F(Form("result%d",j),Form("result%d",j),nx,lx,ux,ny,ly,uy);
    for(int x=1;x<nx;x++){//obj
      double xp = h[j]->GetXaxis()->GetBinCenter(x);
      for(int y=1;y<ny;y++){//afp of x
	double yp = h[j]->GetYaxis()->GetBinCenter(y);
	if(h[j]->GetBinContent(x,y)>0)
	  hres[j]->Fill(xp+result[j]*yp,yp,h[j]->GetBinContent(x,y));
      }
      
    }
    c1->cd(5+j);
    hres[j]->Draw("colz");
  }
  if(afp>0){
    OBJe1 = result[0];
    XFPe1 = result[1];
    TACOBJe1 = result[2];
    TACXFPe1 = result[3];
  }
  OBJ = result[0];
  XFP = result[1];
  TACOBJ = result[2];
  TACXFP = result[3];
 
}
void iicor_pro(int type){
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1); 
  TFile *f[5];
  f[0] = new TFile("../unntruns/hrun001.root");
  f[1] = new TFile("../unntruns/hrun002.root");
  f[2] = new TFile("../unntruns/hrun084.root");
  f[3] = new TFile("../un9runs/hrun003.root");
  f[4] = new TFile("../un100runs/hrun039.root");

  TH2F *h[5];
  TCanvas *c1 = new TCanvas("c","c",0,0,600,400);
  c1->Divide(5,3);
  
  for(int i=0;i<5;i++){
    if(type==0)
      h[i] = (TH2F*)f[i]->Get("dta_vs_aii");
    if(type==1)
      h[i] = (TH2F*)f[i]->Get("ptot_vs_aii_r");
    if(type==2)
      h[i] = (TH2F*)f[i]->Get("ppar_vs_aii_r");
    if(type==3)
      h[i] = (TH2F*)f[i]->Get("etot_vs_aii_r");
  }

  TH2F *hcor;
  TH2F* hres[5];
  TProfile *hpro[5];
  TF1 *wi[5];
  wi[0] = new TF1("wi0","pol1",-50,50);
  wi[1] = new TF1("wi1","pol1",-75,100);
  wi[2] = new TF1("wi2","pol1",-50,50);
  wi[3] = new TF1("wi3","pol1",-75,100);
  wi[4] = new TF1("wi4","pol1",-50,50);
  double min[5];
  double mini[5];
  for(int j=0;j<5;j++){
    c1->cd(1+j);
    h[j]->Draw("colz");
    int nx = h[j]->GetNbinsX();
    float lx = h[j]->GetXaxis()->GetBinLowEdge(1);
    float ux = h[j]->GetXaxis()->GetBinUpEdge(h[j]->GetNbinsX());
    int ny = h[j]->GetNbinsY();
    float ly = h[j]->GetYaxis()->GetBinLowEdge(1);
    float uy = h[j]->GetYaxis()->GetBinUpEdge(h[j]->GetNbinsY());
    hpro[j] = (TProfile*)h[j]->ProfileX(Form("hpro%d",j));
    c1->cd(6+j);
    hpro[j]->Draw();
    wi[j]->SetLineColor(3);
    wi[j]->SetLineWidth(1);
    hpro[j]->Fit(wi[j],"R");
    hres[j] = new TH2F(Form("result%d",j),Form("result%d",j),nx,lx,ux,ny,ly,uy);
    for(int x=1;x<nx;x++){//obj
      double xp = h[j]->GetXaxis()->GetBinCenter(x);
      for(int y=1;y<ny;y++){//afp of x
	double yp = h[j]->GetYaxis()->GetBinCenter(y);
	if(h[j]->GetBinContent(x,y)>0)
	  hres[j]->Fill(xp,yp-wi[j]->GetParameter(1)*xp,h[j]->GetBinContent(x,y));
      }
      
    }
    c1->cd(11+j);
    hres[j]->Draw("colz");
  }
}
void iicor(int type){
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1); 
  TFile *f[5];

  //adjust here!
  f[0] = new TFile("../unntruns/hrun001.root");
  f[1] = new TFile("../unntruns/hrun002.root");
  f[2] = new TFile("../unntruns/hrun084.root");
  f[3] = new TFile("../un9runs/hrun003.root");
  f[4] = new TFile("../un100runs/hrun039.root");

  TH2F *h[5];
  TCanvas *c1 = new TCanvas("c","c",0,0,600,400);
  c1->Divide(5,2);
  c1->cd(1);
  
  for(int i=0;i<5;i++){
    if(type==0)
      h[i] = (TH2F*)f[i]->Get("dta_vs_aii");
    if(type==1)
      h[i] = (TH2F*)f[i]->Get("ptot_vs_aii_r");
    if(type==2)
      h[i] = (TH2F*)f[i]->Get("ppar_vs_aii_r");
    if(type==3)
      h[i] = (TH2F*)f[i]->Get("etot_vs_aii_r");
  }

  TH2F *hcor;
  TH2F* hres[5];
  TH1D *hpro[5][200];
  TF1 *wi[5][200];
  double cor[5][200];
  double width[5][200];
  double min[5];
  double mini[5];
  for(int j=0;j<5;j++){
    min[j] = 1e6;
    mini[j] = -1;
    int nx = h[j]->GetNbinsX();
    float lx = h[j]->GetXaxis()->GetBinLowEdge(1);
    float ux = h[j]->GetXaxis()->GetBinUpEdge(h[j]->GetNbinsX());
    int ny = h[j]->GetNbinsY();
    float ly = h[j]->GetYaxis()->GetBinLowEdge(1);
    float uy = h[j]->GetYaxis()->GetBinUpEdge(h[j]->GetNbinsY());
    for(int i=0;i<steps;i++){
      if(hcor)
	delete hcor;
      hcor = new TH2F(Form("corrected%d%d",j,i),Form("corrected%d%d",j,i),nx,lx,ux,ny,ly,uy);
      cor[j][i] = start+i*inc;
      for(int x=1;x<nx;x++){
	double xp = h[j]->GetXaxis()->GetBinCenter(x);
	for(int y=1;y<ny;y++){
	  double yp = h[j]->GetYaxis()->GetBinCenter(y);
	  if(h[j]->GetBinContent(x,y)>0)
	    hcor->Fill(xp,yp+cor[j][i]*xp,h[j]->GetBinContent(x,y));
	}
	
      }  
      hpro[j][i] = hcor->ProjectionY(Form("hpro%d%d",j,i));
      //cout << j <<"\t" << i<<"\t" << hpro[j][i]->GetMaximumBin()<<endl;
      //cout << j <<"\t" << i<<"\t" << (hpro[j][i]->GetMaximumBin()-3)*hpro[j][i]->GetBinWidth(0)+hpro[j][i]->GetBinCenter(0)<<"\t" << (hpro[j][i]->GetMaximumBin()+3)*hpro[j][i]->GetBinWidth(0)+hpro[j][i]->GetBinCenter(0)<< endl;
      wi[j][i] = new TF1(Form("w%d%d",j,i),"gaus",(hpro[j][i]->GetMaximumBin()-15)*hpro[j][i]->GetBinWidth(0)+hpro[j][i]->GetBinCenter(0),(hpro[j][i]->GetMaximumBin()+15)*hpro[j][i]->GetBinWidth(0)+hpro[j][i]->GetBinCenter(0));
      wi[j][i]->SetLineColor(3);
      wi[j][i]->SetLineWidth(1);
      hpro[j][i]->Fit(wi[j][i],"qR");
      width[j][i] = wi[j][i]->GetParameter(2);
      if(width[j][i]<min[j]){
	min[j] = width[j][i];
	mini[j] = cor[j][i];
      }
    }
  }
  double result[5];
  for(int j=0;j<5;j++){
    int nx = h[j]->GetNbinsX();
    float lx = h[j]->GetXaxis()->GetBinLowEdge(1);
    float ux = h[j]->GetXaxis()->GetBinUpEdge(h[j]->GetNbinsX());
    int ny = h[j]->GetNbinsY();
    float ly = h[j]->GetYaxis()->GetBinLowEdge(1);
    float uy = h[j]->GetYaxis()->GetBinUpEdge(h[j]->GetNbinsY());
    c1->cd(1+j);
    TGraph *g = new TGraph(steps,&cor[j][0],&width[j][0]);
    g->SetMarkerColor(3);
    g->SetMarkerStyle(20);
    g->Draw("AP");
    if(steps>50)
      TF1 *res = new TF1("res",fscheitel,mini[j]-20*inc,mini[j]+20*inc,3);
    else
      TF1 *res = new TF1("res",fscheitel,start,start+steps*inc,3);
    res->SetParameter(1,mini[j]);
    res->SetParameter(2,min[j]);
    //cout << min[j] << "\t" << mini[j] << endl;
    g->Fit(res,"qR");
    result[j] = res->GetParameter(1);
    cout <<h[j]->GetName() <<"\t"<< result[j] << endl;
    hres[j] = new TH2F(Form("result%d",j),Form("result%d",j),nx,lx,ux,ny,ly,uy);
    for(int x=1;x<nx;x++){
      double xp = h[j]->GetXaxis()->GetBinCenter(x);
      for(int y=1;y<ny;y++){
	double yp = h[j]->GetYaxis()->GetBinCenter(y);
	if(h[j]->GetBinContent(x,y)>0)
	  hres[j]->Fill(xp,yp+result[j]*xp,h[j]->GetBinContent(x,y));
      }
      
    }
    c1->cd(6+j);
    hres[j]->Draw("colz");
  }
 
}
Double_t fexpo(Double_t *x, Double_t *par){
  if(x[0]<par[2])
    return par[0]*exp(-par[1]*(par[2]-x[0]));
  return par[0];
}
void iccor(TH2F* h){
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);

  TCanvas *c2 = new TCanvas("c2","c2",0,0,600,400);
  c2->Divide(2,2);
  c2->cd(1);
  h->Draw("colz");
  TProfile *p = h->ProfileX();
  p->GetYaxis()->SetRangeUser(200,500);
  TF1 *fu = new TF1("fu",fexpo,-300,300,3);
  TH2F *hres;

  //adjust here!
  fu->SetParameter(2,200); //x0
  fu->SetParameter(1,0.0002); //p
  fu->SetParameter(0,700); //E
  //see equation in section 6.4.2 of the manual
  
  fu->SetLineColor(2);
  fu->SetLineWidth(1);
  c2->cd(2);
  p->Draw();
  p->Fit(fu,"R");
  
  ICXCorr = fu->GetParameter(1);
  ICX0Corr = fu->GetParameter(2);
  
  int nx = h->GetNbinsX();
  int lx = h->GetXaxis()->GetBinLowEdge(1);
  int ux = h->GetXaxis()->GetBinUpEdge(h->GetNbinsX());
  int ny = h->GetNbinsY();
  int ly = h->GetYaxis()->GetBinLowEdge(1);
  int uy = h->GetYaxis()->GetBinUpEdge(h->GetNbinsY());
  
  hres = new TH2F("res","res",nx,lx,ux,ny,ly,uy);
  for(int x=1;x<nx;x++){//obj
    double xp = h->GetXaxis()->GetBinCenter(x);
    for(int y=1;y<ny;y++){//afp of x
      double yp = h->GetYaxis()->GetBinCenter(y);
      if(h->GetBinContent(x,y)>0){
	if(xp<fu->GetParameter(2))
	  hres->Fill(xp,yp*exp(fu->GetParameter(1)*(fu->GetParameter(2)-xp)),h->GetBinContent(x,y));
	else
	  hres->Fill(xp,yp,h->GetBinContent(x,y));
      }
    }
  }
  
  c2->cd(3);
  hres->Draw("colz");
  
  c2->cd(4);

  cout << "IC.X.Corr:\t" << ICXCorr << endl;
  cout << "IC.X0.Corr:\t" <<ICX0Corr << endl;


  
}
Double_t fscheitel(Double_t *x, Double_t *par){
  return par[0]*(x[0]-par[1])*(x[0]-par[1])+par[2];
  //scheitepunkt x = par[1], y = par[2]
}
