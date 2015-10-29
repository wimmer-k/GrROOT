#include <iostream>
#include <iomanip>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TGraph.h"
#include "TF1.h"
#include "TFile.h"
#include "TH2F.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TImage.h"
#include <TSystem.h>
#include <TGImageMap.h>
#include <TGPicture.h>
#include <TMarker.h>
#include <TEnv.h>

using namespace std;

vector<TMarker*> fMspec;
vector<TMarker*> fMpict;
vector<double> fspecy;
vector<double> fposy;

int frange = 10;
int fselected =-1;
double fxoffset =0;
double fyoffset =0;
double fyslope =1;
int fcrdcNum;

const double fmasky[] = {-25.0,0.0,10.0,20.0,25.0,30.0,50.0};
char* hname = "hcrdcrawrange";
bool cancalibrate = false;

void DynamicHisto();
void DynamicPict();
void Calibrate();
void ChangeName(bool raw){
  if(raw)
    hname = "hcrdcrawrange";
  else
    hname = "hcrdc";
}
void mask(char* fileName, int crdcNum){
  if (crdcNum!=0 && crdcNum!=1){
    cout << "Expecting the crdcs to be number 0 and 1" << endl;
    cout << "quitting" << endl;
    return;
  }
  fcrdcNum = crdcNum;
  
  //gROOT->Reset();
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1); 
  
  TFile *f = new TFile(fileName);
  TH2F* h = (TH2F*)f->Get(Form("%s%dxy",hname,crdcNum));
  TImage *mask = TImage::Open("/GrROOT/scripts/mask.xpm");
  TCanvas *c = new TCanvas("c","c",0,0,600,300);
  c->Divide(2,1);
  

  c->cd(1);
  TPad *ph = new TPad("ph","ph",0,0,1,1);
  ph->Draw();
  ph->cd();
  h->Draw("colz");
  ph->AddExec("histo","DynamicHisto()");
  c->cd(2);
  TPad *pm = new TPad("pm","pm",0,0,1,1);
  pm->Draw();
  pm->cd();
  mask->Draw("xxx");
  pm->AddExec("picture","DynamicPict()");
}

void DynamicPict(){
  TObject *select = gPad->GetSelected();
  if(!select)
    return;

  gPad->GetCanvas()->FeedbackMode(kTRUE);
  //Only act on middle-clicks.
  if(gPad->GetEvent() == 12){

    //Find the location of the click.
    double xp = gPad->GetEventX();
    double yp = gPad->GetEventY();

    //Bail out if you clicked here too often.
    if(fMspec.size()==fMpict.size()-1){
      cout << "Error!! Select a point on the left now" << endl;
      for(UShort_t i=0;i<fMpict.size();i++)
	fMpict[i]->Draw("same");
      gPad->Modified();
      gPad->Update();
      gSystem->ProcessEvents();
      return;
    }
    fselected =-1;

    //Determine which hole on the picture was clicked.
    xp = gPad->AbsPixeltoX((int)xp);
    yp = gPad->AbsPixeltoY((int)yp);
    double min[7] = {0.177, 0.341, 0.447, 0.524, 0.572, 0.610, 0.750};
    double max[7] = {0.245, 0.432, 0.505, 0.572, 0.610, 0.659, 0.8125};
    for(int i=0;i<7;i++){
      if(yp>min[i] && yp<max[i])
	fselected =i;
    }
    if(fselected ==1 && xp>0.518 && xp<0.569)
      fselected = 10;

    cout << "selected " << fselected << endl;
    if(fselected<0){
      for(UShort_t i=0;i<fMpict.size();i++){
	fMpict[i]->Draw("same");
      }
      gPad->Modified();
      gPad->Update();
      gSystem->ProcessEvents();
      return;
    }

    //Store it
    fMpict.push_back(new TMarker(xp,yp,2));
    fMpict.back()->SetMarkerSize(2);
    fMpict.back()->SetMarkerColor(2);
    for(UShort_t i=0;i<fMpict.size()-1;i++){
      fMpict[i]->SetMarkerColor(1);
      fMpict[i]->Draw("same");
    }
    fMpict.back()->Draw("same");

    gPad->Modified();
    gPad->Update();
    gSystem->ProcessEvents();

    if(fselected>-1&&fselected<7)
      fposy.push_back(fmasky[fselected]);
    
    if(fselected ==10)
      fposy.push_back(fmasky[1]);
    
    if(fMspec.size()==fMpict.size() && fselected ==10)
      fxoffset = -fMspec.back()->GetX();

    cancalibrate = false;

    if(fposy.size() == fspecy.size() && fspecy.size()>1)
      cancalibrate = true;

  }
}
void DynamicHisto(){
  TObject *select = gPad->GetSelected();
  if(!select)
    return;
  if(!select->InheritsFrom("TH2")){
    gPad->SetUniqueID(0);
    return;
  }
  TH2F *h = (TH2F*)select->Clone();
  gPad->GetCanvas()->FeedbackMode(kTRUE);

  //Only on middle-clicks
  if(gPad->GetEvent() == 12){

    //Get the click location, convert from pixel location to graph location.
    double xp = gPad->GetEventX();
    double yp = gPad->GetEventY();
    xp = gPad->AbsPixeltoX((int)xp);
    yp = gPad->AbsPixeltoY((int)yp);

    //Zoom in around the position.
    //Use the mean of the zoomed region as the new location.
    h->GetXaxis()->SetRangeUser(xp-frange,xp+frange);
    h->GetYaxis()->SetRangeUser(yp-frange,yp+frange);
    if(h->Integral()<1){
      for(UShort_t i=0;i<fMspec.size();i++)
	fMspec[i]->Draw("same");
      gPad->Modified();
      gPad->Update();
      gSystem->ProcessEvents();
      return;
    }
    xp = h->GetMean(1);
    yp = h->GetMean(2);
      
    cout << xp << "\t "<< yp << endl;

    //You clicked too often on one picture than the other.
    if(fMspec.size()==fMpict.size()+1){
      cout << "Error!! Select a point on the right now" << endl;
      for(UShort_t i=0;i<fMspec.size();i++)
	fMspec[i]->Draw("same");
      gPad->Modified();
      gPad->Update();
      gSystem->ProcessEvents();
      return;
    }
    if(fMspec.size()==fMpict.size() && fselected ==10)
      fxoffset = -fMspec.back()->GetX();
    if(fMspec.size()==fMpict.size()){
      fselected =-1;
    }

    //Make it look pretty by showing you what you clicked on.
    fMspec.push_back(new TMarker(xp,yp,2));
    fMspec.back()->SetMarkerSize(2);
    fMspec.back()->SetMarkerColor(2);
    for(UShort_t i=0;i<fMspec.size()-1;i++){
      fMspec[i]->SetMarkerColor(1);
      fMspec[i]->Draw("same");
    }
    fMspec.back()->Draw("same");
    gPad->Modified();
    gPad->Update();
    gSystem->ProcessEvents();
    fspecy.push_back(yp);

    // //Special case, if the center hole is here, we have our x offset.
    // if(fselected ==10){  // i don't thinkc this is needed anymore
    //   fxoffset = -xp;
    // }
 
    cancalibrate = false;
    if(fposy.size() == fspecy.size() && fspecy.size() >1)
      cancalibrate = true;
    
    fselected =-1;
  }
  delete h;//LR
}
void Calibrate(){
  if(!cancalibrate){
    cout << "Error!! Need more points for calibration!" << endl;
    return;
  }
  if(fspecy.size()!=fposy.size()){
    cout << "Error!! Need same number of points on each side!" << endl;
    cout << "left: " << fspecy.size() << " right: " << fposy.size() << endl;
    return;
  }

  //Do a fit to the calibration points.
  TGraph *g = new TGraph(fspecy.size(),&fspecy[0],&fposy[0]);
  TCanvas* c1 = new TCanvas();
  c1->cd();
  g->Draw("*AP");
  TF1* fit = new TF1("fit","pol1");
  g->Fit(fit);
  fit->Draw("same");
  fyslope = fit->GetParameter(1);
  fyoffset = fit->GetParameter(0);


  cout << "slope " << fyslope << " offset " << fyoffset << endl;

}
void ClearPos(){
  for(UShort_t i=0;i<fMpict.size();i++){
    fMpict[i]->Delete();
  }
  for(UShort_t i=0;i<fMspec.size();i++){
    fMspec[i]->Delete();
  }


  fMpict.clear();
  fMspec.clear();
  fspecy.clear();
  fposy.clear();
  //frowy.clear();
  fselected =-1;
  fxoffset = sqrt(-1.0);
  fyoffset = sqrt(-1.0);
  fyslope = sqrt(-1.0);


  gSystem->ProcessEvents();      

}
void Write(){
  TEnv *env = new TEnv("maskout.dat");
  env->SetValue(Form("Crdc.Y.Offset.%d",fcrdcNum),fyoffset);
  env->SetValue(Form("Crdc.Y.Slope.%d",fcrdcNum),fyslope);
  env->SetValue(Form("Crdc.X.Offset.%d",fcrdcNum),fxoffset);
  env->SetValue(Form("Crdc.X.Slope.%d",fcrdcNum),2.54);
  env->SaveLevel(kEnvLocal);
}
