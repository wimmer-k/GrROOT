#include <iostream>
#include <iomanip>
#include <string>
#include <algorithm>
#include <vector>  

#include "TCanvas.h"
#include "TGraph.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1F.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TImage.h"
#include <TSystem.h>
#include <TMarker.h>
#include "TMath.h"
TCanvas *c;
vector<TMarker*> fmarkers;
vector<double> fpoints;
Double_t fbinw;
TF1 *ffit;
TH1F* fh;
Double_t gaus_lbg(Double_t *x, Double_t *par);
void ClickFit();
void Clear();
void Fit();
void window(){
  fmarkers.clear();
  fpoints.clear();
  c = new TCanvas("window","window",600,400);
  c->AddExec("fit","ClickFit()");
  //TFile *f = new TFile("gamma_new.root");
  // fh = (TH1F*)f->Get("hgamma");
  // fh->Draw();
  fh = NULL;
}
void ClickFit(){
  TObject *select = gPad->GetSelected();
  if(!select)
    return;

  TIter next (c->GetListOfPrimitives());
  TObject *obj;
  while((obj=next())){
    //cout << "reading: " << obj->GetName() << endl;
    if(obj->InheritsFrom("TH1")){
      //cout << "histo: " << obj->GetName() << endl;
      fh = (TH1F*)obj;
    }
  }

  gPad->GetCanvas()->FeedbackMode(kTRUE);
  //Only act on middle-clicks.
  if(gPad->GetEvent() == 12){
    //cout << "middle button " << endl;
    //Find the location of the click.
    double xp = gPad->GetEventX();
    double yp = gPad->GetEventY();
    //cout << xp <<"\t" << yp << endl;
    xp = gPad->AbsPixeltoX((int)xp);
    yp = gPad->AbsPixeltoY((int)yp);
    cout << "point: " << fpoints.size() << " selected at "<< xp << endl;
    fpoints.push_back(xp);
    fmarkers.push_back(new TMarker(xp,yp,2));
    fmarkers.back()->SetMarkerSize(2);
    fmarkers.back()->SetMarkerColor(2);
    for(UShort_t i=0;i<fmarkers.size()-1;i++){
      fmarkers[i]->SetMarkerColor(1);
      fmarkers[i]->Draw("same");
    }
    fmarkers.back()->Draw("same");
    gPad->Modified();
    gPad->Update();
    gSystem->ProcessEvents();

  }
  if(gPad->GetEvent() ==kKeyPress){
    if(fpoints.size()>3 ||fpoints.size()<2)
      return;
	
    switch(gPad->GetEventX()){
    case 'f': //fit
      Fit();
      break;
    default:
      break;
    // case 'z': // ZOOM
    //   mygenerate(1./Z, last_x, last_y);
    //   break;
    // case 'u': // UNZOOM
    //   mygenerate(Z   , last_x, last_y);
    //   break;
    // case 'r': // RESET
    //   mygenerate(-1   , last_x, last_y);
    //   break;
    };
  }
}
void Fit(){
  sort(fpoints.begin(), fpoints.end());
  fbinw = fh->GetBinWidth(1);
  double slope;
  double bg;
  double sigma;
  double content;
  double mean;
  double to;
  if(fpoints.size()==3){
    //cout << fpoints[0] << "\t"<<fpoints[1] << "\t"<<fpoints[2] <<endl;

    slope = fh->GetBinContent(fh->FindBin(fpoints[0])) - fh->GetBinContent(fh->FindBin(fpoints[2]));
    slope/=(fpoints[0]-fpoints[2]);
    bg = fh->GetBinContent(fh->FindBin(fpoints[0])) + fh->GetBinContent(fh->FindBin(fpoints[2]));
    bg-=slope*(fpoints[0]+fpoints[2]);
    bg/=2;

    sigma = 1;
    content = fh->Integral(fh->FindBin(fpoints[0]),fh->FindBin(fpoints[2]));
    mean = fpoints[1];
    to = fpoints[2];
  }
  else if(fpoints.size()==2){
    cout << fpoints[0] << "\t"<<fpoints[1] <<endl;
    slope = fh->GetBinContent(fh->FindBin(fpoints[0])) - fh->GetBinContent(fh->FindBin(fpoints[1]));
    slope/=(fpoints[0]-fpoints[1]);
    bg = fh->GetBinContent(fh->FindBin(fpoints[0])) + fh->GetBinContent(fh->FindBin(fpoints[1]));
    bg-=slope*(fpoints[0]+fpoints[1]);
    bg/=2;

    sigma = 1;
    content = fh->Integral(fh->FindBin(fpoints[0]),fh->FindBin(fpoints[1]));
    mean = fpoints[0]+fpoints[1];
    mean/=2;
    to = fpoints[1];
    cout << slope  << "\t";
    cout << bg  << "\t";
    cout << sigma  << "\t";
    cout << mean  << "\t";
    cout << content  << "\t";
    cout << to  << endl;
  }
  else if(fpoints.size()==1){
    cout << "too few points for fitting! " << endl;
    return;
  }
  else{
    cout << "too many points for fitting! Clearing!" << endl;
    Clear();
    return;
  }
  ffit = new TF1("fit",gaus_lbg,fpoints[0],to,5);
  ffit->SetLineColor(3);
  ffit->SetLineWidth(1);
  ffit->SetParameters(bg,slope,sigma,content,mean);
  ffit->SetParName(0,"BgConstant");
  ffit->SetParName(1,"BgSlope   ");
  ffit->SetParName(2,"Sigma     ");
  ffit->SetParName(3,"Content   ");
  ffit->SetParName(4,"Mean      ");
  
  fh->Fit(ffit,"R");
  //ffit->Draw("same");
  cout << "      Chi Square: " << ffit->GetChisquare() << endl;
  cout << "      FWHM:       " << 2*ffit->GetParameter(2)*sqrt(2*log(2)) << "\t" <<2*ffit->GetParameter(2)*sqrt(2*log(2))/ffit->GetParameter(4) << endl;

  gPad->Modified();
  gPad->Update();
  gSystem->ProcessEvents();
  Clear();

}
void Clear(){
  for(UShort_t i=0;i<fmarkers.size();i++){
    fmarkers[i]->Delete();
  }
  fmarkers.clear();
  fpoints.clear();
  gSystem->ProcessEvents();      

}
Double_t gaus_lbg(Double_t *x, Double_t *par){
/*
  par[0]   background constant
  par[1]   background slope
  par[2]   gauss width
  par[3]   gauss0 constant
  par[4]   gauss0 mean
*/
   static Float_t sqrt2pi = TMath::Sqrt(2*TMath::Pi()), sqrt2 = TMath::Sqrt(2.);
   Double_t arg;
   if (par[2] == 0) 
     par[2]=1; //  force widths /= 0
   arg = (x[0] - par[4])/(sqrt2*par[2]);
   Double_t fitval = par[0] + x[0]*par[1] + fbinw/(sqrt2pi*par[2]) * par[3] * exp(-arg*arg);
   return fitval;
}
