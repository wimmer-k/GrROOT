#include <string>
#include <string.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <iomanip>

#include "TCanvas.h"
#include "TMath.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1.h"
#include "TF1.h"
#include "TEnv.h"
#include "TSpectrum.h"
#include "TVirtualFitter.h"

using namespace std;

Double_t fgammagaussbg(Double_t *x, Double_t *par);
Double_t fgammabg(Double_t *x, Double_t *par);
Double_t fgammastep(Double_t *x, Double_t *par);
Double_t fgammagaus(Double_t *x, Double_t *par);
void fit(){
  int range = 1e4;

  TFile *f = new TFile("hgrroot.root");
  TCanvas *ca = new TCanvas("ca","ca",0,0,1400,900);
  ca->Divide(4,4);
  TH1F* h[7][4];
  TF1* fu[7][4][2];
  TF1* fus[7][4][2][3];
  double res[7*4];
  double det[7*4];
  double slope[7*4];
  double offset[7*4];
  double resolution =5;
  TEnv* output = new TEnv("corecal.dat");

  for(int d=0;d<7;d++){
    for(int c=0;c<4;c++){
      h[d][c] = (TH1F*)f->Get(Form("htraceen_b%d_c%d_cr%d_d%d",0,9,c,d+1));
      ca->cd(1+c*4);
      //h[d][c]->GetXaxis()->SetRangeUser(100,4000);
      TSpectrum *sp = new TSpectrum(2,resolution);
      sp->SetResolution(resolution);
      Int_t nfound = 0;
      nfound = sp->Search(h[d][c],resolution,"nobackground",0.5);
      Float_t *xpeaks = sp->GetPositionX();
      Float_t *ypeaks = sp->GetPositionY();
//       for(int p=0;p<nfound;p++){
// 	cout << xpeaks[p] << "\t" << ypeaks[p] << endl;
//       }
      if(nfound!=2){
	cout << "Found " << nfound << " peaks in spectrum, too many, aborting" << endl;
	continue;
      }
      h[d][c]->DrawCopy();
      //check if first peak is lower in energy, otherwise swap them
      if(xpeaks[0]>xpeaks[1]){
	Float_t temp = xpeaks[1];
	xpeaks[1] = xpeaks[0];
	xpeaks[0] = temp;
	temp = ypeaks[1];
	ypeaks[1] = ypeaks[0];
	ypeaks[0] = temp;
	
      }

      for(int p=0;p<nfound;p++){
	ca->cd(1+c*4+1+p);
	h[d][c]->GetXaxis()->SetRangeUser(xpeaks[p]-range,xpeaks[p]+range);
	h[d][c]->DrawCopy();
	fu[d][c][p] = new TF1(Form("fcore_d%d_c%d_p%d",d,c,p),fgammagaussbg,xpeaks[p]-range,xpeaks[p]+range,6);
	fu[d][c][p]->SetLineColor(3);
	fu[d][c][p]->SetLineWidth(1);
	fu[d][c][p]->SetParameter(0,0);//bg const
	fu[d][c][p]->SetParameter(1,0);//bg slope
	fu[d][c][p]->SetParameter(2,h[d][c]->Integral(xpeaks[p]-range,xpeaks[p]+range));//norm
	fu[d][c][p]->SetParameter(3,xpeaks[p]);//mean
	fu[d][c][p]->SetParLimits(3,xpeaks[p]-500,xpeaks[p]+500);//mean
	fu[d][c][p]->SetParameter(4,100);//sigma
	fu[d][c][p]->SetParLimits(4,0.001,1000);//sigma
	fu[d][c][p]->SetParameter(5,h[d][c]->GetBinContent(h[d][c]->FindBin(xpeaks[p]-range)));//step
	
	h[d][c]->Fit(fu[d][c][p],"Rqn");
	fu[d][c][p]->Draw("same");


	fus[d][c][p][0] = new TF1(Form("fcore_d%d_c%d_p%d_bg",d,c,p),fgammabg,xpeaks[p]-range,xpeaks[p]+range,6);
	fus[d][c][p][1] = new TF1(Form("fcore_d%d_c%d_p%d_st",d,c,p),fgammastep,xpeaks[p]-range,xpeaks[p]+range,6);
	fus[d][c][p][2] = new TF1(Form("fcore_d%d_c%d_p%d_ga",d,c,p),fgammagaus,xpeaks[p]-range,xpeaks[p]+range,6);
	

	fus[d][c][p][0]->SetLineColor(5);
	fus[d][c][p][1]->SetLineColor(4);
	fus[d][c][p][2]->SetLineColor(2);
	for(int k=0;k<3;k++){
	  fus[d][c][p][k]->SetLineWidth(1);
	  for(int l=0;l<6;l++)
	    fus[d][c][p][k]->SetParameter(l,fu[d][c][p]->GetParameter(l));
	  fus[d][c][p][k]->Draw("same");
	}

      }//peaks
      //res[d*4+c] = 2.35*fu[d][c][1]->GetParameter(4)*(1332.492-1173.228)/(fu[d][c][1]->GetParameter(3)-fu[d][c][0]->GetParameter(3));

      slope[d*4+c] = (1332.492-1173.228)/(fu[d][c][1]->GetParameter(3)-fu[d][c][0]->GetParameter(3));
      offset[d*4+c] = (1332.492+1173.228)-slope[d*4+c]*(fu[d][c][1]->GetParameter(3)+fu[d][c][0]->GetParameter(3));
      offset[d*4+c]*=0.5;
      
      //cout << fu[d][c][0]->GetParameter(3)*slope[d*4+c] + offset[d*4+c] << "\t" << fu[d][c][1]->GetParameter(3)*slope[d*4+c] + offset[d*4+c] << endl;
      
      res[d*4+c] = fu[d][c][1]->GetParameter(2);
      det[d*4+c] = d*5+c;

      output->SetValue(Form("Slope.d%d.c%d",d,c),slope[d*4+c]);
      output->SetValue(Form("Offset.d%d.c%d",d,c),offset[d*4+c]);


      //cout << det[d*4+c] <<"\t" << res[d*4+c] << endl;
      //cout << fu[d][c][0]->GetParameter(4) <<"\t" << fu[d][c][1]->GetParameter(4) << endl;
      //ca->cd(1+c*4+3);
      TPad *p1 = (TPad *)(ca->cd(1+c*4+3)); 
      p1->SetLogy();
      h[d][c]->GetXaxis()->SetRangeUser((int)xpeaks[0]-range,(int)xpeaks[1]+range);
      h[d][c]->DrawCopy();
      for(int p=0;p<nfound;p++){
	fu[d][c][p]->Draw("same");
	for(int k=0;k<3;k++){
	  fus[d][c][p][k]->Draw("same");
	}
      }

    }//crystals
    if(d==0)
      ca->SaveAs("corefits.ps(");
    else
      ca->SaveAs("corefits.ps");
  }//detectors
  output->SaveLevel(kEnvLocal);

  //TGraph *gres = new TGraph(7*4,det,res);
  TGraph *gres = new TGraph(7*4,slope,offset);
  TCanvas *ca2 = new TCanvas("ca2","ca2",0,0,1400,900);
  ca2->cd();
  gres->Draw("AP");
  gres->SetMarkerStyle(20);
  ca2->SaveAs("corefits.ps)");
}


Double_t fgammagaussbg(Double_t *x, Double_t *par){
  static Float_t sqrt2pi = TMath::Sqrt(2*TMath::Pi()), sqrt2 = TMath::Sqrt(2.);
  Double_t arg;
  Double_t result = par[0] + par[1]*x[0];
  
  Double_t norm  = par[2];
  Double_t mean  = par[3];
  Double_t sigma = par[4];

  Double_t step = par[5];

  arg = (x[0]-mean)/(sqrt2*sigma);
  result += 1/(sqrt2pi*sigma) * norm * exp(-arg*arg);

  result += step/pow(1+exp(sqrt2*arg),2);

  return result;

}
Double_t fgammabg(Double_t *x, Double_t *par){
  Double_t result = par[0] + par[1]*x[0]; 
  return result;

}
Double_t fgammastep(Double_t *x, Double_t *par){
  static Float_t sqrt2pi = TMath::Sqrt(2*TMath::Pi()), sqrt2 = TMath::Sqrt(2.);
  Double_t arg;
  Double_t result = 0;
  
  Double_t norm  = par[2];
  Double_t mean  = par[3];
  Double_t sigma = par[4];

  Double_t step = par[5];
  arg = (x[0]-mean)/(sqrt2*sigma);
  result += step/pow(1+exp(sqrt2*arg),2);

  return result;

}
Double_t fgammagaus(Double_t *x, Double_t *par){
  static Float_t sqrt2pi = TMath::Sqrt(2*TMath::Pi()), sqrt2 = TMath::Sqrt(2.);
  Double_t arg;

  Double_t norm  = par[2];
  Double_t mean  = par[3];
  Double_t sigma = par[4];


  arg = (x[0]-mean)/(sqrt2*sigma);
  Double_t result = 1/(sqrt2pi*sigma) * norm * exp(-arg*arg);

  return result;

}
