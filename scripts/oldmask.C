//#include <iostream>
//#include <iomanip>
//#include <string>
//
//#include "TGraph.h"
//#include "TCanvas.h"
//#include "TPad.h"
//#include "TStyle.h"
//#include "TROOT.h"
//#include "TMath.h"
{
  gROOT->Reset();
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1); 
  float yreal[5];
  float y0[5];
  float y1[5];
  
  //ifstream a("mask.dat");
  //ifstream a("scripts/masklater.dat");
  ifstream a("masklater.dat");
  if(a.bad()){
    cerr<<"Unable to open "<<a<<"!\nexiting ... \n";
    exit(2);
  }
  a.ignore(1000, '\n');
  for(int t=0;t<5;t++){
    a >> yreal[t] >> y0[t] >> y1[t];
    a.ignore(1000, '\n');
  }
  float x[2] = {1700,3300};
  float y[2] = {-5,35};
  TGraph *g = new TGraph(2,x,y);
  TCanvas *c = new TCanvas("c","c",0,0,600,400);
  TGraph *g0 = new TGraph(5,y0,yreal);
  g0->SetMarkerColor(2);
  g0->SetMarkerStyle(20);
  
  TGraph *g1 = new TGraph(5,y1,yreal);
  g1->SetMarkerColor(3);
  g1->SetMarkerStyle(20);
  c->cd();
  g->Draw("AP");
  g->GetXaxis()->SetRangeUser(1800,3200);
  g0->Draw("P");
  g0->Fit("pol1");
  g1->Draw("P");
  g1->Fit("pol1");
}
