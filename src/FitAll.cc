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
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#include "TStyle.h"
#include "TROOT.h"

#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TGraphErrors.h"
#include "TLine.h"
#include "TF1.h"
#include "TEnv.h"
#include "TMath.h"
#include "TCanvas.h"
#include "Peaks.hh"
#include "CommandLineInterface.hh"
using namespace std;
int alphanum(int b, int c){
  int a = -1;
  int s = -1;
  if(b==0||b==2){
    if(c<6){
      s = c;
      a = b/2*3;
    }
    else{
      s = c-6;
      a = b/2*3+1;
    }
  }
  if(b==1){
    if(c<3){
      a = 1;
      s = c+3;
    }
    else{
      a = 2;
      s = c-3;
    }
  }
  if(b==3){
    if(c<3){
      a = 4;
      s = c+3;
    }
    else{
      a = 5;
      s = c-3;
    }
  }
  
  
  return a*10+s;
}
int Q1map(int b, int c, int cr){
  if(cr==1){
    if(b==2)
      return alphanum(3,c);
    if(b==3)
      return alphanum(2,c);
    return alphanum(b,c);
  }
  if(cr==3)
    return alphanum(b,c);

  int map[4][9] = {{1,  8,  3,  4, 35,  0,  7, 14,  9}, 
		   {10, 29, 30, 13, 20, 15, 16, 23, 24}, 
		   {19, 26, 21, 22, 17, 18, 25, 32, 27},
		   {28, 11, 12, 31,  2, 33, 34,  5,  6}};

  if(cr==0)
    return map[b][c]/6*10 + map[b][c]%6;

  if(cr==2){
    if(b==2)
      return map[3][c]/6*10 + map[3][c]%6;
    if(b==3)
      return map[2][c]/6*10 + map[2][c]%6;
    return map[b][c]/6*10 + map[b][c]%6;
  }
  return -1;
}



int main(int argc, char* argv[]){
  //vector<int> range;
  char* InputFile = NULL;
  char* OutputFile = NULL;
  CommandLineInterface* interface = new CommandLineInterface();

  interface->Add("-i", "inputfile", &InputFile);
  interface->Add("-o", "outputfile", &OutputFile);
  interface->CheckFlags(argc, argv);

  if(InputFile == NULL || OutputFile == NULL ){
    cerr<<"You have to provide both the input file and the output file!"<<endl;
    exit(1);
  }
  gROOT->Reset();
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1); 
  TFile* File = new TFile(InputFile);
  TFile* HFile = new TFile(Form("%s.root",OutputFile),"RECREATE");
  //TEnv* output = new TEnv(Form("%s.dat",OutputFile));
  TH1S* h[4][10][4][7];//board 4, ch 10, cry 4, detector 7


  TH1F* hres[7][4];// det cry

  TCanvas * ca= new TCanvas("ca","ca",800,600);
  TCanvas * cs= new TCanvas("cs","cs",800,600);
  int ctr =0, ectr=0;
  int crangelow[4] = {2e5,6e4,12e5,4e5};
  int crangehig[4] = {5e5,18e4,20e5,10e5};
  int srangelow = 3e5;
  int srangehig = 6e5;
  
  TF1 *func[2];
  for(int d=0;d<7;d++){
    ca->Divide(4,4);
    for(int cr=0;cr<4;cr++){
      hres[d][cr] = new TH1F(Form("hres_d%d_c%d",d+1,cr),Form("hres_d%d_c%d",d+1,cr),50,0,50);
      cs->Divide(6,6);
      //cout << "----------------------------------------------" << endl;
      cout << "detector " << d << " crystal " << cr << endl; 
      for(int b=0;b<4;b++){
	for(int c=0;c<10;c++){
	  
	  h[b][c][cr][d] = (TH1S*)File->Get(Form("henergy_b%d_c%d_cr%d_d%d",b,c,cr,d+1));
	  if(h[b][c][cr][d]->Integral()>0){
	    ctr++;
	  }
	  else{
	    //cout << Form("henergy_b%d_c%d_cr%d_d%d",b,c,cr,d+1) << endl;
	    ectr++;
	    continue;
	  }
	  if(c!=9){
	    int a = alphanum(b,c)/10;
	    int s = alphanum(b,c)%10;
	    if(d==1){ //Q2 adjust here
	      if(alphanum(b,c)==24)
		a = 3;
	      if(alphanum(b,c)==34)
		a = 2;
	    }
	    if(d==0){ //Q1 adjust here
	      a = Q1map(b,c,cr)/10;
	      s = Q1map(b,c,cr)%10;
	    }
	    //cout << "b " << b << " c " << c << "\ta " << a << " s " << s << endl;
	    cs->cd(a*6+s+1);
	    h[b][c][cr][d]->GetXaxis()->SetRangeUser(srangelow,srangehig);
	    TSpectrum *sp = new TSpectrum(2,5);
	    sp->SetResolution(5);
	    Int_t nfound = sp->Search(h[b][c][cr][d],5,"nobackground",0.1);
	    //cout << "Found " << nfound << " peaks in spectrum" << endl;
	    Float_t *xpeaks = sp->GetPositionX();
	    if(nfound!=2){
	      cout << "problem in " << Form("henergy_b%d_c%d_cr%d_d%d",b,c,cr,d+1) << endl;
	      continue;
	    }
	    
	    for(int f=0;f<nfound;f++){
	      
	      Float_t xp = xpeaks[f];
	      Int_t bin = h[b][c][cr][d]->GetXaxis()->FindBin(xp);
	      Float_t yp = h[b][c][cr][d]->GetBinContent(bin);
	      
	      func[f] = new TF1(Form("fseg_%d_b%d_c%d_cr%d_d%d",f,b,c,cr,d+1),fonegausbg,xp-10000,xp+10000,5);
	      func[f]->SetLineColor(3);
	      func[f]->SetLineWidth(1);
	      func[f]->SetParameter(0,(h[b][c][cr][d]->GetBinContent((int)xp-10000)+ h[b][c][cr][d]->GetBinContent((int)xp+10000) )/2 );
	      func[f]->SetParameter(1,0);
	      func[f]->SetParameter(2,yp*1.4*500);
	      func[f]->SetParameter(3,xp);
	      func[f]->SetParameter(4,500);
	      h[b][c][cr][d]->Fit(func[f],"q0R");
	      HFile->cd();
	      func[f]->Write();
	      
	    }
	    h[b][c][cr][d]->Draw();
	    h[b][c][cr][d]->GetXaxis()->SetRangeUser(srangelow,srangehig);
	    HFile->cd();
	    h[b][c][cr][d]->Write();
	    func[0]->Draw("same");
	    func[1]->Draw("same");
	    double res =1;
 	    if(func[0]->GetParameter(3)<func[1]->GetParameter(3))
 	      res = func[1]->GetParameter(4)*(1332.492-1173.228)/(func[1]->GetParameter(3)-func[0]->GetParameter(3));
 	    else
	      res = func[0]->GetParameter(4)*(1332.492-1173.228)/(func[0]->GetParameter(3)-func[1]->GetParameter(3));

	    hres[d][cr]->SetBinContent(1+s+7*a,res); 
	    
	    
	    
	    h[b][c][cr][d]->Draw();
	  }
	}//chans

	ca->cd(1+cr*4+b);
	if(h[b][9][cr][d]->Integral()>0){
	  h[b][9][cr][d]->GetXaxis()->SetRangeUser(crangelow[b],crangehig[b]);
	  if(d==4 && cr==3 && b==1)
	    h[b][9][cr][d]->GetXaxis()->SetRangeUser(crangelow[0],crangehig[0]);
	    
	  TSpectrum *sp = new TSpectrum(2,5);
	  sp->SetResolution(5);
	  Int_t nfound = sp->Search(h[b][9][cr][d],5,"nobackground",0.5);
	  //cout << "Found " << nfound << " peaks in spectrum" << endl;
	  Float_t *xpeaks = sp->GetPositionX();
	  if(nfound!=2){
	    cout << "problem in " << Form("henergy_b%d_c%d_cr%d_d%d",b,9,cr,d+1) << endl;
	    continue;
	  }
	  
	  for(int f=0;f<nfound;f++){
	    
	    Float_t xp = xpeaks[f];
	    Int_t bin = h[b][9][cr][d]->GetXaxis()->FindBin(xp);
	    Float_t yp = h[b][9][cr][d]->GetBinContent(bin);
	    
	    func[f] = new TF1(Form("fcore_%d_b%d_cr%d_d%d",f,b,cr,d+1),fonegausbg,xp-10000,xp+10000,5);
	    func[f]->SetLineColor(3);
	    func[f]->SetLineWidth(1);
	    func[f]->SetParameter(0,(h[b][9][cr][d]->GetBinContent((int)xp-10000)+ h[b][9][cr][d]->GetBinContent((int)xp+10000) )/2 );
	    func[f]->SetParameter(1,0);
	    func[f]->SetParameter(2,yp*1.4*500);
	    func[f]->SetParameter(3,xp);
	    func[f]->SetParameter(4,500);
	    h[b][9][cr][d]->Fit(func[f],"q0R");
	    HFile->cd();
	    func[f]->Write();
	  }
	  h[b][9][cr][d]->Draw();
	  h[b][9][cr][d]->GetXaxis()->SetRangeUser(crangelow[b],crangehig[b]);
	  if(d==4 && cr==3 && b==1)
	    h[b][9][cr][d]->GetXaxis()->SetRangeUser(crangelow[0],crangehig[0]);
	  HFile->cd();
	  h[b][9][cr][d]->Write();
	  func[0]->Draw("same");
	  func[1]->Draw("same");
	  double res =0;
	  if(func[0]->GetParameter(3)<func[1]->GetParameter(3))
	    res = func[1]->GetParameter(4)*(1332.492-1173.228)/(func[1]->GetParameter(3)-func[0]->GetParameter(3));
	  else
	    res = func[0]->GetParameter(4)*(1332.492-1173.228)/(func[0]->GetParameter(3)-func[1]->GetParameter(3));
	  hres[d][cr]->SetBinContent(1+7*6+4+b,res); 
	}
      }//boards  
      if(d==0&&cr==0)
	cs->SaveAs(Form("%s_seg.ps(",OutputFile));
      else if(d==6&&cr==3)
	cs->SaveAs(Form("%s_seg.ps)",OutputFile));
      else
	cs->SaveAs(Form("%s_seg.ps",OutputFile));
      cs->Clear();
      
       HFile->cd();
       hres[d][cr]->Write();
    }//cyrstals
    if(d==0)
      ca->SaveAs(Form("%s.ps(",OutputFile));
    else if(d==6)
      ca->SaveAs(Form("%s.ps)",OutputFile));
    else
      ca->SaveAs(Form("%s.ps",OutputFile));
    ca->Clear();
  }//detectors	 
  cout << ctr << "\t" << ectr<< endl;
  
  return 0;
}
