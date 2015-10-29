#include "TFile.h"
#include "TH1F.h"
#include "TF1.h"
#include "TH2F.h"
#include "TH2S.h"
#include "TH1S.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TMath.h"
#include "TLine.h"
#include "TCutG.h"
using namespace TMath;
using namespace std;

void PressEnterToContinue(){
  int c;
  printf( "Press ENTER to continue... " );
  fflush( stdout );
  do c = getchar(); 
  while ((c != '\n') && (c != EOF));
}
void plot(char *name, double y0, double y1, bool pos = false){
  TFile *f = new TFile("hist/hist_time_0204.root");
  TH2F* h2d;
  h2d = (TH2F*)f->Get(name);
  if(h2d==NULL)
    return;
  
  h2d->GetYaxis()->SetRangeUser(y0,y1);
  TProfile *hp = (TProfile*)h2d->ProfileX();
  int nx = h2d->GetNbinsX();
  int lx = h2d->GetXaxis()->GetBinLowEdge(1);
  int ux = h2d->GetXaxis()->GetBinUpEdge(h2d->GetNbinsX());
  int ny = h2d->GetNbinsY();
  int ly = h2d->GetYaxis()->GetBinLowEdge(1);
  int uy = h2d->GetYaxis()->GetBinUpEdge(h2d->GetNbinsY());

  cout << nx << "\t" << lx << "\t" << ux << "\t" << h2d->GetXaxis()->GetBinWidth(1) << endl;
  cout << ny << "\t" << ly << "\t" << uy << "\t" << h2d->GetYaxis()->GetBinWidth(1) << endl;


  TH1F * hpp = new TH1F("hpp","hpp",nx,lx,ux);
  TH2F* result = new TH2F("result","resul",nx,lx,ux,ny,ly,uy);

  double setvalue;
  setvalue = hp->GetMean(2);

  //adjust here when your mask runs where!
  double x[3] = {0,47793559,105931206};
  TF1 *fu[3];
  if(pos){
    setvalue =0;
    for(int i=0;i<3;i++){
      fu[i]= new TF1(Form("fu%d",i),"pol0",x[i],x[i]+4e6);
      if(i==1)
	fu[i]= new TF1(Form("fu%d",i),"pol0",x[i]-2e6,x[i]+2e6);
      hp->Fit(fu[i],"Rq");
      cout << fu[i]->GetParameter(0) <<"\t";
      setvalue += fu[i]->GetParameter(0);
    }
    cout << endl;
    setvalue /=3;
  }

  for(int b = 0;b<hpp->GetNbinsX()+1;b++){
  //for(int b = 0;b<100+1;b++){
    hpp->SetBinContent(b,-(hp->GetBinContent(b)- setvalue));
    for(int y = 0;y<h2d->GetNbinsY()+1;y++){
      double yp = result->GetYaxis()->GetBinCenter(y);
      //cout << y <<"\t"<< yp <<"\t";
      yp -=(hp->GetBinContent(b)- setvalue);
      //cout << yp <<"\t" << result->GetYaxis()->FindBin(yp) << endl;
      result->SetBinContent(b,result->GetYaxis()->FindBin(yp),h2d->GetBinContent(b,y));
    }
  }

  TCanvas *c = new TCanvas("c","c",0,0,800,700);
  c->Divide(2,2);
  c->cd(1);
  h2d->RebinX(100);
  h2d->GetYaxis()->SetRangeUser(y0,y1);
  h2d->DrawCopy("colz");
  TLine *l[3];
  if(pos){
    for(int i=0;i<3;i++){
      l[i] = new TLine(x[i],-150,x[i],150);
      l[i]->SetLineColor(2);
      l[i]->Draw();
    }
  }
  c->cd(2);
  hp->Draw();
  if(pos){
    for(int i=0;i<3;i++){
      l[i]->Draw();
      fu[i]->Draw("same");
    }
  }
  c->cd(3);
  hpp->Draw();
  c->cd(4);
  result->RebinX(100);
  result->GetYaxis()->SetRangeUser(y0,y1);
  result->Draw("colz");
  cout << hp->GetMean(1) <<"\t"<< hp->GetMean(2) <<"\t"<< setvalue << endl;


  hpp->SetName(Form("%s_cor",name));
  
  // TFile *corrections = new TFile("corrections.root","UPDATE");
  // hpp->Write("",TObject::kOverwrite);
  // corrections->ls();
  // corrections->Close();
  c->Modified();
  c->Update();
  //PressEnterToContinue();
}
void icplot(double y0, double y1){
  TFile *f = new TFile("hist/hist_time_0204.root");
  //TFile *f = new TFile("hist/hist_test_0128.root");
  TH2F* h2d;
  h2d = (TH2F*)f->Get("IC_vs_time");
  if(h2d==NULL)
    return;

  TFile *cut = new TFile("ICcut.root");
  TCutG *gc = (TCutG*)cut->Get("mainpart");
  TCutG *gc2 = (TCutG*)cut->Get("lowpart");
  
  h2d->GetYaxis()->SetRangeUser(y0,y1);
  TProfile *hp = (TProfile*)h2d->ProfileX("prof",0,750,"[mainpart]");
  TProfile *hp2 = (TProfile*)h2d->ProfileX("prof2",0,750,"[lowpart]");
  int nx = h2d->GetNbinsX();
  int lx = h2d->GetXaxis()->GetBinLowEdge(1);
  int ux = h2d->GetXaxis()->GetBinUpEdge(h2d->GetNbinsX());
  int ny = h2d->GetNbinsY();
  int ly = h2d->GetYaxis()->GetBinLowEdge(1);
  int uy = h2d->GetYaxis()->GetBinUpEdge(h2d->GetNbinsY());

  cout << nx << "\t" << lx << "\t" << ux << "\t" << h2d->GetXaxis()->GetBinWidth(1) << endl;
  cout << ny << "\t" << ly << "\t" << uy << "\t" << h2d->GetYaxis()->GetBinWidth(1) << endl;

  TH1F * hpp = new TH1F("hpp","hpp",nx,lx,ux);
  TH1F * hpp2 = new TH1F("hpp2","hpp2",nx,lx,ux);
  TH2F* result = new TH2F("result","resul",nx,lx,ux,ny,ly,uy);

  double setvalue[2];
  setvalue[0] = hp->GetMean(2);
  setvalue[1] = hp2->GetMean(2);
  cout << hp->GetMean(1) <<"\t"<< hp->GetMean(2) <<"\t"<< setvalue[0] << endl;
  cout << hp2->GetMean(1) <<"\t"<< hp2->GetMean(2) <<"\t"<< setvalue[1] << endl;
  for(int b = 0;b<hpp->GetNbinsX()+1;b++){
  //for(int b = 0;b<100;b++){
    //cout << hp->GetBinContent(b)<<"\t"<<hp2->GetBinContent(b)<<endl;
    double g =0;
    if((hp->GetBinContent(b)-hp2->GetBinContent(b)) !=0)
      g = (setvalue[0]-setvalue[1])/(hp->GetBinContent(b)-hp2->GetBinContent(b));
    double o = (setvalue[0]+setvalue[1])-g*(hp->GetBinContent(b)+hp2->GetBinContent(b));
    o/=2;
    //cout << g<<"\t"<<o<<endl;
    hpp->SetBinContent(b,g);
    hpp2->SetBinContent(b,o);
    for(int y = 0;y<h2d->GetNbinsY()+1;y++){
      double yp = result->GetYaxis()->GetBinCenter(y);
      yp = yp*g+o;
      result->SetBinContent(b,result->GetYaxis()->FindBin(yp),h2d->GetBinContent(b,y));
    }
  }

  TCanvas *c = new TCanvas("c","c",0,0,800,700);
  c->Divide(3,2);
  c->cd(1);
  //h2d->RebinX(100);
  h2d->GetYaxis()->SetRangeUser(y0,y1);
  h2d->DrawCopy("colz");
  gc->Draw();
  gc2->Draw();
  c->cd(2);
  hp->Draw();
  c->cd(3);
  hp2->Draw();
  c->cd(4);
  hpp->Draw();
  c->cd(5);
  hpp2->Draw();
  c->cd(6);
  //result->RebinX(100);
  result->GetYaxis()->SetRangeUser(y0,y1);
  result->Draw("colz");
 

  hpp->SetName("IC_vs_time_cor");
  
  // TFile *corrections = new TFile("corrections.root","UPDATE");
  // hpp->Write("",TObject::kOverwrite);
  // corrections->ls();
  // corrections->Close();
  c->Modified();
  c->Update();
  //PressEnterToContinue();
}
void noplot(char *name, double y0, double y1, bool pos = false){
  TFile *f = new TFile("hist/hist_time_0204.root");
  TH2F* h2d;
  h2d = (TH2F*)f->Get(name);
  if(h2d==NULL)
    return;
  
  h2d->GetYaxis()->SetRangeUser(y0,y1);
  TProfile *hp = (TProfile*)h2d->ProfileX();
  int nx = h2d->GetNbinsX();
  int lx = h2d->GetXaxis()->GetBinLowEdge(1);
  int ux = h2d->GetXaxis()->GetBinUpEdge(h2d->GetNbinsX());
  int ny = h2d->GetNbinsY();
  int ly = h2d->GetYaxis()->GetBinLowEdge(1);
  int uy = h2d->GetYaxis()->GetBinUpEdge(h2d->GetNbinsY());


  TH1F * hpp = new TH1F("hpp","hpp",nx,lx,ux);
 
  double setvalue;
  setvalue = hp->GetMean(2);
  //adjust here when your mask runs where!
  double x[3] = {0,47793559,105931206};
  TF1 *fu[3];
  if(pos){
    setvalue =0;
    for(int i=0;i<3;i++){
      fu[i]= new TF1(Form("fu%d",i),"pol0",x[i],x[i]+4e6);
      if(i==1)
	fu[i]= new TF1(Form("fu%d",i),"pol0",x[i]-2e6,x[i]+2e6);
      hp->Fit(fu[i],"Rq");
      cout << fu[i]->GetParameter(0) <<"\t";
      setvalue += fu[i]->GetParameter(0);
    }
    cout << endl;
    setvalue /=3;
  }

  for(int b = 0;b<hpp->GetNbinsX()+1;b++){
    hpp->SetBinContent(b,-(hp->GetBinContent(b)- setvalue));
  }

  hpp->SetName(Form("%s_cor",name));
  
  TFile *corrections = new TFile("corrections.root","UPDATE");
  hpp->Write("",TObject::kOverwrite);
  corrections->ls();
  corrections->Close();
}
void noicplot(double y0, double y1){
  TFile *f = new TFile("hist/hist_time_0204.root");
  TH2F* h2d;
  h2d = (TH2F*)f->Get("IC_vs_time");
  if(h2d==NULL)
    return;

  TFile *cut = new TFile("ICcut.root");
  TCutG *gc = (TCutG*)cut->Get("mainpart");
  TCutG *gc2 = (TCutG*)cut->Get("lowpart");
  
  h2d->GetYaxis()->SetRangeUser(y0,y1);
  TProfile *hp = (TProfile*)h2d->ProfileX("prof",0,750,"[mainpart]");
  TProfile *hp2 = (TProfile*)h2d->ProfileX("prof2",0,750,"[lowpart]");
  int nx = h2d->GetNbinsX();
  int lx = h2d->GetXaxis()->GetBinLowEdge(1);
  int ux = h2d->GetXaxis()->GetBinUpEdge(h2d->GetNbinsX());
  int ny = h2d->GetNbinsY();
  int ly = h2d->GetYaxis()->GetBinLowEdge(1);
  int uy = h2d->GetYaxis()->GetBinUpEdge(h2d->GetNbinsY());

  cout << nx << "\t" << lx << "\t" << ux << "\t" << h2d->GetXaxis()->GetBinWidth(1) << endl;
  cout << ny << "\t" << ly << "\t" << uy << "\t" << h2d->GetYaxis()->GetBinWidth(1) << endl;

  TH1F * hpp = new TH1F("hpp","hpp",nx,lx,ux);
  TH1F * hpp2 = new TH1F("hpp2","hpp2",nx,lx,ux);
  hpp->SetName("IC_vs_time_corg");
  hpp2->SetName("IC_vs_time_coro");

  double setvalue[2];
  setvalue[0] = hp->GetMean(2);
  setvalue[1] = hp2->GetMean(2);
  for(int b = 0;b<hpp->GetNbinsX()+1;b++){
    double g =0;
    if((hp->GetBinContent(b)-hp2->GetBinContent(b)) !=0)
      g = (setvalue[0]-setvalue[1])/(hp->GetBinContent(b)-hp2->GetBinContent(b));
    double o = (setvalue[0]+setvalue[1])-g*(hp->GetBinContent(b)+hp2->GetBinContent(b));
    o/=2;
    hpp->SetBinContent(b,g);
    hpp2->SetBinContent(b,o);
  }

  TFile *corrections = new TFile("corrections.root","UPDATE");
  hpp->Write("",TObject::kOverwrite);
  hpp2->Write("",TObject::kOverwrite);
  corrections->ls();
  corrections->Close();

}
void run(){
  noicplot(450,800);
  noplot("obj_vs_time",-1350,-1260);
  noplot("xfp_vs_time",-920,-800);
  noplot("xfptac_vs_time",2780,2820);
  noplot("objtac_vs_time",2060,2130);
  noplot("y0_vs_time",-150,150,true);
  noplot("y1_vs_time",-150,150,true);
}
