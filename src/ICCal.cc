#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>

#include "TFile.h"
#include "TString.h"
#include "TKey.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TMath.h"
#include "TCanvas.h"
#include "Peaks.hh"
#include "CommandLineInterface.hh"
using namespace std;
int main(int argc, char* argv[]){
  //vector<int> range;
  char* InputFile = NULL;
  char* OutputFile = NULL;
  int nmatch =-1;
  int nblobs =2;
  CommandLineInterface* interface = new CommandLineInterface();

  //The commands to be read from the command line.
  interface->Add("-i", "inputfile", &InputFile);
  interface->Add("-o", "outputfile", &OutputFile);
  interface->Add("-b", "how many blobs default is 2", &nblobs);
  interface->Add("-m", "channel for matching", &nmatch);
  interface->CheckFlags(argc, argv);

  if(InputFile == NULL || OutputFile == NULL ){
    cerr<<"You have to provide both the input file and the output file!"<<endl;
    exit(1);
  }
  if(nblobs<2){
    cerr<<"The number of blobs has to be 2 or greater!"<<endl;
    exit(1);
  }

  vector<char*> name;
  vector<vector<double> > mean;
  vector<vector<TH1F*> > h;
  vector<vector<TF1*> > f;

    
  name.resize(nblobs);
  mean.resize(nblobs);
  h.resize(nblobs);
  f.resize(nblobs);

  cout<<"input file: "<<InputFile<< endl;
  TFile* File = new TFile(InputFile);
  TIter next(File->GetListOfKeys());
  TKey *key;
  char* Name = NULL;
  int ctr = 0;

  //Loop over the histograms in the input file to find the name of each blob.
  //Note that these cannot have any underscores in the cut name, otherwise the string manipulation fails.
  //It matches "ic_ch0_.*_[^_]*" instead of "ic_ch0_.*"
  while ((key=(TKey*)next())){
    if(strcmp(key->GetClassName(),"TH1F") == 0 && strstr((char*)key->GetName(), "ic_ch0_")){
      Name = (char*)key->GetName();
      int len = strlen(Name);
      int pos = len - 1;
      while(pos >= 0){
	if(Name[pos] == '_'){
	  name[ctr] = (Name + pos+1);
	  mean[ctr].resize(16);
	  h[ctr].resize(16);
	  f[ctr].resize(16);
	  for(int c=0;c<16;c++){
	    f[ctr][c] = new TF1(Form("f_%d_%d",ctr,c),fonegaus,0,2000,3);
	    f[ctr][c]->SetLineColor(3);
	    f[ctr][c]->SetLineWidth(1);
	  }

	  ctr++;
	  break;
	}
	pos--;
      }
      cout << "name " << name[ctr-1] << endl;
      if(ctr == nblobs)
	break;
    }
  }
  ofstream output(Form("%s.dat",OutputFile));
  if(nmatch<0)
    output << "name\t\t\tcont\tmean\tfwhm\tchi" << endl;
  TCanvas *c = new TCanvas("fitresult","fitresult",600,400);
  vector<TH1F*> htest;
  vector<int> max;
  TFile *fout = new TFile(Form("%s.root",OutputFile),"RECREATE");
  //Make the canvases in preparation for showing the fits.
  if(nblobs>0 && nblobs<7){
    if(nblobs==2)
      c->Divide(2,1);
    else if(nblobs<5)
      c->Divide(2,2);
    else 
      c->Divide(3,2);

    htest.resize(nblobs);
    max.resize(nblobs);
    for(int b=0;b<nblobs;b++){
      max[b] =0;
      htest[b] = new TH1F(Form("blob_%d",b),Form("blob_%d",b),2500,0,2500);
      c->cd(1+b);
      htest[b]->Draw();
    }
  }

  for(UShort_t blob=0;blob<nblobs;blob++){
    for(UShort_t ch=0;ch<16;ch++){
      h[blob][ch] = (TH1F*) File->Get(Form("ic_ch%d_%s",ch,name[blob]));
      if(h[blob][ch]==NULL || h[blob][ch]->Integral()==0){
	mean[blob][ch] =0;
	continue;
      }
      else{
	f[blob][ch]->SetParameter(0,h[blob][ch]->Integral());
	f[blob][ch]->SetParameter(1,h[blob][ch]->GetMean());
	f[blob][ch]->SetParameter(2,h[blob][ch]->GetRMS());
	
	h[blob][ch]->Fit(f[blob][ch],"q0");
	mean[blob][ch] = f[blob][ch]->GetParameter(1);
	
	if(nblobs<7){
	  c->cd(blob+1);
	  if(h[blob][ch]->GetMaximum()>max[blob])
	    max[blob] = h[blob][ch]->GetMaximum();
	  h[blob][ch]->SetLineColor(ch+1);
	  f[blob][ch]->SetLineColor(2);
	  h[blob][ch]->Draw("same");
	  f[blob][ch]->Draw("same");
	  fout->cd();
	  h[blob][ch]->Write();
	  f[blob][ch]->Write();
	}
	
	
	if(nmatch<0)
	  output << Form("ic_ch%d_%s",ch,name[blob]) <<"\t" << f[blob][ch]->GetParameter(0) << "\t" << f[blob][ch]->GetParameter(1) << "\t" << f[blob][ch]->GetParameter(2) << "\t" << f[blob][ch]->GetChisquare() << endl;
      }
    }//channels
  }//blobs
  
  if(nblobs<7){
    for(int b=0;b<nblobs;b++)
      htest[b]->GetYaxis()->SetRangeUser(0,max[b]*1.1);
  }
  c->SaveAs(Form("%s.pdf",OutputFile));
  if(nmatch<0)
    return 1;
  double a[16], b[16];
  double sumx =0;
  double sumxs =0;
  double sumy =0;
  for(int c=0;c<16;c++){
    a[c]=0;
    sumx =0;
    sumxs =0;
    sumy =0;
    for(int bl=0;bl<nblobs;bl++){
      a[c]+=mean[bl][c]*mean[bl][nmatch];
      sumx +=mean[bl][c];
      sumy +=mean[bl][nmatch];
      sumxs +=mean[bl][c]*mean[bl][c];
    }
    if(sumx>0&&sumy>0&&sumxs>0){
      a[c] = a[c] - sumx*sumy/nblobs;
      a[c] = a[c]/(sumxs-sumx*sumx/nblobs);
      b[c] = sumy - a[c]*sumx;
      b[c] /= nblobs;
    }
    else{
	a[c] = 0;
	b[c] = 0;
    }
    //cout << "a " << a[c] << "\tb " << b[c] << endl;
    output << "IonChamber.Offset."<<c<<":\t"<<b[c]<<endl;
    output << "IonChamber.Slope."<<c<<":\t"<<a[c]<<endl;
  }
  output << "IonChamber.Offset.DE:\t0.0"<<endl;
  output << "IonChamber.Slope.DE:\t1.0"<<endl;
  cout << "done with fitting " << nblobs << " matching channel "<<nmatch<<endl << " output written to " << Form("%s.dat",OutputFile) << endl;
  fout->Close();
  return 0 ;
}
