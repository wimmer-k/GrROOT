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

#include "TFile.h"
#include "TString.h"
#include "TKey.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TGraph.h"
#include "TF1.h"
#include "TEnv.h"
#include "TMath.h"
#include "TCanvas.h"
#include "Peaks.hh"
#include "CommandLineInterface.hh"
using namespace std;

double GetFitMean(TH1F* hst){
  if(hst->Integral()<2){
    return sqrt(-1);
  }
  TF1* fit = new TF1("fit",fonegaus,0,1500,3);
  fit->SetLineColor(3);
  fit->SetLineWidth(1);

  fit->SetParameter(0,hst->Integral());

  double center;
  if(hst->Integral()>100){
    center = hst->GetMaximumBin();
  } else {
    TH1F* temp = new TH1F(*hst);
    temp->Rebin(16);
    center = 16*temp->GetMaximumBin();
  }
  fit->SetParameter(1,center);
  fit->SetRange(center-25,center+25);
  fit->SetParLimits(1,center-25,center+25);

  fit->SetParameter(2,hst->GetRMS());
  hst->Fit(fit,"RqL");

  return fit->GetParameter(1);
}

int main(int argc, char* argv[]){
  //vector<int> range;
  char* InputFile = NULL;
  char* OutputFile = NULL;
  int HOutput = 0;
  int nmatch =-1;
  int niter =1;
  const int blobs =4;  //can be adjusted, depending on the case you might not have 4 good isotopes that are spread over the whole focal plane.
  CommandLineInterface* interface = new CommandLineInterface();

  interface->Add("-i", "inputfile", &InputFile);
  interface->Add("-o", "outputfile", &OutputFile);
  interface->Add("-n", "iteration n=0 new calfile is created, != 0 old calibration is updated", &niter);
  interface->Add("-m", "channel for matching", &nmatch);
  interface->Add("-h", "also produce histogramm output file", &HOutput);
  interface->CheckFlags(argc, argv);

  if(InputFile == NULL || OutputFile == NULL ){
    cerr<<"You have to provide both the input file and the output file!"<<endl;
    exit(1);
  }
  if(HOutput == 0){
    cout<<"No Histogram file created."<<endl;
  }

  TFile* File = new TFile(InputFile);

  //Suffixes of spectra to be used in calibration.
  string hname[blobs];
  string Name;
  TIter next(File->GetListOfKeys());
  TKey *key;
  int namesfound = 0;
  string prefix = Form("crdcpad0_0_%d_",nmatch);
  while((key=(TKey*)next()) && namesfound<blobs){
    Name = key->GetName();
    //Found a spectrum with the right name.
    if(Name.compare(0,prefix.length(),prefix)==0){
      string newblob = Name.substr(prefix.length());
      //Check if I've already found that name.
      bool alreadyUsed = false;
      for(int i=0; i<blobs; i++){
	if(hname[i] == newblob){
	  alreadyUsed = true;
	}
      }
      if(!alreadyUsed && newblob!="all" && newblob.find("coinc")==string::npos){
	cout << "Found outgoing cut: " << newblob << endl;
	hname[namesfound] = newblob;
	namesfound++;
      }
    }
  }
  if(namesfound<blobs){
    cout << "Histogram file must have at least "<<blobs<<" outgoing gates" << endl;
    cout << "change line 56\"   const int blobs =4; \" to your number of outgoing gates"<<endl;
    exit(2);
  }

  //Set the initial gains and offsets as determined in the previous iteration.
  //If it is the first iteration, assume 0 offset and 1 gain.
  double mean[2][224][blobs];
  double gain[2][224];
  double offset[2][224];
  cout<<"input file: "<<InputFile<< endl;
  TEnv* output = new TEnv(Form("%s.dat",OutputFile));
  for(int c=0;c<2;c++){
    for(int p=0;p<224;p++){
      if(niter==0){
	gain[c][p] = 1.0;
	offset[c][p] = 0.0;
      }
      else{
	gain[c][p] = output->GetValue(Form("Crdc.%d.Slope.%03d",c,p),1.0);
	offset[c][p] = output->GetValue(Form("Crdc.%d.Offset.%03d",c,p),0.0);
      }
    }
  }

  //Prepare the output file, if requested.
  TFile* HFile = NULL;
  if(HOutput)
    HFile = new TFile(Form("%s.root",OutputFile),"RECREATE");


  TH1F *h[blobs];
  TF1 *f;
  f = new TF1("f",fonegaus,0,1500,3);
  f->SetLineColor(3);
  f->SetLineWidth(1);

  //For each crdc
  for(int c=0;c<2;c++){
    //Grab the histograms to be used in calibrating the matching channel.
    for(int hi=0;hi<blobs;hi++){
      if(niter<2)
	h[hi]=(TH1F*) File->Get(Form("crdcpad0_%d_%d_%s",c,nmatch,hname[hi].c_str()));
      else
	h[hi]=(TH1F*) File->Get(Form("crdcpad_%d_%d_%s",c,nmatch,hname[hi].c_str()));
      if(h[hi]->Integral()==0){
	cout << "error matching channel not found in all blobs" << endl;
	cout << "missing in " << hname[hi].c_str() << endl;
	return 7;
      }
    }

    //Perform the fit for the matching channel.
    cout << "Matching channel: ";
    for(int i=0;i<blobs;i++){
      mean[c][nmatch][i] = GetFitMean(h[i]);
      if(HOutput){
	HFile->cd();
	h[i]->Write();
      }
      cout << mean[c][nmatch][i] << "\t";
    }
    cout << endl;

    //For each other channel
    for(int p=0;p<224;p++){
      if(p==nmatch)
	continue;

      //Grab the histograms.
      for(int hi=0;hi<blobs;hi++){
	if(niter<2)
	  h[hi]=(TH1F*) File->Get(Form("crdcpad0_%d_%d_%s",c,p,hname[hi].c_str()));
	else
	  h[hi]=(TH1F*) File->Get(Form("crdcpad_%d_%d_%s",c,p,hname[hi].c_str()));
      }

      int goodPeaks = 0;
      for(int i=0;i<blobs;i++){
	if(h[i]==NULL){
	  //This spectrum is empty, don't bother calibrating.
	  mean[c][p][i] = 0;
	  continue;
	} else {
	  goodPeaks++;
	  if(HOutput){
	    HFile->cd();
	    h[i]->Write();
	  }
	  mean[c][p][i] = GetFitMean(h[i]);
	}
      }

      if(goodPeaks==blobs || goodPeaks==blobs-1){
	//Make a linear fit between the current channel and the matching channel.
	//Use this to determine the new slope and offset.
	TGraph* g = new TGraph(goodPeaks,mean[c][p],mean[c][nmatch]);
	g->SetName(Form("line_%d_%d",c,p));
	//Some of the pads have very bad statistics, and so the fit can give negative gains.
	//Using "[0]+[1]*x" instead of "pol1" so that I can specify a limit on the gain.
	TF1* line = new TF1("line","[0]+[1]*x");
	line->SetParameter(0,0.0);
	line->SetParameter(1,1.0);
	line->SetParLimits(1,0.1,2.0);
	g->Fit(line,"q");
	double newoffset = offset[c][p]*line->GetParameter(1) + line->GetParameter(0);
	double newgain = gain[c][p]*line->GetParameter(1);

	cout << "CRDC: " << c << "\tPAD: " << p << endl
	     << "par0 " << line->GetParameter(0) << " par1 " << line->GetParameter(1) << endl
	     << "old offset: " << offset[c][p]
	     << "\tnew offset: " << newoffset << endl
	     << "old gain: " << gain[c][p]
	     << "\tnew gain: " << newgain << endl;

	if(newoffset>-250 && newoffset<250 &&
	   newgain>0.5 && newgain<1.5){
	  offset[c][p] = newoffset;
	  gain[c][p] = newgain;
	  output->SetValue(Form("Crdc.%d.Slope.%03d",c,p),gain[c][p]);
	  output->SetValue(Form("Crdc.%d.Offset.%03d",c,p),offset[c][p]);
	} else {
	  cout << "New gain/offset are unreasonable, keeping old values" << endl;
	  output->SetValue(Form("Crdc.%d.Slope.%03d",c,p),Form("%f # %f",gain[c][p],newgain));
	  output->SetValue(Form("Crdc.%d.Offset.%03d",c,p),Form("%f # %f",offset[c][p],newoffset));
	}

	cout << endl;

	if(HOutput){
	  HFile->cd();
	  g->Write();
	}
      }
    }//pad
  }//crdc


  output->SaveLevel(kEnvLocal);
  if(HOutput){
    HFile->Close();
  }
  return 0;

}
