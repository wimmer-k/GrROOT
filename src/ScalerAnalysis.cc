#include <iostream>
#include <iomanip>
#include <string.h>
#include <sys/time.h>
#include "Riostream.h"

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TMath.h"

#include "TKey.h"

#include "CommandLineInterface.hh"
#include "Scaler.hh"
#include "Scalerdefs.h"

using namespace TMath;
using namespace std;

//ClassImp(Pad)

double get_time(){
    struct timeval t;
    gettimeofday(&t, NULL);
    double d = t.tv_sec + (double) t.tv_usec/1000000;
    return d;
}

void removeDot(string &str){
  std::string temp;
  for(unsigned int i=0; i<str.length(); i++){
    if(str[i] == '.')
      temp += '_';
    else
      temp += str[i];
  }
  str = temp;
}

int main(int argc, char* argv[]){
  double time_start = get_time();
  vector<char*> InputFiles;
  char* OutputFile = NULL;
  char* SetFile = NULL;

  CommandLineInterface* interface = new CommandLineInterface();

  //The settings file here refers to the file relating scaler number and scaler name.
  //This is NOT the settings file used in calibration.
  interface->Add("-i", "inputfiles", &InputFiles);
  interface->Add("-o", "outputfile", &OutputFile);
  interface->Add("-s", "settingsfile", &SetFile);
  interface->CheckFlags(argc, argv);

  if(InputFiles.size() == 0 || OutputFile == NULL){
    cerr<<"You have to provide at least one input file and the output file!"<<endl;
    exit(1);
  }
  if(SetFile==NULL){
    exit(2);
  }


  //Open the input files, build the TChain that will be looped over.
  //Open the output file.
  cout<<"input file(s):"<<endl;
  for(unsigned int i=0; i<InputFiles.size(); i++){
    cout<<InputFiles[i]<<endl;
  }
  cout<<"output file: "<<OutputFile<< endl;
  TChain* tr;
  tr = new TChain("sc");
  for(unsigned int i=0; i<InputFiles.size(); i++){
    tr->Add(InputFiles[i]);
  }
  if(tr == NULL){
    cout << "could not find tree sc in file " << endl;
    for(unsigned int i=0; i<InputFiles.size(); i++){
      cout<<InputFiles[i]<<endl;
    }
    return 3;
  }
  Scaler *sc = new Scaler;
  tr->SetBranchAddress("scaler",&sc);
  TFile* outfile = new TFile(OutputFile,"recreate");
  if(outfile->IsZombie()){
    return 4;
  }


  cout << "creating histograms " << endl;
  Double_t nentries = tr->GetEntries();
  ifstream file(SetFile, ios::in);
  char temp[900];
  int nr;
  char title[80];
  vector<TGraph*> hrate;
  vector<TGraph*> hvalue;
  vector<long long> previous;
  vector<vector<long long> > addme;
  ULong64_t previous_ts = 0;
  vector<int> scalerID;
  int scalers = 0;

  //files
  addme.resize(InputFiles.size());

  //Read through the entirety of the settings file.
  //Make graphss for each of the scalers found.
  while(!file.eof()){
    file.getline(temp,900);
    if (temp[0]=='#'){
      continue;
    }
    int filled = sscanf(temp, "%s %d",title,&nr);
    if (filled!=2){
      continue;
    }
    string stitle = title;
    removeDot(stitle);
    scalerID.push_back(nr);
    TGraph* newRate = new TGraph(nentries);
    newRate->SetName(Form("rate_%s",stitle.c_str()));
    newRate->SetTitle(Form("rate_%s",stitle.c_str()));
    hrate.push_back(newRate);
    TGraph* newVal = new TGraph(nentries);
    newVal->SetName(Form("value_%s",stitle.c_str()));
    newVal->SetTitle(Form("value_%s",stitle.c_str()));
    hvalue.push_back(newVal);
    previous.push_back(0);
    scalers++;
    for(unsigned int i=0; i<InputFiles.size(); i++){
      addme[i].push_back(0);
    }


    cout << "Found scaler " << stitle << " at channel " << nr << endl;
  }
  TGraph* ratio[3];
  ratio[0] = new TGraph(nentries);
  ratio[0]->SetName("ratio_Trigger");
  ratio[0]->SetTitle("ratio_Trigger");
  ratio[1] = new TGraph(nentries);
  ratio[1]->SetName("ratio_Clock");
  ratio[1]->SetTitle("ratio_Clock");
  ratio[2] = new TGraph(nentries);
  ratio[2]->SetName("tsdiff");
  ratio[2]->SetTitle("tsdiff");

  Int_t nbytes = 0;
  Int_t status;

  int filectr =0;

  cout << nentries << " entries in tree " << endl;
  for(int i=0; i<nentries;i++){
    status = tr->GetEvent(i);
    if(status == -1){
      cerr<<"Error occured, couldn't read entry "<<i<<" from tree "<<tr->GetName()<<" in file "<<tr->GetFile()->GetName()<<endl;
      return 5;
    }
    else if(status == 0){
      cerr<<"Error occured, entry "<<i<<" in tree "<<tr->GetName()<<" in file "<<tr->GetFile()->GetName()<<" doesn't exist"<<endl;
      return 6;
    }
    nbytes += status;

    if(sc->GetTS()<previous_ts){
      cout << "ts\tprev "<< previous_ts << "\tnow " << sc->GetTS() << "\tdiff" << previous_ts-sc->GetTS() << " new file "<< endl;
      filectr++;
    }
    double savelifetimes[2][2] ={ {0,0},{0,0}};
    for(int j=0;j<scalers;j++){
      long long id = scalerID[j];
      long long value = sc->GetValue(id);
      if(sc->GetTS()<previous_ts){
	addme[filectr][j] = previous[j];
	//cout << "addme["<<filectr<<"]["<<j<<"] " << addme[filectr][j] << "\tprevious["<<j<<"] " << previous[j] << endl;
      }
      //if the file changes, scalers are reset, therefore we need to add the final value of the last run
      if(value+addme[filectr][j]-previous[j]>1e7)
	addme[filectr][j] -= pow(2,24);
      value += addme[filectr][j];
      //filling
      hvalue[j]->SetPoint(i,i,value);
      hrate[j]->SetPoint(i,i,value-previous[j]);

      if(j==10)
	savelifetimes[0][0] = value-previous[j];
      if(j==9)
	savelifetimes[0][1] = value-previous[j];
      if(j==12)
	savelifetimes[1][0] = value-previous[j];
      if(j==11)
	savelifetimes[1][1] = value-previous[j];
      previous[j] = value;

    }
    ratio[0]->SetPoint(i,i,(double)savelifetimes[0][0]/savelifetimes[0][1]);
    ratio[1]->SetPoint(i,i,(double)savelifetimes[1][0]/savelifetimes[1][1]);

    ratio[2]->SetPoint(i,i,sc->GetTS()-previous_ts);
    previous_ts= sc->GetTS();
    // if(i%1000 == 0){
    //   double time_end = get_time();
    //   cout<<setw(5)<<setiosflags(ios::fixed)<<setprecision(1)<<(100.*i)/nentries<<" % done\t"<<(Float_t)i/(time_end - time_start)<<" events/s \r"<<flush;
    // }
  }
  cout << endl;


  cout << "writing histograms" << endl;

  long long rawclock = 0, liveclock = 0, rawtrigger = 0, livetrigger=0;

  for(int i=0;i<scalers;i++){
    hrate[i]->Write("",TObject::kOverwrite);
    hvalue[i]->Write("",TObject::kOverwrite);
    if(strcmp(hrate[i]->GetName(),"rate_OBJ_Scint") == 0){
      cout << "OBJ_Scint\t" << previous[i] << endl;
    } else if(strcmp(hrate[i]->GetName(),"rate_XFP_Scint") == 0){
      cout << "XFP_Scint\t" << previous[i] << endl;
    } else if(strcmp(hrate[i]->GetName(),"rate_Raw_Clock") == 0){
      rawclock = previous[i];
      cout << "Raw_Clock\t" << previous[i] << endl;
    } else if(strcmp(hrate[i]->GetName(),"rate_Live_Clock") == 0){
      liveclock = previous[i];
      cout << "Live_Clock\t" << previous[i] << endl;
    } else if(strcmp(hrate[i]->GetName(),"rate_Raw_Trigger") == 0){
      rawtrigger = previous[i];
      cout << "Raw_Trigger\t" << previous[i] << endl;
    } else if(strcmp(hrate[i]->GetName(),"rate_Live_Trigger") == 0){
      livetrigger = previous[i];
      cout << "Live_Trigger\t" << previous[i] << endl;
    }
  }
  ratio[0]->Write("",TObject::kOverwrite);
  ratio[1]->Write("",TObject::kOverwrite);
  ratio[2]->Write("",TObject::kOverwrite);

  cout << "ratios " << endl;
  if(rawclock>0)
    cout << "clock " <<setprecision(5)<< (double)liveclock/rawclock << endl;
  else
    cout << "clock infinity!!" << endl;
  if(rawtrigger>0)
    cout << "trigger " <<setprecision(5)<< (double)livetrigger/rawtrigger << endl;
  else
    cout << "trigger infinity!!!" << endl;

  outfile->Close();
  delete tr;
  for(int j=0;j<scalers;j++){
    cout << hvalue[j]->GetName() << "\t" << previous[j] << endl;

  }

  double time_end = get_time();
  cout << "Run time " << time_end - time_start << " s." << endl;
  return 0;

}




