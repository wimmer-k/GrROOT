#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <sys/time.h>
#include <signal.h>
#include <cstdlib>

#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TCutG.h"
#include "TKey.h"

#include "CommandLineInterface.hh"
#include "Settings.hh"
#include "RunInfo.hh"
#include "Gretina.hh"
#include "GretinaCalc.hh"
#include "Trace.hh"
#include "Mode3Calc.hh"
#include "S800.hh"
#include "S800Calc.hh"
#include "Calibration.hh"

using namespace TMath;
using namespace std;

bool signal_received = false;
void signalhandler(int sig){
  if (sig == SIGINT){
    signal_received = true;
  }
}

double get_time(){
    struct timeval t;
    gettimeofday(&t, NULL);
    double d = t.tv_sec + (double) t.tv_usec/1000000;
    return d;
}
int main(int argc, char* argv[]){
  signal(SIGINT,signalhandler);

  vector<char*> InputFiles;
  char* OutputLoc = NULL;
  char* CutFile = NULL;
  vector<char*> SettingsFile;
  int LastEvent = -1;
  int tac = 1;
  //bool rawfiles = false;
  CommandLineInterface* interface = new CommandLineInterface();
  interface->Add("-i", "input file(s)", &InputFiles);
  interface->Add("-o", "output file", &OutputLoc);
  interface->Add("-c", "cutfile", &CutFile);
  interface->Add("-t", "0 for OBJ TAC, 1 for OBJ TDC, 2 for XFP TAC, 3 for XFP TDC (default = 1)", &tac);
  interface->Add("-n", "last event", &LastEvent);
  //interface->Add("-rf","split raw files instead of cal files", &rawfiles);
  interface->Add("-s","settings file (used only for splitting raw files)",&SettingsFile);
  interface->CheckFlags(argc, argv);

  //Check inputs
  if(InputFiles.size()==0){
    cout << "No input file given " << endl;
    return 1;
  }
  if(OutputLoc == NULL){
    cout << "No output ROOT file given " << endl;
    return 1;
  }
  if(CutFile == NULL){
    cout << "No cut file given " << endl;
    return 1;
  }
  if(tac<0 || tac>3){
    cout << "Unknown tac specification" << endl;
    return 1;
  }
  string output_s = OutputLoc;
  bool splitfiles = !(output_s.length() > 5 &&
		      output_s.substr(output_s.length()-5) == ".root");

  bool rawfiles;
  {
    TFile* tf = new TFile(InputFiles[0]);
    rawfiles = (tf->Get("gtr")!=NULL);
    delete tf;
  }


  if(rawfiles && SettingsFile.size()==0){
    cout << "Must supply a settings file to split raw files" << endl;
    return 1;
  }

  //Tell the user the settings
  cout << "input file(s): ";
  for(vector<char*>::iterator it = InputFiles.begin(); it!=InputFiles.end(); it++){
    cout << *it << " ";
  }
  cout << endl;
  if(splitfiles){
    cout << "output dir: " << OutputLoc << endl;
  } else {
    cout << "output file: " << OutputLoc << endl;
  }
  cout << "cuts file: " << CutFile << endl;
  cout << "Cuts were applied on ";
  if(tac==0){
    cout << "OBJ TAC";
  } else if (tac==1){
    cout << "OBJ TDC";
  } else if (tac==2){
    cout << "XFP TAC";
  } else if (tac==1){
    cout << "XFP TDC";
  }
  cout << " parameters" << endl;
  cout << "------------------------------------" << endl;

  S800* s800 = new S800;
  S800Calc* s800Calc = new S800Calc;
  Gretina* gretina = new Gretina;
  GretinaCalc* gretinaCalc = new GretinaCalc;
  Mode3Event* mode3 = new Mode3Event;
  Mode3Calc* mode3Calc = new Mode3Calc;
  Settings* set = NULL;
  Calibration* cal = NULL;
  //RunInfo* info = NULL;

  //Set up reading of input file(s)
  TChain* input;
  if(rawfiles){
    input = new TChain("gtr");
    input->SetBranchAddress("mode3Event",&mode3);
    input->SetBranchAddress("s800",&s800);
    input->SetBranchAddress("gretina",&gretina);
    set = new Settings(SettingsFile);
    cal = new Calibration(set);
  } else {
    input = new TChain("ctr");
    input->SetBranchAddress("gretinacalc",&gretinaCalc);
    input->SetBranchAddress("s800calc",&s800Calc);
    input->SetBranchAddress("mode3calc",&mode3Calc);
    //Look through each root file to find the settings file.
    for(vector<char*>::iterator it = InputFiles.begin();
	it!=InputFiles.end() && set==NULL; it++){
      TFile* tf = new TFile(*it);
      set = (Settings*)tf->Get("settings");
      delete tf;
    }
  }
  for(vector<char*>::iterator it = InputFiles.begin(); it!=InputFiles.end(); it++){
    input->Add(*it);
  }

  //Read in the cuts file for incoming and outgoing particle ID
  vector<TCutG*> InPartCut;
  vector<vector<TCutG*> > OutPartCut;
  TFile *cFile = new TFile(CutFile);
  TKey *key;

  //Loop over, finding all incoming cuts.  (Cuts with "in", but not "out")
  TIter iterIn(cFile->GetListOfKeys());
  while((key=(TKey*)iterIn())){
    string name = key->GetName();
    string classname = key->GetClassName();
    if(classname=="TCutG" &&
       name.find("in")!=string::npos &&
       name.find("out")==string::npos){
      cout << "incoming cut found " << name << endl;
      InPartCut.push_back((TCutG*)cFile->Get(name.c_str()));
    }
  }
  OutPartCut.resize(InPartCut.size());

  //Loop again, finding all outgoing cuts.  (Cuts with "in" and "out")
  //Match the outgoing cut to an incoming cut (outgoing cut must start with name of incoming cut)
  TIter iterOut(cFile->GetListOfKeys());
  while((key=(TKey*)iterOut())){
    string name = key->GetName();
    string classname = key->GetClassName();
    if(classname=="TCutG" &&
       name.find("in")!=string::npos &&
       name.find("out")!=string::npos){
      for(unsigned short i=0;i<InPartCut.size();i++){
	string inName = InPartCut[i]->GetName();
	if(name.find(inName)==0){
	  OutPartCut[i].push_back((TCutG*)cFile->Get(name.c_str()));
	  cout << "outgoing cut found "<< name << endl;
	}
      }
    }
  }
  cFile->Close();

  //Set up output file for single file output.
  TFile* ofile = NULL;
  if(splitfiles){
    string command = "mkdir -p " + output_s;
    system(command.c_str());
  } else {
    ofile = new TFile(OutputLoc,"RECREATE");
    set->Write("settings",TObject::kOverwrite);
  }

  //Make a tree for each outgoing cut.
  //If splitting into multiple files, make each output file
  vector<vector<TTree*> > splittree;
  splittree.resize(InPartCut.size());
  vector<vector<TFile*> > output_files;
  output_files.resize(InPartCut.size());

  for(UShort_t in=0;in<InPartCut.size();in++){
    for(UShort_t ou=0;ou<OutPartCut[in].size();ou++){
      string cutname = OutPartCut[in][ou]->GetName();
      string treename;
      if(splitfiles){
	string filename = Form("%s/%s.root",OutputLoc,cutname.c_str());
	output_files[in].push_back(new TFile(filename.c_str(),"RECREATE"));
	set->Write("settings",TObject::kOverwrite);
	treename = rawfiles ? "gtr" : "ctr";
      } else {
	treename = rawfiles ? Form("gtr_%s",cutname.c_str()) : Form("ctr_%s",cutname.c_str());
      }
      TTree* newTree = new TTree(treename.c_str(),"Gretina/S800");
      if(rawfiles){
	newTree->Branch("s800",&s800, 320000);
	newTree->Branch("gretina",&gretina, 320000);
	newTree->Branch("mode3Event",&mode3, 320000);
      } else {
	newTree->Branch("s800calc",&s800Calc, 320000);
	newTree->Branch("gretinacalc",&gretinaCalc, 320000);
	newTree->Branch("mode3calc",&mode3Calc, 320000);
      }
      newTree->BranchRef();
      splittree[in].push_back(newTree);
    }
  }

  Double_t nentries = input->GetEntries();
  cout << "Nr of Events " << nentries << endl;
  if (LastEvent!=-1){
    nentries = LastEvent;
    cout << "reading only first " << nentries << " events " << endl;
  }

  double time_start = get_time();
  long events_written = 0;
  for(int i=0; i<nentries && !signal_received; i++){
    s800->Clear();
    s800Calc->Clear();
    gretina->Clear();
    gretinaCalc->Clear();
    mode3->Clear();
    mode3Calc->Clear();

    int status = input->GetEvent(i);

    if(status == -1){
      cerr << "IO error reading event" << endl;
      return 5;
    } else if(status == 0){
      cerr << "Error occured, entry " << i << " doesn't exist" << endl;
      return 6;
    }

    if(rawfiles){
      Gretina* grClone = (Gretina*)gretina->Clone();
      S800* sClone = (S800*)s800->Clone();
      Mode3Event* m3Clone = (Mode3Event*)mode3->Clone();
      cal->BuildAllCalc(sClone,grClone,m3Clone,
			s800Calc,gretinaCalc,mode3Calc);
      delete grClone;
      delete sClone;
      delete m3Clone;
    }

    TOF* tof = s800Calc->GetTOF();
    for(UShort_t in=0;in<InPartCut.size();in++){

      if( (tac%2==1 && InPartCut[in]->IsInside(tof->GetOBJ(),tof->GetXFP())) ||
	  (tac%2==0 && InPartCut[in]->IsInside(tof->GetTACOBJ(),tof->GetTACXFP())) ){
	double icDE = s800Calc->GetIC()->GetDE();
	for(UShort_t ou=0;ou<OutPartCut[in].size();ou++){
	  if((tac == 0 && OutPartCut[in][ou]->IsInside(tof->GetTACOBJC(),icDE)) ||
	     (tac == 1 && OutPartCut[in][ou]->IsInside(tof->GetOBJC(),icDE)) ||
	     (tac == 2 && OutPartCut[in][ou]->IsInside(tof->GetTACXFPC(),icDE)) ||
	     (tac == 3 && OutPartCut[in][ou]->IsInside(tof->GetXFPC(),icDE))
	     ){

	      splittree[in][ou]->Fill();
	      events_written++;

	  }
	}
      }
    }

    //Progress meter
    if(i%1000 == 0){
      double time_end = get_time();
      cout << setw(5) << setiosflags(ios::fixed) << setprecision(1)
	   << (100.*i)/nentries<<" % done\t"
	   << (Float_t)i/(time_end - time_start) << " events/s "
	   << (nentries-i)*(time_end - time_start)/(Float_t)i << "s to go \r" << flush;
    }
  }
  cout <<  endl;


  //Close all trees.
  cout << "closing file ..." << endl;
  for(UShort_t in=0;in<InPartCut.size();in++){
    for(UShort_t ou=0;ou<OutPartCut[in].size();ou++){
      if(splitfiles){
	output_files[in][ou]->cd();
      }
      splittree[in][ou]->Write("",TObject::kOverwrite);
      if(splitfiles){
	output_files[in][ou]->Close();
      }
    }
  }

  if(!splitfiles){
    ofile->Close();
  }
  double time_end = get_time();
  cout << "Program Run time " << time_end - time_start << " s." << endl;
  cout << "Calculated " << nentries/(time_end - time_start) << " events/s." << endl;
  cout << events_written << "/" << nentries << " ("
       << 100*double(events_written)/double(nentries) << "%) written" << endl;
  return 0;
}
