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
#include <string>
#include <vector>
#include <sys/time.h>
#include <signal.h>

#include "TFile.h"
#include "TTree.h"
#include "TStopwatch.h"
#include "Gretina.hh"
#include "Gretinadefs.h"

#include "CommandLineInterface.hh"
#include "UnpackedEvent.hh"
#include "Settings.hh"
#include "Calibration.hh"
#include "RunInfo.hh"

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
  double time_start = get_time();
  TStopwatch timer;
  timer.Start();
  signal(SIGINT,signalhandler);

  int LastBuffer =-1;
  char *InputFile = NULL;
  char *RootFile = NULL;
  vector<char*> SettingFile;
  int denom = 10000;
  bool wrawtree = false;
  bool wrawhist = false;
  bool wcaltree = false;
  bool wcalhist = false;
  bool noHFC = false;

  //Read in the command line arguments
  CommandLineInterface* interface = new CommandLineInterface();
  interface->Add("-lb", "last buffer to be read", &LastBuffer);
  interface->Add("-i", "input file", &InputFile);
  interface->Add("-o", "output file", &RootFile);
  interface->Add("-s", "settingsfile", &SettingFile);
  interface->Add("-d", "output every x buffers/event/block/headers...", &denom);
  interface->Add("-rt", "write raw tree", &wrawtree);
  interface->Add("-rh", "write raw histos", &wrawhist);
  interface->Add("-ct", "write cal tree", &wcaltree);
  interface->Add("-ch", "write cal histos", &wcalhist);
  interface->Add("--no-HFC","do not use HFC as an intermediate step",&noHFC);
  interface->CheckFlags(argc, argv);


  //Tell the user what will be written.
  if(wrawhist)
    cout<<"writing raw histograms " << endl;
  if(wrawtree)
    cout<<"writing raw tree " << endl;
  if(wcalhist)
    cout<<"writing cal histograms " << endl;
  if(wcaltree)
    cout<<"writing cal tree " << endl;

  if(!wrawtree && !wrawhist && !wcaltree && !wcalhist){
    cout << "No output format specified." << endl;
    cout << "Assuming that the raw tree is to be written." << endl;
    wrawtree = true;
  }

  //Complain about missing mandatory arguments
  if(InputFile == NULL){
    cout << "No input file given " << endl;
    return 1;
  }
  if(RootFile == NULL){
    cout << "No output ROOT file given " << endl;
    return 2;
  }
  if(SettingFile.size() == 0){
    cout << "No settings file given " << endl;
    return 3;
  }

  //get the run number from the filename
  int run;
  TString ifname(InputFile);
  ifname.Remove(0,ifname.Length()-21); // Last 21 characters: RunXXXX/Global.dat.gz but works also for Global.dat
  sscanf(ifname.Data(),"%*sRun%04d/%*s",&run);

  //Open the input and output files.
  TFile *ofile = new TFile(RootFile,"RECREATE");
  cout<<"input file: "<<InputFile<< endl;
  cout<<"run number " << run << endl;
  cout<<"writing to output file: "<<RootFile<< endl;
  cout << "------------------------------------" << endl;
  FILE *infile;

  string infilestr = string(InputFile);
  char extension[20];
  strcpy(extension, infilestr.substr(infilestr.find_last_of(".")+1).c_str());
  if((strcmp(extension,"gz")==0 || strcmp(extension,"gzip")==0)
     && !noHFC){
    infile = popen(Form("zcat %s | GEB_HFC -p ", InputFile),"r");
  } else if (strcmp(extension,"dat")==0 && !noHFC){
    infile = popen(Form("GEB_HFC -p %s ", InputFile),"r");
  } else if ((strcmp(extension,"gz")==0 || strcmp(extension,"gzip")==0)
	     && noHFC){
    infile = popen(Form("zcat %s ", InputFile),"r");
  } else if (strcmp(extension,"dat")==0 && noHFC){
    infile = fopen(InputFile,"r");
  } else {
    cout << "Unknown file type.  Will try to open as a .dat file, using HFC" << endl;
    infile = popen(Form("GEB_HFC -p %s ", InputFile),"r");
  }

  if(infile == NULL){
    cout << "Sorry I couldn't find the file: " << InputFile << ". Aborting ..." << endl;
    ofile->Close();
    return 3;
  }

  //Read the settings
  Settings* set = new Settings(SettingFile);
  ofile->cd();
  int vl = set->VLevel();

  //store important information
  RunInfo* info = new RunInfo();
  info->SetRunNr(run);

  //Initialize the data structures for the event building.
  int buffers = 0, sbuffers = 0;
  long long int bytes_read = 0;
  struct crys_ips_abcd1234 inbuf_abcd1234[1];
  struct crys_ips_abcd5678 inbuf_abcd5678[1];
  UnpackedEvent *evt = new UnpackedEvent(set);
  evt->SetWrite(wrawtree, wrawhist, wcaltree, wcalhist);
  evt->SetVL(vl);
  Calibration *cal = new Calibration(set);
  evt->SetCalibration(cal);
  evt->Init();

  //Create the scaler tree in the output file.
  TTree *sc = new TTree("sc", "Scaler Tree");
  Scaler *ScalerBranch = evt->GetScaler();
  sc->Branch("scaler",&ScalerBranch,32000);
  sc->BranchRef();
  if(vl>0)
    cout << "start " << endl;

  //Loop over the entirety of the input file.
  while(!feof(infile) && !signal_received){

    //Finish reading if you have read as many buffers as the user has requested.
    if(LastBuffer > 0 && buffers >= LastBuffer)
      break;

    size_t bsize;
    unsigned short buffer[4096], evtlength[1], *pevent;

    //Read the header of the new buffer.
    unsigned int header[4];
    bsize = fread(header, sizeof(unsigned int), 4, infile);
    long long int ts = (long long int)header[3] << 32;
    ts += header[2];



    if(vl>1){
      for(int i=0;i<4;i++)
	cout <<(hex)<< header[i] <<"\t"<<(dec)<< header[i] <<endl;
    }
    bytes_read += 4*sizeof(unsigned int);

    //Decode the event differently depending on which type of data is identified in the header.
    //For each, pass the data from the file into the UnpackedEvent to be read.
    if(header[0]==GRETINA_ID){
      if(vl>0)
	cout << "gret timestamp:\t"<< ts << "\tlength: "<< header[1] << "\thex: " <<(hex) << ts << "\tlength: "<< header[1] <<(dec)<< endl;
      Crystal* crys;
      if(header[1]==sizeof(crys_ips_abcd1234)){
	bsize = fread(&inbuf_abcd1234[0], sizeof(crys_ips_abcd1234), 1, infile);
	bytes_read += sizeof(struct crys_ips_abcd1234);
	crys = new Crystal(inbuf_abcd1234[0]);
      } else if (header[1]==sizeof(crys_ips_abcd5678)){
	bsize = fread(&inbuf_abcd5678[0], sizeof(crys_ips_abcd5678), 1, infile);
	bytes_read += sizeof(struct crys_ips_abcd5678);
	crys = new Crystal(inbuf_abcd5678[0]);
      } else {
	cout << "Unknown size for mode2 data" << endl;
	break;
      }
      int error = evt->DecodeGretina(crys,ts);
      if(error){
	cout << "An error ("<<error<<") occured in DecodeGretina() while processing file: " << InputFile << ". Continuing ..." << endl;
	continue;
      }
    }
    else if(header[0]==S800_ID){
      if(vl>0)
	cout << "s800 timestamp:\t"<< ts << "\tlength: "<< header[1] << "\thex: " <<(hex) << ts << "\tlength: "<< header[1] <<(dec)<< endl;
      bsize = fread(evtlength, sizeof(unsigned short), 1, infile);
      if(bsize == 0){
	cout << "error s800 header but no data" << endl;
	break;
      }
      bytes_read += sizeof(unsigned short);
      evtlength[0]--;
      bsize = fread(buffer, sizeof(unsigned short), evtlength[0], infile);
      bytes_read += evtlength[0]*sizeof(unsigned short);
      if(bsize == 0){
	cout << "error s800 expected " << evtlength[0] << " words but found none" << endl;
	break;
      }
      if(vl>3){
	for(int b=0;b<evtlength[0];b++){
	  cout << b <<"\t" <<hex<< buffer[b] << "\t" <<dec<< buffer[b] << endl;
	}
      }
      pevent = buffer;
      int error = evt->DecodeS800(pevent, ts, evtlength[0]);
      if(error){
	cout << "An error ("<<error<<") occured in DecodeS800() while processing file: " << InputFile << ". Continuing ..." << endl;
	continue;
      }

    }
    else if(header[0]==S800_CARD29_ID){
      char cBuf[16382*4];
      bsize = fread(cBuf, sizeof(char), header[1], infile);
      bytes_read += header[1]*sizeof(char);
      int error = evt->DecodeMode3(cBuf, header[1], ts, true);
      if(error){
	cout << "An error ("<<error<<") occured at buffer nr " << buffers << " in DecodeMode3() while processing file: " << InputFile << ". Continuing ..." << endl;
	//continue;
      }
    }
    else if(header[0]==TRACE_ID){
      if(vl>1){
	cout << "---------------------starting trace------------------------------- " << endl;
	cout << "trace timestamp:\t"<< ts << "\tlength: "<< header[1] << "\thex: " <<(hex) << ts << "\tlength: "<< header[1] <<(dec)<< endl;
      }
      char cBuf[16382*4];
      if(header[1]>4*16382)
	cout << "buffer nr " << buffers << ", size of cBuf not sufficient " << endl;
      bsize = fread(cBuf, sizeof(char), header[1], infile);
      bytes_read += header[1]*sizeof(char);
      int error = evt->DecodeMode3(cBuf, header[1], ts);
      if(error){
	cout << "An error ("<<error<<") occured at buffer nr " << buffers << " in DecodeMode3() while processing file: " << InputFile << ". Continuing ..." << endl;
	//continue;
      }

    }
    else if(header[0]==SCALER_ID || header[0]==NONEVENT_ID){
      int length = (header[1])/2;
      if(vl>1){
	cout << "scaler timestamp:\t"<< ts << "\tlength: "<< header[1] << "\thex: " <<(hex) << ts << "\tlength: "<< header[1] <<(dec)<< endl;
	if(wrawtree){
	  cout << "current event " << evt->GetTree()->GetEntries() << endl;
	}
      }
      bsize = fread(buffer, sizeof(unsigned short), length, infile);
      bytes_read += length*sizeof(unsigned short);

      if(bsize == 0){
	cout << "error scaler header but no data" << endl;
	cout << "scaler timestamp:\t"<< ts << "\tlength: "<< header[1] << "\thex: " <<(hex) << ts << "\tlength: "<< header[1] <<(dec)<< endl;
	continue;
      }
      if(vl>3){
	for(int b=0;b<length;b++){
	  cout << b <<"\t" <<hex<< buffer[b] << "\t" <<dec<< buffer[b] << endl;
	}
      }
      pevent = buffer;
      int error = evt->DecodeScaler(pevent, ts);
      if(error == 1){
	//There was a PhysicsEventCountItem, which we are currently ignoring.
	continue;
      } else if (error) {
	cout << "An error ("<<error<<") occured in DecodeS800() while processing file: " << InputFile << ". Continuing ..." << endl;
	continue;
      }

      if(evt->HasScaler()){
	sbuffers++;
 	sc->Fill();
      }
      //break;
    }
    else{
      cout << "unidentified header " << header[0] << "\thex: " <<(hex) << header[0] <<(dec)<< endl;
      cout << "unidentified timestamp:\t"<< ts << "\tlength: "<< header[1] << "\thex: " <<(hex) << ts << "\tlength: "<< header[1] <<(dec)<< endl;
      cout << "at buffer nr " << buffers << endl;
      if(buffers==0)
	continue;
      break;
    }

    //Write the trees out to disk every denom events.
    buffers++;
    if(buffers % denom == 0){
      if(wrawtree)
	evt->GetTree()->AutoSave();
      if(wcaltree)
	evt->GetCalTree()->AutoSave();
      sc->AutoSave();
      double time_end = get_time();
      cout << "\r" << buffers << " buffers read... "<<bytes_read/(1024*1024)<<" MB... "<<buffers/(time_end - time_start) << " buffers/s" << flush;
    }
  }

  //Finish reading the last event and close it out.
  evt->WriteLastEvent();
  for(int j=0;j<NSCALER;j++)
    info->Integral(j,evt->GetScaler()->GetValue(j));
  info->SetRunTime(evt->GetRunTime());
  info->SetBuffer(buffers,sbuffers);
  info->SetEvents(evt->NrOfEvents());

  cout << "Total of " << buffers << " data buffers (" << "run time " << evt->GetRunTime()<<" s and "<<bytes_read/(1024*1024)<<" MB) and" << endl;
  cout << "         " << sbuffers << " scaler buffers read." << endl;

  if(wrawtree||wrawhist){
    cout << "Total of " << evt->NrOfEvents() << " raw events";
    if(wrawtree)
      cout <<"("<<evt->GetTree()->GetZipBytes()/(1024*1024)<<" MB)";
    cout <<" written."  << endl;
    cout << evt->NrOfHits() << " hits and " << evt->NrOfStrangeHits() << " strange hits (bad ip) "<<setprecision(2)<< (float)evt->NrOfStrangeHits()/evt->NrOfHits()*100.<< " %"  <<setw(5)<< endl;
  }
  if(wcaltree){
    info->SetEntries(evt->NrOfCalEvents());
    cout << "Total of " << evt->NrOfCalEvents() << " cal events ("<<evt->GetCalTree()->GetZipBytes()/(1024*1024)<<" MB) written."  << endl;
  }
  double time_end = get_time();
  cout << "Program Run time " << time_end - time_start << " s." << endl;
  cout << "Unpacked " << buffers/(time_end - time_start) << " buffers/s." << endl;
  timer.Stop();
  cout << "\n CPU time: " << timer.CpuTime() << "\tReal time: " << timer.RealTime() << endl;

  ofile->cd();
  //Final cleanup and writing of files.
  if(wrawtree){
    evt->GetTree()->Write("",TObject::kOverwrite);
  }

  if(wcaltree){
    evt->GetCalTree()->Write("",TObject::kOverwrite);
  }

  sc->Write("",TObject::kOverwrite);
  info->Write("runinfo",TObject::kOverwrite);
  if(wcaltree){
    set->Write("settings",TObject::kOverwrite);
  }
  //Cannot simply close ofile, even though it was the one opened at the beginning.
  //If the tree spills out into a new file after 1.8 GB, it closes the original file and opens a new one.
  //ofile then has a dead pointer.
  //However, sc->GetCurrentFile() still points to the currently open file.
  sc->GetCurrentFile()->Close();
  //ofile->Close();
  return 0;
}
