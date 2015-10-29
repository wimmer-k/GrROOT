#include "Simulation.hh"
using namespace std;


Simulation::Simulation(SimSettings* settings): UnpackedEvent(settings){
  cout <<"contstructing simulation "
       << endl;
  fSett = settings;
}
void Simulation::Init(){
  UnpackedEvent::Init();
  fRand = new TRandom();
  fGretinaSim = new GretinaSim;
  if(fwsimtree){
    cout << "setting up simulation tree " << endl;
    fsimtr = new TTree("str","Geant4 emitted gamma rays");
    fsimtr->Branch("GretinaSim",&fGretinaSim, 320000);
    fsimtr->BranchRef();
    cout << "Simulation: " << "done setting up simulation tree" << endl;
  }
  ReadSimResolution(fSett->SimResolutionFile());
  ReadSimThresholds(fSett->SimThresholdFile());
  fcal->SetHasMap(true);
  fGretinaSim->Clear();
}
int Simulation::DecodeGretinaG4Sim(G4SIM_EGS* g4Sim, long long int ts){

  fGretinaSim->SetTimeStamp(ts);
  for(Int_t i = 0; i < g4Sim->num; i++){
    fGretinaSim->AddEmittedGamma(&g4Sim->gammas[i]);
  }

  if(fvl>2){
    cout << "Simultaion: GretinaSim: " << endl;
    for(Int_t i = 0; i < fGretinaSim->GetMult(); i++)
      cout << "e = " << fGretinaSim->GetEmittedGamma(i)->GetEnergy()
	   << " (x, y, z) = (" 
	   << fGretinaSim->GetEmittedGamma(i)->GetPos().X() << ", "
	   << fGretinaSim->GetEmittedGamma(i)->GetPos().Y() << ", "
	   << fGretinaSim->GetEmittedGamma(i)->GetPos().Z() 
	   << ")  (phi, theta) = ("
	   << fGretinaSim->GetEmittedGamma(i)->GetPhi() << ", "
	   << fGretinaSim->GetEmittedGamma(i)->GetTheta() << ")" << endl;

    cout << "Simulation: GretinaSim event found with timestamp " << ts << endl;
  }

  // These events go into a separate tree as singles. No need to pay 
  // attention to the time window.
  
  fsimtr->Fill();
  fnsimentries++;
  fGretinaSim->Clear();

  return 0;
}
int Simulation::DecodeS800PhysicsData(S800_PHYSICSDATA* s800pd, long long int ts){

  if(fvl>2)
    cout << "Simulation: S800_PHYSICSDATA: ata = " << s800pd->ata 
	 << " yta = "  << s800pd->yta 
	 << " bta = "  << s800pd->bta 
	 << " dta = "  << s800pd->dta
	 << " ts = "  << ts << endl;

  fS800Calc->SetTimeStamp(ts);
  fS800Calc->GetTRACK()->SetATA(s800pd->ata);
  fS800Calc->GetTRACK()->SetYTA(s800pd->yta);
  fS800Calc->GetTRACK()->SetBTA(s800pd->bta);
  fS800Calc->GetTRACK()->SetDTA(s800pd->dta*100);//KW: GrROOT takes dta in %

  if(fvl>2){
    cout << "UnpackedEvent: S800PhysicsData event found with timestamp " << ts << endl;
  }
  return 0;
}
void Simulation::ReadSimResolution(const char* filename){
  TEnv* env = new TEnv(filename);
   if(fvl>1){
     env->Print();
   }
  double globA = env->GetValue("Detector.All.Crystal.All.A",0.0);
  double globB = env->GetValue("Detector.All.Crystal.All.B",0.0);
  double globC = env->GetValue("Detector.All.Crystal.All.C",0.0);
  for(int det=0;det<7;det++){
    for(int cr=0;cr<4;cr++){
      int toRes = env->GetValue(Form("Detector.%d.Crystal.%d.SimResolution",det,cr),0);
      if(toRes){
	simresolution newSimres;
	newSimres.Detector = det;
	newSimres.Crystal = cr;
	newSimres.A = env->GetValue(Form("Detector.%d.Crystal.%d.A",det,cr),0.0);
	newSimres.B = env->GetValue(Form("Detector.%d.Crystal.%d.B",det,cr),0.0);
	newSimres.C = env->GetValue(Form("Detector.%d.Crystal.%d.C",det,cr),0.0);
	fSimResolutions.push_back(newSimres);
      } else {
	simresolution newSimres;
	newSimres.Detector = det;
	newSimres.Crystal = cr;
	newSimres.A = globA;
	newSimres.B = globB;
	newSimres.C = globC;
	fSimResolutions.push_back(newSimres);
      }
    }
  }
}

void Simulation::ReadSimThresholds(const char* filename){
  TEnv* env = new TEnv(filename);
  for(int det=0;det<7;det++){
    for(int cr=0;cr<4;cr++){
      int toThresh = env->GetValue(Form("Detector.%d.Crystal.%d.SimThreshold",det,cr),0);
      if(toThresh){
	simthreshold newSimthresh;
	newSimthresh.Detector = det;
	newSimthresh.Crystal = cr;
	newSimthresh.E = env->GetValue(Form("Detector.%d.Crystal.%d.E",det,cr),0.0);
	newSimthresh.dE = env->GetValue(Form("Detector.%d.Crystal.%d.dE",det,cr),0.001);
	fSimThresholds.push_back(newSimthresh);
      }
    }
  }
}
bool Simulation::SimResolution(Gretina* gr){
  //cout << __PRETTY_FUNCTION__ << " with resolution " << fSett->SimGretinaPositionResolution() << endl;
  for(int i=0; i<gr->GetMult(); i++){
    Crystal* crys = gr->GetHit(i);
    //Position resolution
    if(fSett->SimGretinaPositionResolution()>0){
      for(int j=0; j<crys->GetMult(); j++){
	crys->GetIPoint(j)->SetPosition(fRand->Gaus(crys->GetIPoint(j)->GetPosition().X(),
						    fSett->SimGretinaPositionResolution()),
					fRand->Gaus(crys->GetIPoint(j)->GetPosition().Y(),
						    fSett->SimGretinaPositionResolution()),
					fRand->Gaus(crys->GetIPoint(j)->GetPosition().Z(),
						    fSett->SimGretinaPositionResolution()) );
      }
    }
    //Energy resolutions
    int det = fSett->Hole2Det(crys->GetHole());
    int crysnum = crys->GetCrystal();
    for(vector<simresolution>::iterator it = fSimResolutions.begin(); it!=fSimResolutions.end(); it++){
      if(det==it->Detector && crysnum==it->Crystal){
	crys->SetEnergy(fRand->Gaus(crys->GetEnergy(),
				    it->A*sqrt(1.0+crys->GetEnergy()*it->B) + it->C*crys->GetEnergy()));
	for(int j=0; j<crys->GetMult(); j++){
	  IPoint* ipoint = crys->GetIPoint(j);
	  ipoint->SetEnergy(fRand->Gaus(ipoint->GetEnergy(),
					it->A*sqrt(1.0+ipoint->GetEnergy()*it->B) + it->C*ipoint->GetEnergy()));
	}
	break;
      }
    }
  }

  return true;
}

bool Simulation::SimThresholds(Gretina* gr){
  //cout << __PRETTY_FUNCTION__ << endl;
  for(int i=0; i<gr->GetMult(); i++){
    Crystal* crys = gr->GetHit(i);
    int det = fSett->Hole2Det(crys->GetHole());
    int crysnum = crys->GetCrystal();
    for(vector<simthreshold>::iterator it = fSimThresholds.begin(); it!=fSimThresholds.end(); it++){
      if(det==it->Detector && crysnum==it->Crystal){
	if( fRand->Uniform(0,1) >
	    0.5*(1.0 + tanh( (crys->GetEnergy() - it->E) / it->dE ) ) ){
	  crys->SetEnergy(0.);
	}
	break;
      }
    }
  }

  return true;
}
void Simulation::CloseEvent(){
  //cout << __PRETTY_FUNCTION__ << endl;
  SimResolution(fGretina);
  SimThresholds(fGretina);
  UnpackedEvent::CloseEvent();
}
