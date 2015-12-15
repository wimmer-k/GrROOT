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
#include "Calibration.hh"

#include "S800.hh"
#include "Gretina.hh"
#include "GretinaCalc.hh"
#include "Gretinadefs.h"
#include "Trace.hh"

using namespace std;

#define uint unsigned int

struct recal{
  int Detector;
  int Crystal;
  double Slope;
  double Offset;
};

Calibration::Calibration(){
  ResetCtrs();
  ftracking = new Tracking;
}
Calibration::Calibration(Settings* setting, int event){
  ResetCtrs();
  fhasmap = false;
  fSett = setting;
  fRand = new TRandom();
  fped.resize(2);
  fslope.resize(2);
  foffset.resize(2);
  fbad.resize(2);
  for(int i=0;i<2;i++){
    fped[i].resize(fSett->CrdcChannels());
    fslope[i].resize(fSett->CrdcChannels());
    foffset[i].resize(fSett->CrdcChannels());
  }
  fcrdccal.resize(fSett->CrdcChannels());
  ReadCalibration(fSett->CalFile(),fSett->PedestalFile());
  ReadGretinaRecalibration(fSett->GretinaRecalFile());
  ReadBadPads(fSett->BadFile());
  fICoffset.resize(fSett->IonChamberChannels());
  fICslope.resize(fSett->IonChamberChannels());


  fobj_cor = NULL;
  fxfp_cor = NULL;
  fobjtac_cor = NULL;
  fxfptac_cor = NULL;
  fy_cor[0] = NULL;
  fy_cor[1] = NULL;
  fIC_cor[0] = NULL;
  fIC_cor[1] = NULL;

  ReadTimeCorrections(fSett->TimeCorFile());
  fevent =event;

  ReadICCalibration(fSett->CalFileIC());
  ReadHodoCalibration(fSett->HodoCalFile());

  ReadNeighbors(fSett->NeighborFile());
  ReadGretinaCal(fSett->GretinaCalFile());

  ReadMatrix(fSett->MatrixFile());


  fverbose = fSett->VLevel();
  if(fverbose>=1)
    PrintCalibration();
  ReadMap(fSett->MapFile());

  fAddBackType = fSett->AddBackType();
  fbeta = fSett->TargetBeta();
  fmode3 = new Mode3Calc();

  ftracking = new Tracking((TrackSettings*)setting);
}

Calibration::~Calibration(){}

void Calibration::ReadICCalibration(const char *filename){
  TEnv *iccal = new TEnv(filename);
  for(int i=0;i<fSett->IonChamberChannels();i++){
    fICoffset[i] = iccal->GetValue(Form("IonChamber.Offset.%d",i),0.0);
    fICslope[i] = iccal->GetValue(Form("IonChamber.Slope.%d",i),1.0);
  }
  fde_slope = iccal->GetValue("IonChamber.Slope.DE",1.0);
  fde_offset = iccal->GetValue("IonChamber.Offset.DE",0.0);
}
void Calibration::ReadHodoCalibration(const char *filename){
  TEnv *hodocal = new TEnv(filename);
  for(int i=0;i<32;i++){
    fhodooffset[i] = hodocal->GetValue(Form("Hodoscope.Offset.%d",i),0.0);
    fhodogain[i] = hodocal->GetValue(Form("Hodoscope.Gain.%d",i),1.0);
  }
}
void Calibration::ReadBadPads(const char *filename){
  TEnv *bad = new TEnv(filename);
  for(UShort_t i=0;i<2;i++){
    fbad[i].resize(bad->GetValue(Form("Crdc.%d.Nofbadpads",i),0));
    for(UShort_t k=0;k<fbad[i].size();k++){
      fbad[i][k] = bad->GetValue(Form("Crdc.%d.badpad.%d",i,k),0);
    }
  }
}
void Calibration::ReadCalibration(const char *filename, const char *pedfile){
  TEnv *crdcpedestals = new TEnv(pedfile);
  TEnv *crdccal = new TEnv(filename);
  for(int c=0;c<2;c++){//crdc
    for(int p=0;p<fSett->CrdcChannels();p++){//pads
      fped[c][p] = crdcpedestals->GetValue(Form("Crdc.%d.Ped.%03d",c,p),0.0);
      fslope[c][p] = crdccal->GetValue(Form("Crdc.%d.Slope.%03d",c,p),1.0);
      foffset[c][p] = crdccal->GetValue(Form("Crdc.%d.Offset.%03d",c,p),0.0);
    }
  }
}

void Calibration::ReadGretinaRecalibration(const char* filename){
  TEnv* env = new TEnv(filename);
  for(int det=0;det<7;det++){
    for(int cr=0;cr<4;cr++){
      int toRecal = env->GetValue(Form("Detector.%d.Crystal.%d.Recal",det,cr),0);
      //The file says to recalibrate this hit, so record the values for later use.
      if(toRecal){
	recal newRecal;
	newRecal.Detector = det;
	newRecal.Crystal = cr;
	newRecal.Slope = env->GetValue(Form("Detector.%d.Crystal.%d.Slope",det,cr),1.0);
	newRecal.Offset = env->GetValue(Form("Detector.%d.Crystal.%d.Offset",det,cr),0.0);
	fGretRecal.push_back(newRecal);
      }
    }
  }
}

void Calibration::ReadMatrix(const char *filename){

  if(fSett->VLevel()>0)
    cout << filename << endl;
  int in, sizr;
  char str[132];
  sprintf(str, "%s", filename);
  in = open(str, O_RDONLY, 0);
  if(fSett->VLevel()>0){
    if(in != 0){
      printf("%s is open (input) binary format\n", str);
    }
    else{
      printf("Could not open %s\n", str);
    }
  }
  sizr = read(in, (char *) fcrmat, sizeof(fcrmat));
  if(fSett->VLevel()>0)
    printf("Read %i bytes into fcrmat\n", sizr);
  close(in);
  if(fSett->VLevel()>1){
    for(int i=0;i<MAXDETPOS;i++){
      for(int j=0;j<MAXCRYSTALNO;j++){
	cout << "Hole : " << i << "\tCrystal: " << j << endl;
	for(int k=0;k<4;k++){
	  for(int m=0;m<4;m++){
	    cout << fcrmat[i][j][k][m] << "\t";
	  }
	  cout << endl;
	}
	cout <<"-------------------------------------"<< endl;
      }
      cout <<"-------------------------------------"<< endl;
    }
  }
}

void Calibration::PrintCalibration(){
  for(UShort_t i=0;i<2;i++){
    cout << "CRDC " << i << endl;
    cout << fbad[i].size() << " bad pads " << endl;
    for(UShort_t chan=0; chan<fSett->CrdcChannels(); chan++){
      for(UShort_t b=0;b<fbad[i].size();b++){
	if(chan==fbad[i][b])
	  cout << "BAD " << endl;
      }
      cout << "Ped " << chan<< ": " << fped[i][chan] << endl;
      cout << "Slope " << chan<< ": " << fslope[i][chan] << endl;
      cout << "Offset " << chan<< ": " << foffset[i][chan] << endl;
    }
  }
  cout << "Ion Chamber" << endl;
  for(int i=0;i<fSett->IonChamberChannels();i++){
    cout << "Ped " << i << ": " << fICoffset[i] << endl;
    cout << "Slope " << i << ": " << fICslope[i] << endl;
  }
  cout << "DE Ped: " << fde_offset << endl;
  cout << "DE Slope: " << fde_slope << endl;

  cout << "Hodoscope" << endl;
  for(int i=0;i<32;i++){
    cout << "Ped " << i << ": " << fhodooffset[i] << endl;
    cout << "Gain " << i << ": " << fhodogain[i] << endl;
  }
}
void Calibration::ReadMap(const char* filename){
  fmaxcoefficient.resize(fSett->TrackParams());
  forder.resize(fSett->TrackParams());
  fexponent.resize(fSett->TrackParams());
  fcoefficient.resize(fSett->TrackParams());
  for(UShort_t i=0;i<fSett->TrackParams();i++){
    forder[i].resize(fSett->TrackCoefs());
    fcoefficient[i].resize(fSett->TrackCoefs());
    fexponent[i].resize(fSett->TrackParams());
    for(UShort_t k=0;k<fSett->TrackParams();k++){
      fexponent[i][k].resize(fSett->TrackCoefs());
    }
  }
  fmaxorder = 0;
  FILE* file;
  char line[80];
  int index, par, exp[6];
  Int_t ord;
  float co;
  char title[120];
  char *ret;
  file = fopen(filename, "r");
  if(file == NULL){
    cout << "Sorry I couldn't find the map file: " << filename << endl;
    cout << "Will continue without the map file" << endl;
    return;
  }
  ret = fgets(title, 120, file);
  sscanf(title, "S800 inverse map - Brho=%g - M=%d - Q=%d", &fbrho, &fmass, &fcharge);
  if(fSett->VLevel()>0)
    cout << "brho " << fbrho << " mass " << fmass << " charge " << fcharge << endl;
  while(strstr(line, "COEFFICIENT") == NULL)
    ret = fgets(line, 80, file);
  par = 0;
  while(!feof(file)){
    ret = fgets(line, 80, file);
    while (strstr(line, "------------------") == NULL){
      sscanf(line, "%d %g %d %d %d %d %d %d %d", &index, &co, &ord, &exp[0], &exp[1], &exp[2], &exp[3], &exp[4], &exp[5]);
      if(index > fSett->TrackCoefs()){
	cout << "Too many coefficients in map.  Increase TS800_TRACK_COEFFICIENTS." << endl;
	break;
      }
      if(par > fSett->TrackParams()){
	cout << "Too many parameters in map.  Increase TS800_TRACK_PARAMETERS." << endl;
	break;
      }
      fmaxcoefficient[par] = index;
      forder[par][index-1] = ord;
      fcoefficient[par][index-1] = co;
      //cout << "max coef " << fmaxcoefficient[par] << " order " << forder[par][index-1] << " coef " << fcoefficient[par][index-1] << endl;
      for(int k=0; k<fSett->TrackParams(); k++){
	fexponent[par][k][index-1] = exp[k];
	//cout << "exp["<<k<<"] " << exp[k] << " fexponent["<<par<<"]["<<k<<"]["<<index-1<<"] " << fexponent[par][k][index-1] << endl;
      }
      ret = fgets(line, 80, file);
    }
    if(ord > fmaxorder)
      fmaxorder = ord;
    par++;
  }
  //cout << "Done reading map from " << filename << "." << endl;
  //cout << "Title: " << title << endl;
  //cout << "Order: " << fmaxorder << endl;
  fclose(file);
  fhasmap = true;
}
void Calibration::ReadTimeCorrections(const char *filename){
  //Remember the current directory so that we can revert to it later.
  TDirectory* outfile = gDirectory;

  TFile *cor = new TFile(filename);
  if(cor->IsZombie()){
    cerr << "ignore previous warning!!" << endl;
    cerr << "File " << filename << " not existing, cannot perform time dependent corrections!" << endl;
    fobj_cor = NULL;
    fxfp_cor = NULL;
    fobjtac_cor = NULL;
    fxfptac_cor = NULL;
    fy_cor[0] = NULL;
    fy_cor[1] = NULL;
    fIC_cor[0] = NULL;
    fIC_cor[1] = NULL;
    outfile->cd();
  }
  else{
    fobj_cor = (TH1F*)cor->Get("obj_vs_time_cor");
    fxfp_cor = (TH1F*)cor->Get("xfp_vs_time_cor");
    fobjtac_cor = (TH1F*)cor->Get("objtac_vs_time_cor");
    fxfptac_cor = (TH1F*)cor->Get("xfptac_vs_time_cor");
    fy_cor[0] = (TH1F*)cor->Get("y0_vs_time_cor");
    fy_cor[1] = (TH1F*)cor->Get("y1_vs_time_cor");
    fIC_cor[0] = (TH1F*)cor->Get("IC_vs_time_corg");
    fIC_cor[1] = (TH1F*)cor->Get("IC_vs_time_coro");
    outfile->cd();
    fobj_cor->Write("",TObject::kOverwrite);
    fxfp_cor->Write("",TObject::kOverwrite);
    fobjtac_cor->Write("",TObject::kOverwrite);
    fxfptac_cor->Write("",TObject::kOverwrite);
    if(fy_cor[0]!=NULL)
      fy_cor[0]->Write("",TObject::kOverwrite);
    if(fy_cor[1]!=NULL)
      fy_cor[1]->Write("",TObject::kOverwrite);
    if(fIC_cor[0]!=NULL)
      fIC_cor[0]->Write("",TObject::kOverwrite);
    if(fIC_cor[1]!=NULL)
      fIC_cor[1]->Write("",TObject::kOverwrite);
  }

}

Float_t Calibration::MapCalc(int calcorder, int parameter, Float_t *input){
  Float_t cumul=0;
  Float_t multiplicator;
  for(int index=0; index<fmaxcoefficient[parameter]; index++){
    if (calcorder < forder[parameter][index]) break;
    multiplicator = 1;
    for(int nex=0; nex<fSett->TrackParams(); nex++){
      if(fexponent[parameter][nex][index] != 0){
	multiplicator *= pow(input[nex], fexponent[parameter][nex][index]);
      }
    }
    cumul += multiplicator * fcoefficient[parameter][index];
  }
  return cumul;
}


void Calibration::BuildAllCalc(S800* inS800, Gretina* inGret, Mode3Event* inMode3,
			       S800Calc* outS800, GretinaCalc* outGret, Mode3Calc* outMode3){

  //Determine which of the components are present in this event.
  bool hass800 = inS800->GetTS()>0;
  bool hasgret = inGret->GetMult()>0;
  bool hasmode3 = false;
  Trace* card29trace = NULL;
  //inMode3 will also have card29 data, but we shouldn't bother calibrating it.
  for(int i=0; i<inMode3->GetMult(); i++){
    if(inMode3->GetHit(i)->GetTrace(0)->GetHole()!=31)
      hasmode3 = true;
    else
      card29trace = inMode3->GetHit(i)->GetTrace(0);
  }

  //in case of simulation, outS800 already contains the information needed.
  //outS800->Clear();
  outGret->Clear();
  outMode3->Clear();

  if(hass800)
    BuildS800Calc(inS800,outS800);
  if(hasgret){
    BuildGretinaCalc(inGret,outGret);
    fgretactr++;
    if(outS800->GetTimeStamp()>-1){
      if(fhasmap)
    	outGret->DopplerCorrect(fSett,outS800->GetTRACK());
      else
    	outGret->DopplerCorrect(fbeta);
    }
    if(card29trace != NULL){
      fcard29ctr++;
      outGret->CorrectTime(card29trace);
    }
  }
  if(hasmode3)
    BuildMode3Calc(inMode3,outMode3);

  fevent++;
}

void Calibration::BuildS800Calc(S800* in, S800Calc* out){
  //s800
  PAD pad[2];
  PPAC ppac[2];
  IC ich;
  TOF tof;
  SCINT scint[3];
  TRACK track;
  IITRACK iitrack;
  HODO hodo;
  out->Clear();
  ich.Clear();
  tof.Clear();
  track.Clear();
  iitrack.Clear();
  hodo.Clear();


  out->SetRegistr(in->GetTrigger()->GetRegistr());
  out->SetTimeS800(in->GetTrigger()->GetS800());
  out->SetTimeStamp(in->GetTS());
  SetTS800(in->GetTrigger()->GetS800());

  //CRDCs
  for(UShort_t k=0;k<2;k++){
    pad[k].Clear();
    SetPad(in->GetCrdc(k));
    pad[k] = GetPad();
  }

  if (out->GetRegistr()==2){
    static double ypos0offset = -fSett->YSlope(0)*fSett->YCoincHeadstart(0);
    pad[0].SetXY(pad[0].GetX(),pad[0].GetY()+ypos0offset);
    static double ypos1offset = -fSett->YSlope(1)*fSett->YCoincHeadstart(1);
    pad[1].SetXY(pad[1].GetX(),pad[1].GetY()+ypos1offset);
  }

  //PPACs
  PpacSort(in->GetTPpac()->GetChannels(), in->GetTPpac()->GetData(),&ppac[0],&ppac[1]);

  //IITRACK
  BuildIITrackCalc(ppac[0].GetX(), ppac[0].GetY(), ppac[1].GetX(), ppac[1].GetY(),
		   &iitrack);

  //TRACK
  BuildTrackCalc(pad[0].GetX(), pad[0].GetY(), pad[1].GetX(), pad[1].GetY(),
		 &track);

   //IC
  bool icgood = BuildIonChamber(in->GetIonChamber(),&track,&ich);


 //Hodo
  BuildHodoscope(in->GetHodoscope(), &hodo);

  //TOF
  BuildTof(in->GetTimeOfFlight(),&track,
	   &tof);

  //SCINTILLATOR
  bool scintgood[3] ={false,false,false};
  for(UShort_t s=0; s<3; s++){
    scintgood[s] = BuildScint(in->GetScintillator(s),&scint[s]);
  }

  //set counters
  //if(scintgood[0]){ //set the counters if the scintillator is good (above threshold)
  if(icgood){ //set the counters if the IC is good 
    for(UShort_t k=0;k<2;k++){
      if(!isnan(pad[k].GetX()))
	fpadctr[k]++;
      if(!isnan(ppac[k].GetX()) && !isnan(ppac[k].GetY()))
	fppacctr[k]++;
    }
    if(!isnan(iitrack.GetAII()))
      fiitrackctr++;
    if(!isnan(track.GetAFP()))
      ftrackctr++;
    if(tof.GetTACOBJ()>0)
      fobjctr++;
    if(tof.GetTACXFP()>0){
      fxfpctr++;
      if(tof.GetTACOBJ()>0)
	ftofctr++;
    }

  }

  //set Calculated S800
  out->SetPAD(pad[0],0);
  out->SetPAD(pad[1],1);
  out->SetPPAC(ppac[0],0);
  out->SetPPAC(ppac[1],1);
  out->SetIC(ich);
  out->SetTOF(tof);
  out->SetSCINT(scint[0],0);
  out->SetSCINT(scint[1],1);
  out->SetSCINT(scint[2],2);
  out->SetTRACK(track);
  out->SetIITRACK(iitrack);
  out->SetHODO(hodo);
}

void Calibration::BuildTrackCalc(Float_t x0, Float_t y0, Float_t x1, Float_t y1,
				 TRACK* out){
  out->Clear();
  if(isnan(x0)||isnan(y0)||isnan(x1)||isnan(y1))
    return;

  Float_t gap = fSett->Gap();
  Float_t pi = TMath::Pi();
  Float_t xsin, ysin, input[4], ratio, betagam, gamma, p0, e0;
  Short_t order = 5; //why?
  Float_t zfp = fSett->FPShift();
  betagam = fbrho/3.107/fmass*fcharge; //3.107 = 3.33547*0.9315
  //logbook page 4 Brho = 3.335 p in GeV/c
  gamma = sqrt(betagam*betagam+1);
  ratio = gamma/(gamma+1);
  p0 = fbrho*fcharge/3.33547; // in GeV/c   //3.33547 = 3.107/931.5*1000
  //p0 = betagam*m = betagam*fmass*0.9315
  e0 = 931.5016*fmass*(gamma-1); // in MeV
  if(order > fmaxorder)
    order = fmaxorder;

  //cout << "x0 " << x0 << "\ty0 " << y0 << "\tx1 " << x1 << "\ty1 " << y1 << endl;
  out->SetAFP(atan((x1 - x0)/gap) * 1000); // in mrad
  out->SetXFP(x0 + zfp * tan(out->GetAFP()/1000)); // in mm
  out->SetBFP(atan((y1 - y0) / gap) * 1000); // in mrad
  out->SetYFP(y0 + zfp * tan(out->GetBFP()/1000)); // in mm

  out->SetYFP(out->GetYFP()+fSett->ShiftY()); // in mm

  // Prepare input for COSY inverse map
  input[0] = -out->GetXFP()/1000;
  input[1] = -out->GetAFP()/1000;
  input[2] = out->GetYFP()/1000;
  input[3] = out->GetBFP()/1000;

  if(order > 0){
    out->SetATA(MapCalc(order, 0, input) * 1000); // in mrad
    out->SetYTA(MapCalc(order, 1, input) * 1000); // in mm
    out->SetBTA(MapCalc(order, 2, input) * 1000); // in mrad
    out->SetDTA(MapCalc(order, 3, input) * 100); // in % dE/E

    //cout << "out->fata " << out->fata << " out->fyta " << out->fyta << " out->fbta " << out->fbta << " out->fdta " << out->fdta << endl;
    if(fabs(out->GetATA()) > 100)
      out->SetATA(sqrt(-1.0));
    if(fabs(out->GetYTA()) > 100)
      out->SetYTA(sqrt(-1.0));
    if(fabs(out->GetBTA()) > 200)
      out->SetBTA(sqrt(-1.0));
    if(fabs(out->GetDTA()) > 10)
      out->SetDTA(sqrt(-1.0));

    out->SetATA(out->GetATA()+fSett->AngleA());
    out->SetBTA(out->GetBTA()+fSett->AngleB());
    xsin = sin(out->GetATA()/1000);
    ysin = sin(out->GetBTA()/1000);

    //azita azimutal angle at target = phi
    double s800_phi;
    //this accounts for the rotation of the coordinate system gretina - s800
    if (xsin > 0 && ysin > 0) {
      s800_phi = -atan(ysin/xsin) ;
    } else if (xsin < 0 && ysin > 0) {
      s800_phi = -pi + atan(ysin/fabs(xsin)) ;
    } else if (xsin < 0 && ysin < 0) {
      s800_phi = -pi - atan(fabs(ysin)/fabs(xsin));
    } else if (xsin > 0 && ysin < 0) {
      s800_phi = -2*pi + atan(fabs(ysin)/xsin);
    } else {
      s800_phi = 0.0 ;
    }
    out->SetPhi(s800_phi);


    //scatter polar angle at traget = theta
    out->SetTheta(asin(sqrt(xsin*xsin + ysin*ysin)));

    out->SetPtot(p0*(1.0 + out->GetDTA()/100*ratio));
    out->SetPpar(out->GetPtot()*cos(out->GetTheta()));
    out->SetPtra(out->GetPtot()*sin(out->GetTheta()));
    out->SetEtot(e0*(1.0 + out->GetDTA()/100));
  }
}
void Calibration::BuildIITrackCalc(Float_t x0, Float_t y0, Float_t x1, Float_t y1,
				   IITRACK* out){
  out->Clear();
  if(isnan(x0)||isnan(x1))
    return;

  Float_t gap = fSett->PpacGap();
  Float_t zii = fSett->PpacZ();
  Float_t pi = TMath::Pi();
  Float_t xsin, ysin;

  out->SetAII(atan((x1 - x0)/gap) * 1000); // in mrad
  out->SetXII(x0 + zii * tan(out->GetAII()/1000)); // in mm
  if(isnan(y0)||isnan(y1))
    return;
  out->SetBII(atan((y1 - y0) / gap) * 1000); // in mrad
  out->SetYII(y0 + zii * tan(out->GetBII()/1000)); // in mm
  xsin = sin(out->GetAII()/1000);
  ysin = sin(out->GetBII()/1000);

  //azita azimutal angle at intermediate image = phi
  if(xsin > 0 && ysin > 0)
    out->SetPhi(atan(ysin/xsin));
  else if(xsin < 0 && ysin > 0)
    out->SetPhi((pi - atan(ysin/fabs(xsin))));
  else if(xsin < 0 && ysin < 0)
    out->SetPhi((pi + atan(fabs(ysin)/fabs(xsin))));
  else if(xsin > 0 && ysin < 0)
    out->SetPhi((2*pi - atan(fabs(ysin)/xsin)));
  else
    out->SetPhi(0.0);

  //scatter polar angle at intermediate image = theta
  out->SetTheta(asin(sqrt(xsin*xsin + ysin*ysin)));
}

void Calibration::PpacSort(vector<Short_t>* channel, vector<Short_t>* data,
			   PPAC* out1, PPAC* out2){
  vector<vector<Short_t> > sum;
  sum.resize(4);
  for(UShort_t i=0;i<4;i++)
    sum[i].resize(fSett->PpacChannels());
  if(channel->size() != data->size()){
    cerr << " channel ("<<channel->size()<<") and data ("<<data->size()<<" have different sizes in PPAC" << endl;
    return;
  }
  for(UShort_t f=0;f<channel->size();f++){
    int det = channel->at(f)/64;
    int ch = channel->at(f) -det*64;
    if(det%2==0)
      ch = tppacmapx[ch];

    else
      ch = tppacmapy[ch];
    sum[det][ch] += data->at(f);
  }
  out1->Clear();
  out2->Clear();
  out1->SetStrips(sum[0],sum[1]);
  out2->SetStrips(sum[2],sum[3]);

  Short_t mult[4] = {0,0,0,0};
  Short_t maxvalue[4] = {0,0,0,0};
  Short_t maxstrip[4] = {-1,-1,-1,-1};

  for(Int_t i=0; i<fSett->PpacChannels(); i++){
    for(Int_t j=0; j<4; j++){
      if(sum[j][i]>0)
	mult[j]++;
      if(sum[j][i]>maxvalue[j]){
	maxvalue[j] = sum[j][i];
	maxstrip[j] = i;
      }
    }
  }
  out1->SetMult(mult[0],mult[1]);
  out1->SetMax(maxstrip[0],maxstrip[1]);
  out2->SetMult(mult[2],mult[3]);
  out2->SetMax(maxstrip[2],maxstrip[3]);

  float mean[4];

  for(Int_t j=0; j<4; j++){
    int low = (int)maxstrip[j] - (int)fSett->PpacWidth()/2;
    int high = low + (int)fSett->PpacWidth();
    if(low<0)
      low = 0;
    if(high>63)
      high = 63;
    float tmom =0;
    int tsum =0;
    for(int i=low;i<high;i++){
      tsum += sum[j][i];
      tmom += i*sum[j][i];
    }
    mean[j] = (tmom/tsum-32)*fSett->PpacPitch();
  }
  out1->SetXY(mean[0],mean[1]);
  out2->SetXY(mean[2],mean[3]);
}

void Calibration::CrdcCal(vector<Short_t>* channel, vector<Short_t>* data, int nr){
  Short_t index;
  vector<Float_t> sum;
  vector<Short_t> samples;
  sum.clear();
  sum.resize(fSett->CrdcChannels());
  samples.clear();
  samples.resize(fSett->CrdcChannels());
  fcrdccal.clear();
  fcrdccal.resize(fSett->CrdcChannels());
  if(channel->size() != data->size()){
    cerr << " channel ("<<channel->size()<<") and data ("<<data->size()<<" have different sizes " << endl;
    return;
  }

  for(UShort_t f=0;f<channel->size();f++){
    index = channel->at(f);
    sum[index] += data->at(f) - fped[nr][index];
    samples[index]++;
    if(fSett->VLevel()>6)
      cout << "sum["<<index<<"] " << sum[index] << " samples["<<index<<"] " << samples[index] << endl;
  }
  for(UShort_t ch=0; ch<fSett->CrdcChannels(); ch++){
    if(samples[ch] > 0){//read minsamples from file?
#ifdef SPECTCL_MODE
      fcrdccal[ch]=sum[ch]/fSett->SampleWidth();  // spectcl compatability
#else
      fcrdccal[ch]=sum[ch]/samples[ch];
#endif
      fcrdccal[ch]*=fslope[nr][ch];
      fcrdccal[ch]+=foffset[nr][ch];
    }
    else
      fcrdccal[ch]=sqrt(-1.0);
  }
}

void Calibration::SetPad(GCrdc* in){
  int id = in->GetID();
  fpad.Clear();
  fpad.SetID(id);
  //cout << "processing pad " << fpad.GetID() << endl;
  this->CrdcCal(in->GetChannels(),in->GetData(),id);
  fpad.SetCal(fcrdccal);
  double x = sqrt(-1.0);
  if(fSett->Method(id)==0){
    x = fSett->XOffset(id) + fSett->XSlope(id) * this->CalcGravity();
  } else if(fSett->Method(id)==1){
    x = fSett->XOffset(id) + fSett->XSlope(id) * this->CalcGravityClassic();
  } else if(fSett->Method(id)==2){
    x = fSett->XOffset(id) + fSett->XSlope(id) * this->CalcFit();
  } else if(fSett->Method(id)==3){
    x = fSett->XOffset(id) + fSett->XSlope(id) * this->CalcGravity() + 0.0* this->CalcFit();
  } else {
    cout << "unknown method number" << endl;
  }

  double y =  sqrt(-1.0);
  float tac = in->GetTAC();
  if(tac>0 && tac <4000){
    y= fSett->YOffset(id) + fSett->YSlope(id)*in->GetTAC();
    //y depends a bit on x, this can be corrected
    if(!isnan(x))
      y -= (x*fSett->YCorSlope(id)+fSett->YCorOffset(id));
    double y_cor = 0;
    //cout << "fevent\t" << fevent << " id " << id << " inte "<< fy_cor[id]->Integral() << endl;
    if(fy_cor[id]!=NULL)
      y_cor = fy_cor[id]->GetBinContent(fevent/10000+1);
    //cout << "crdc " << id << " y_cor " << y_cor << endl;
    y+=y_cor;

  }
  fpad.SetXY(x,y);

  fpad.SetTAC(tac);
}

Float_t Calibration::CalcGravity(){
  if(fSett->VLevel()>2)
    cout << " start of " << __PRETTY_FUNCTION__  << endl;
  int crdcnum = fpad.GetID();
  int maxpad = -1; // maximum pad
  Float_t padmax = -1000; //value of that
  //Find pad with maximum value.
  //Both adjacent pads must have valid entries.
  for(UShort_t i=0;i<fSett->CrdcChannels();i++){
    if(!IsBad(crdcnum,i) && !isnan(fcrdccal[i]) && fcrdccal[i] > padmax){
      if( fcrdccal[i-1]>0 && fcrdccal[i+1]>0 && i!=fpad.GetChan().front() && i!=fpad.GetChan().back()){

	maxpad = i;
	padmax = fcrdccal[i];
      } // left and right pad ok
    }
  }
  //Verbose printouts.
  if(fSett->VLevel()>2){
    cout <<"pad " << crdcnum << "determined maximum " << maxpad << " with value " << padmax << endl;
    //cout << crdcnum << " values found for chan " << fpad.GetChan().front() << " to " << fpad.GetChan().back() << endl;
    if(fSett->VLevel()>6){
      for(UShort_t i=0;i<fpad.GetChan().size();i++)
	cout << "chan " << fpad.GetChan(i) << " value " << fpad.GetCal(i) << endl;
    }
  }

  //Bail out if no maximum has been found.
  if(maxpad<0)
    return sqrt(-1.0);

  int l = maxpad;
  int r = maxpad;
  Float_t padleft = sqrt(-1.0);
  Float_t padright = sqrt(-1.0);
  if(!IsBad(crdcnum,maxpad-1)&& !isnan(fcrdccal[maxpad-1]))
    padleft = fcrdccal[maxpad-1];
  if(!IsBad(crdcnum,maxpad+1)&& !isnan(fcrdccal[maxpad+1]))
    padright = fcrdccal[maxpad+1];

  //saturated pads
  while(1){
    while(l>fpad.GetChan().front() && fcrdccal[l] >= fSett->Saturation(crdcnum)*fSett->MSaturation(crdcnum) ){
      l--;
    }
    if(IsBad(crdcnum,l)&& l>fpad.GetChan().front() && fcrdccal[l-1] >= fSett->Saturation(crdcnum)*fSett->MSaturation(crdcnum))
      l--;
    else
      break;
  }
  //cout << " end while(1) l " << l << endl;
  while(1){
    while(r<fpad.GetChan().back() && fcrdccal[r] >= fSett->Saturation(crdcnum)*fSett->MSaturation(crdcnum) ){
      r++;
    }
    if(IsBad(crdcnum,r) && r<fpad.GetChan().back() && fcrdccal[r+1] >= fSett->Saturation(crdcnum)*fSett->MSaturation(crdcnum)) //(dirk : r+1<fpad.GetChan()[fpad.GetChan().size()-1] ??)
      r++;
    else
      break;
  }
  if(l!=r && fSett->VLevel()>2)
    cout << "l " << l << " r " << r << endl;

  maxpad = (l+r)/2;
  fpad.SetMaxPad(maxpad);
  fpad.SetPadMax(0,padleft);
  fpad.SetPadMax(1,padmax);
  fpad.SetPadMax(2,padright);

  int lowpad = (int)maxpad - (int)fSett->GravityWidth(crdcnum)/2;
  int highpad = lowpad + (int)fSett->GravityWidth(crdcnum);
  if(fSett->VLevel()>2)
    cout << lowpad << "\t" << highpad << endl;

  if(lowpad < fpad.GetChan().front())
    lowpad = fpad.GetChan().front();
  if(highpad >= fpad.GetChan().back()+1)
    highpad = fpad.GetChan().back()+1;
  if(fSett->VLevel()>2)
    cout << "corrected to " << lowpad << "\t" << highpad << endl;

  Float_t sum = 0, mom = 0;
  Short_t padused =0;
  for(UShort_t i=lowpad;i<highpad;i++){
    Float_t padamp;
    if(IsBad(crdcnum,i)){
      padamp = 0;
      if(!IsBad(crdcnum,i-1) && i>fpad.GetChan().front() && !isnan(fcrdccal[i-1]) && fcrdccal[i-1]>0)
	padamp += (int)(0.5*fcrdccal[i-1]);
      if(!IsBad(crdcnum,i+1) && i<fpad.GetChan().back() && !isnan(fcrdccal[i+1]) && fcrdccal[i+1]>0)
	padamp += (int)(0.5*fcrdccal[i+1]);
    }
    else if(!isnan(fcrdccal[i]) && fcrdccal[i]>0)
      padamp = fcrdccal[i];
    else
      padamp =0;

    if(padamp>fSett->Saturation(crdcnum)){
      //cout << "never here " << endl;
      padamp = fSett->Saturation(crdcnum);
    }
    sum +=     (int) padamp; //why that?
    mom += i * (int) padamp;

    if(padamp>1)
      padused++;

  } //lowpad .. highpad
  //if( crdcnum==0)
  Float_t x = (Float_t)mom/sum;

  if(x < 0 || x > fSett->CrdcChannels()){ // no gravity center found
    x = sqrt(-1.0);
  }
  fpad.SetXGravity(x);
  return x;

}


Float_t Calibration::CalcGravityClassic(){
  int crdcnum = fpad.GetID();
  int maxpad = -1; // maximum pad
  Float_t padmax = -1000; //value of that
  for(UShort_t i=0;i<fSett->CrdcChannels();i++){
    //find pad with maximum
    if(!isnan(fcrdccal[i]) && fcrdccal[i] > padmax){
      //neighboring pad also valid
      if( (!isnan(fcrdccal[i-1])&&fcrdccal[i-1]>0) && (!isnan(fcrdccal[i+1])&&fcrdccal[i+1]>0) ){

	maxpad = i;
	padmax = fcrdccal[i];
      } // left and right pad ok
    }
  }
  /*
    //without left right check
  for(UShort_t i=0;i<fpad.GetChan().size();i++){
    //cout <<padmax << "\t" << fpad.GetChan()[i] << "\t" << fpad.GetCal()[i] << endl;
    if(!isnan(fpad.GetCal()[i]) && fpad.GetCal()[i] > padmax){
      maxpad = fpad.GetChan()[i];
      padmax = fpad.GetCal()[i];
    }
  }
  */
  if(maxpad<0)// no max pad found
    return sqrt(-1.0);
  if(fSett->VLevel()>2)
    cout << "determined maximum " << maxpad << " with value " << padmax << endl;
  int lowpad = (int)maxpad - (int)fSett->GravityWidth(crdcnum)/2;
  int highpad = lowpad + (int)fSett->GravityWidth(crdcnum) + 1;
  if(lowpad < 0)
    lowpad = 0;
  if(highpad >= fSett->CrdcChannels())
    highpad = fSett->CrdcChannels() - 1;
  if(fSett->VLevel()>2)
    cout << lowpad << "\t" << highpad << endl;
  Float_t sum = 0, mom = 0;// mom2 = 0, mom3 = 0;
  for(UShort_t i=0;i<fpad.GetChan().size();i++){
    if(fpad.GetChan()[i] >= lowpad && fpad.GetChan()[i] <= highpad){
      sum+=fpad.GetCal()[i];
      mom+=fpad.GetChan()[i]*fpad.GetCal()[i];
    }//get the ones between low and highpad
  }// loop over all channels with hits
  Float_t x = (Float_t)mom/sum;

  if(x < 0 || x > fSett->CrdcChannels()){ // no gravity center found
    x = sqrt(-1.0);
  }
  fpad.SetXGravity(x);
  return x;
}
Float_t Calibration::CalcFit(){
  if(fSett->VLevel()>2)
    cout << " start of " << __PRETTY_FUNCTION__  << endl;
  int crdcnum = fpad.GetID();
  int maxpad = -1; // maximum pad
  Float_t padmax = -1000; //value of that
  //find pad with maximum
  for(UShort_t i=0;i<fSett->CrdcChannels();i++){
    if(!IsBad(crdcnum,i) && !isnan(fcrdccal[i]) && fcrdccal[i] > padmax){
      if( fcrdccal[i-1]>0 && fcrdccal[i+1]>0 && i!=fpad.GetChan().front() && i!=fpad.GetChan().back()){

	maxpad = i;
	padmax = fcrdccal[i];
      } // left and right pad ok
    }
  }
  if(fSett->VLevel()>2){
    cout <<"pad " << crdcnum << "determined maximum " << maxpad << " with value " << padmax << endl;
    //cout << crdcnum << " values found for chan " << fpad.GetChan()[0] << " to " << fpad.GetChan()[fpad.GetChan().size()-1] << endl;
    for(UShort_t i=0;i<fpad.GetChan().size();i++)
      cout << "chan " << fpad.GetChan()[i] << " value " << fpad.GetCal()[i] << endl;
  }

  if(maxpad<0)// no max pad found
    return sqrt(-1.0);
  int l = maxpad;
  int r = maxpad;
  Float_t padleft = sqrt(-1.0);
  Float_t padright = sqrt(-1.0);
  if(!IsBad(crdcnum,maxpad-1)&& !isnan(fcrdccal[maxpad-1]))
    padleft = fcrdccal[maxpad-1];
  if(!IsBad(crdcnum,maxpad+1)&& !isnan(fcrdccal[maxpad+1]))
    padright = fcrdccal[maxpad+1];

  //saturated pads
  while(1){
    while(l>fpad.GetChan().front() && fcrdccal[l] >= fSett->Saturation(crdcnum)*fSett->MSaturation(crdcnum) ){
      l--;
    }
    //cout << " l " << l << endl;
    if(IsBad(crdcnum,l)&& l>fpad.GetChan().front() && fcrdccal[l-1] >= fSett->Saturation(crdcnum)*fSett->MSaturation(crdcnum))
      l--;
    else
      break;
    //cout << "bad pad check l " << l << endl;

  }
  //cout << " end while(1) l " << l << endl;
  while(1){
    while(r<fpad.GetChan().back() && fcrdccal[r] >= fSett->Saturation(crdcnum)*fSett->MSaturation(crdcnum) ){
      r++;
    }
    if(IsBad(crdcnum,r) && r<fpad.GetChan().back() && fcrdccal[r+1] >= fSett->Saturation(crdcnum)*fSett->MSaturation(crdcnum)) //(dirk : r+1<fpad.GetChan()[fpad.GetChan().size()-1] ??)
      r++;
    else
      break;
  }
  if(l!=r && fSett->VLevel()>2)
    cout << "l " << l << " r " << r << endl;

  maxpad = (l+r)/2;
  fpad.SetMaxPad(maxpad);
  fpad.SetPadMax(0,padleft);
  fpad.SetPadMax(1,padmax);
  fpad.SetPadMax(2,padright);

  int lowpad = (int)maxpad - (int)fSett->GravityWidth(crdcnum)/2;
  int highpad = lowpad + (int)fSett->GravityWidth(crdcnum);
  if(fSett->VLevel()>2)
    cout << lowpad << "\t" << highpad << endl;

  if(lowpad < fpad.GetChan().front())
    lowpad = fpad.GetChan().front();
  if(highpad >= fpad.GetChan().back()+1)
    highpad = fpad.GetChan().back()+1;
  if(fSett->VLevel()>2)
    cout << "corrected to " << lowpad << "\t" << highpad << endl;


  const int MA=4;
  int j, mfit=MA;
  double chisq, *x, *y, *sig, *a;
  int NPT = highpad - lowpad + 1; // Maximum number of points
  a = (double*)malloc((MA+1)*sizeof(double));
  x = (double*)malloc((NPT+1)*sizeof(double));
  y = (double*)malloc((NPT+1)*sizeof(double));
  sig = (double*)malloc((NPT+1)*sizeof(double));
  j = 0;
  for(UShort_t i=lowpad;i<highpad;i++){
    Float_t padamp;
    if(IsBad(crdcnum,i)){
      padamp = 0;
      if(!IsBad(crdcnum,i-1) && i>fpad.GetChan().front() && !isnan(fcrdccal[i-1]) && fcrdccal[i-1]>0)
	padamp += (int)(0.5*fcrdccal[i-1]);
      if(!IsBad(crdcnum,i+1) && i<fpad.GetChan().back() && !isnan(fcrdccal[i+1]) && fcrdccal[i+1]>0)
	padamp += (int)(0.5*fcrdccal[i+1]);
    }
    else if(!isnan(fcrdccal[i]) && fcrdccal[i]>0)
      padamp = fcrdccal[i];
    else
      padamp =0;

    if(padamp>fSett->Saturation(crdcnum)){
      //cout << "never here " << endl;
      padamp = fSett->Saturation(crdcnum);
    }
    if(padamp>0){
      j++;
      x[j] = (double)i;
      y[j] = padamp;
      sig[j] = 50.0;

    }


  }

  // reset number of points
  int Npoints = j;
  // forget it if there are too few points
  if(Npoints < MA){
    free(a);
    free(x);
    free(y);
    free(sig);
    return sqrt(-1.0);
  }
// guess for fit starting point
  a[1] = padmax;
  a[2] = maxpad;
  a[3] = 3.0;
  a[4] = 0;

  Float_t xf;
  if (Fit(x, y, sig, a, mfit, Npoints, &chisq)) {
    //cout << "success! " << a[2] << "\t" << chisq<< endl;
    xf = a[2];
    //x_chi2 = chisq;
  }
  else{
    cout << "fitting error " << a[2] << "\t" << chisq<< endl;
    xf = sqrt(-1.0);
  }
  free(a);
  free(x);
  free(y);
  free(sig);
  fpad.SetXFit(xf);
  return xf;

}

bool Calibration::BuildIonChamber(GIonChamber* in, TRACK* track, IC* out){
  out->Clear();
  out->SetCal(ICCal(in));
  double sum = ICSum(out->GetCal());
  double IC_corg = 1;
  double IC_coro = 0;
  if(fIC_cor[0]!=NULL)
    IC_corg = fIC_cor[0]->GetBinContent(fevent/10000+1);
  if(fIC_cor[1]!=NULL)
    IC_coro = fIC_cor[1]->GetBinContent(fevent/10000+1);
  sum *= IC_corg;
  sum += IC_coro;
  out->SetSum(sum);
  out->SetDE(ICDE(sum,track));
  if(!isnan(sum) && sum>fSett->ICThresh()){
    fichctr++;
    return true;
  }
  return false;
}

void Calibration::BuildHodoscope(GHodoscope* in, HODO* out){
  out->Clear();
  out->SetTime(in->GetTime());
  bool fired = false;
  vector<Double_t> en;
  vector<Double_t> ch;
  for(UShort_t c=0;c<in->GetData()->size();c++){
    if(in->GetData()->at(c)>4050)
      continue;
    double cal = in->GetData()->at(c)-0.5+fRand->Uniform(0,1);
    UShort_t chv = in->GetChannels()->at(c);
    cal*=fhodogain[chv];
    cal+=fhodooffset[chv];
    if(cal>fSett->HodoThresh()){
      out->Set(cal,chv);
      en.push_back(cal);
      ch.push_back(chv);
      fired = true;
    }
  }
  if(fired) fhodoctr++;
  //addback for hodoscope
  while(ch.size()>0){
    //first is good
    out->AddFirstAB(en.back(),ch.back());
    //cout << "added first " << ch.back() << "\t" << en.back() << endl;
    vector<Short_t> usedchannels;
    usedchannels.push_back(ch.back());
    en.pop_back();
    ch.pop_back();
    while(usedchannels.size()>0){
      Short_t neigh = -1;
      for(UShort_t i=0;i<ch.size();i++){
	if(HodoNeighbor(usedchannels.back(),ch[i])){
	  neigh = i;
	  break;
	}
      }
      if(neigh == -1)
	usedchannels.pop_back();
      else{
	out->AddAB(en[neigh],ch[neigh]);
	usedchannels.push_back(ch[neigh]);
	en.erase(en.begin()+neigh);
	ch.erase(ch.begin()+neigh);
      }
    }
  }
}
bool Calibration::HodoNeighbor(Short_t ch1, Short_t ch2){
  if(ch1-4 == ch2 || ch1+4 == ch2 )//above below
    return true;
  if((ch1+1 == ch2) && ch2/4 == (ch1+1)/4) // right/samerow
    return true;
  if((ch1-1 == ch2) && ch2/4 == (ch1-1)/4) // left/samerow
    return true;
  if((ch1-3 == ch2) && ch2/4 == (ch1-3)/4) // right/uprow
    return true;
  if((ch1-5 == ch2) && ch2/4 == (ch1-5)/4) // left/uprow
    return true;
  if((ch1+5 == ch2) && ch2/4 == (ch1+5)/4) // right/lowrow
    return true;
  if((ch1+3 == ch2) && ch2/4 == (ch1+3)/4) // left/lowrow
    return true;
  
  return false;
  
}
vector<Float_t> Calibration::ICCal(GIonChamber* in){
  vector<int>* chan = in->GetChannels();
  vector<float>* raw = in->GetData();
  vector<Float_t> cal;
  cal.resize(fSett->IonChamberChannels());
  if(chan->size() != raw->size()){
    cerr << " channel ("<<chan->size()<<") and data ("<<raw->size()<<" have different sizes " << endl;
    return cal;
  }
  for(unsigned int f=0;f<chan->size();f++){
    if( (chan->at(f)>-1) && (chan->at(f)<fSett->IonChamberChannels()) ){
      cal[chan->at(f)] = raw->at(f)*fICslope[chan->at(f)]+fICoffset[chan->at(f)];
    }
    else{
      cerr << " channel "<<chan->at(f)<<" not found!" << endl;

    }
  }
  return cal;
}
Float_t Calibration::ICSum(vector<Float_t> cal){
  Short_t ch = 0;
  Float_t sum =0;
  for(UShort_t j=0; j<cal.size(); j++){
    if(cal[j]>0){
      sum += cal[j];
      ch++;
    }
  }
  if(ch > 0)
    sum/= ch;
  else
    sum = sqrt(-1.0);

  return sum;

}
Float_t Calibration::ICDE(Float_t sum, TRACK* track){
  Float_t x = track->GetXFP();
  Float_t y = track->GetYFP();
  if(!isnan(sum) && !isnan(track->GetAFP())){
    if(!isnan(y))
      sum += sum*fSett->dE_ytilt()*y;
    if(!isnan(x) && x < fSett->dE_x0tilt())
      sum *= exp(fSett->dE_xtilt()* (fSett->dE_x0tilt() -x) );
    fs800valid = 0;
    return sum * fde_slope + fde_offset;
  } else {
    return sqrt(-1.0);
  }
}

bool Calibration::BuildScint(GScintillator* in, SCINT* out){
  out->Clear();
  out->SetTime(TimeOffset(in->GetTime_up()),
	       TimeOffset(in->GetTime_down()));
  out->SetDE(in->GetDE_up(),
	     in->GetDE_down());
  if(out->GetDE()>fSett->GetScintThresh()){
    fscintctr++;
    return true;
  }
  return false;
}

Float_t Calibration::TimeOffset(Float_t time1, Float_t time2){
  return time1 - time2;
}
Float_t Calibration::TimeOffset(Float_t time){
  return time - fts800;
}
Float_t Calibration::RFCorrection(TRACK* tr){
  return fSett->RFe1C()*tr->GetAFP() + fSett->RFC()*tr->GetXFP();
}
Float_t Calibration::OBJCorrection(TRACK* tr){
  return fSett->OBJe1C()*tr->GetAFP() + fSett->OBJC()*tr->GetXFP();
}
Float_t Calibration::XFPCorrection(TRACK* tr){
  return fSett->XFPe1C()*tr->GetAFP() + fSett->XFPC()*tr->GetXFP();
}
Float_t Calibration::TACOBJCorrection(TRACK* tr){
  return fSett->TACOBJe1C()*tr->GetAFP() + fSett->TACOBJC()*tr->GetXFP();
}
Float_t Calibration::TACXFPCorrection(TRACK* tr){
  return fSett->TACXFPe1C()*tr->GetAFP() + fSett->TACXFPC()*tr->GetXFP();
}

void Calibration::BuildTof(GTimeOfFlight* in, TRACK* track, TOF* out){
  out->Clear();
  double obj_cor = 0;
  if(fobj_cor!=NULL)
    obj_cor = fobj_cor->GetBinContent(fevent/10000+1);
  double xfp_cor = 0;
  if(fxfp_cor!=NULL)
    xfp_cor = fxfp_cor->GetBinContent(fevent/10000+1);
  double objtac_cor = 0;
  if(fobjtac_cor!=NULL)
    objtac_cor = fobjtac_cor->GetBinContent(fevent/10000+1);
  double xfptac_cor = 0;
  if(fxfptac_cor!=NULL)
    xfptac_cor = fxfptac_cor->GetBinContent(fevent/10000+1);

  out->Set(TimeOffset(in->GetRF()), TimeOffset(in->GetOBJ())+obj_cor, TimeOffset(in->GetXFP())+xfp_cor);
  out->SetCorr(RFCorrection(track), OBJCorrection(track), XFPCorrection(track));
  out->SetTAC(in->GetTACOBJ()+objtac_cor, in->GetTACXFP()+xfptac_cor);
  out->SetTACCorr(TACOBJCorrection(track), TACXFPCorrection(track));


}

void Calibration::ReadNeighbors(const char *filename){
  TEnv *file = new TEnv(filename);
  for(int d=0;d<7;d++){
    for(int c=0;c<4;c++){
      fnrofneigh[d][c] = file->GetValue(Form("Gretina.Neighbors.NumNeigh.%d",fSett->Det2Hole(d)*4+c),0);
      for(int n=0;n<fnrofneigh[d][c];n++){
	fneighbors[d][c][n] = file->GetValue(Form("Gretina.Neighbors.%d.%d",fSett->Det2Hole(d)*4+c,n),0);
      }
    }
  }
  if(fSett->VLevel()>0){
    for(int d = 0; d<7; d++){
      for (int c = 0; c<4; c++){
	cout << "Detector " << d << " " << c << " has " << fnrofneigh[d][c] << " neighbors." << endl;
	cout << "They are ";
	for (int n=0;n<6; n++){
	  cout << fneighbors[d][c][n] << " ";
	}
	cout << endl;
      }
    }
  }
  return;
}
void Calibration::ReadGretinaCal(const char *filename){
  TEnv *cal = new TEnv(filename);
  for(int d=0;d<7;d++){
    for(int c=0;c<4;c++){
      fgslope[d][c] = cal->GetValue(Form("Slope.d%d.c%d",d,c),1.0);
      fgoffset[d][c] = cal->GetValue(Form("Offset.d%d.c%d",d,c),0.0);
    }
  }
}

bool Calibration::IsNeighbor(int d1, int c1, int d2, int c2){
  //check if d2 neighbors d1
  for(int i=0;i<fnrofneigh[d1][c1];i++){
    if(fneighbors[d1][c1][i] == fSett->Det2Hole(d2)*4 + c2)
      return true;
  }
  return false;
}
bool Calibration::IsNeighbor(int ID1, int ID2){
  return IsNeighbor(fSett->Hole2Det(ID1/4), ID1%4, fSett->Hole2Det(ID2/4), ID2%4);
}
int Calibration::NumNeighbors(int det, int cry){
  return fnrofneigh[fSett->Det2Hole(det)][cry];
}

bool Calibration::BuildMode3Calc(Mode3Event* in, Mode3Calc* out){
  out->Clear();
  for(int i=0; i<in->GetMult(); i++){
    Trace* trace = in->GetHit(i)->GetCoreTrace();
    if(trace!=NULL){
      bool status = Mode3Calibration(trace, out);
      if (!status){
	return false;
      }
    }
  }

  if (fAddBackType == 0){
    AddBackMode3Hole(out);
  } else if (fAddBackType == 1){
    AddBackMode3Neighbors(out);
  } else if (fAddBackType == 2){
    AddBackMode3Everything(out);
  }

  return true;
}

bool Calibration::Mode3Calibration(Trace* trace, Mode3Calc* out){
  int hole = trace->GetHole();
  int cry = trace->GetCrystal();
  int id = hole*4+cry;

  long long int time = trace->GetLED();
  float en = trace->GetEnergy()+fRand->Uniform(0,1);

  if(en<0){
    if(fverbose >= 1){
      cout << "trace with negative uncalibrated energy " << en << " being thrown out" << endl;
    }
    return false;
  }
  en = fgslope[fSett->Hole2Det(hole)][cry]*en+fgoffset[fSett->Hole2Det(hole)][cry];
  if(fverbose >= 1){
    cout << "making hit with energy " << en << ", time " << time << " in ID " << id << " (Hole: " << id/4 << ", Detector: " << fSett->Hole2Det(id/4) << ", crys: " << id%4 << ")" << endl;
  }

  //Card29 events should not be written to calibrated file.
  //These are for calibrating the S800, not for their own merits.
  if(trace->GetHole()!=31){
    out->AddHit(new Mode3Detection(en, time, id));
  }
  return true;
}

float Calibration::CalibratedEnergy(Trace* trace){
  int hole = trace->GetHole();
  int cry = trace->GetCrystal();
  float en = trace->GetEnergy()+fRand->Uniform(0,1);
  en = fgslope[fSett->Hole2Det(hole)][cry]*en+fgoffset[fSett->Hole2Det(hole)][cry];
  return en;
}

void Calibration::AddBackMode3Hole(Mode3Calc* mode3){
  //All hits within a cluster = holenumber are summed.
  for (int i=0; i<mode3->GetMult(); i++){
    bool addbacked = false;
    for (int j=0; j<mode3->GetMultAB(); j++){
      if (mode3->GetHit(i)->GetHole() == mode3->GetHitAB(j)->GetHole() ){
	mode3->GetHitAB(j)->AddBackHit(mode3->GetHit(i));
	addbacked = true;
      }
    }
    if (!addbacked){
      mode3->AddHitAB(new Mode3Detection(mode3->GetHit(i)));
    }
  }
}

void Calibration::AddBackMode3Neighbors(Mode3Calc* mode3){
  //Depth-first search through crystals to find all neighbors.
  //All hits within a region of adjacent hits should be summed.
  vector<Mode3Detection*> unusedHits = mode3->GetHits();
  while (unusedHits.size() > 0){
    Mode3Detection* currentHit = new Mode3Detection(unusedHits.back());
    unusedHits.pop_back();
    mode3->AddHitAB(currentHit);
    vector<int> IDsToCheck;
    IDsToCheck.push_back(currentHit->GetID());
    while(IDsToCheck.size() > 0){
      int IDchecking = IDsToCheck.back();
      int neighbor = -1;
      for (uint i=0; i<unusedHits.size(); i++){
	if (IsNeighbor(IDchecking, unusedHits[i]->GetID())){
	  neighbor = i;
	  break;
	}
      }
      if (neighbor == -1){
	IDsToCheck.pop_back();
      } else {
	currentHit->AddBackHit(unusedHits[neighbor]);
	IDsToCheck.push_back(unusedHits[neighbor]->GetID());
	unusedHits.erase(unusedHits.begin() + neighbor);
      }
    }
  }
}

void Calibration::AddBackMode3Everything(Mode3Calc* mode3){
  //All hits within Gretina are summed.
  if (mode3->GetMult() < 1)
    return;
  mode3->AddHitAB(new Mode3Detection(mode3->GetHit(0)));
  for (int i=1; i<mode3->GetMult(); i++){
    mode3->GetHitAB(0)->AddBackHit(mode3->GetHit(i));
  }
}

void Calibration::BuildGretinaCalc(Gretina* in, GretinaCalc* out){
  //Perform any recalibrations necessary.
  GretinaRecalibrate(in);

  //For no add-back, we don't need to do anything beyond copying relevant parameters.
  out->Clear();
  if (in->GetMult()==0){
    return;
  }
  for(int i=0; i<in->GetMult(); i++){
    if(in->GetHit(i)->GetError()>0)
      continue;
    CalibrateIPoints(in->GetHit(i));
    out->AddHit(new HitCalc(in->GetHit(i)));
  }
  //Perform the addback, which fills the add-backed vector.
  if (fAddBackType == -1){
    AddBackGretinaCrystal(out);
  } else if (fAddBackType == 0){
    AddBackGretinaHole(out);
  } else if (fAddBackType == 1){
    AddBackGretinaNeighbors(out);
  } else if (fAddBackType == 2){
    AddBackGretinaEverything(out);
  } else if (fAddBackType == 3){
    AddBackGretinaNeighborsNoSearch(out);
  } else if (fAddBackType == 4){
    ClusterGretina(out,in);
  } else {
    cout << "unknown addback type: " << fAddBackType << endl;
  }
  if(fSett->StoreAllIPoints())
    AllGretinaHits(out,in);

}

bool Calibration::GretinaRecalibrate(Gretina* gr){
  //Check to see if any hits were in crystals marked to be recalibrated.
  //If so, apply the recalibration given.
  for(int i=0; i<gr->GetMult(); i++){
    Crystal* crys = gr->GetHit(i);
    int det = fSett->Hole2Det(crys->GetHole());
    int crysnum = crys->GetCrystal();
    for(vector<recal>::iterator it = fGretRecal.begin(); it!=fGretRecal.end(); it++){
      if(det==it->Detector && crysnum==it->Crystal){
	crys->SetEnergy(crys->GetEnergy()*it->Slope + it->Offset);
	for(int j=0; j<crys->GetMult(); j++){
	  IPoint* ipoint = crys->GetIPoint(j);
	  ipoint->SetEnergy(ipoint->GetEnergy()*it->Slope + it->Offset);
	}
	break;
      }
    }
  }

  return true;
}

TVector3 Calibration::TransformCoordinates(int hole, int cry, TVector3 local){
  /* Need to convert from mm to cm for this to actually work properly. */
  double x = local.X()/10;
  double y = local.Y()/10;
  double z = local.Z()/10;
  double xt = fcrmat[hole][cry][0][0] * x + fcrmat[hole][cry][0][1] * y + fcrmat[hole][cry][0][2] * z + fcrmat[hole][cry][0][3];
  double yt = fcrmat[hole][cry][1][0] * x + fcrmat[hole][cry][1][1] * y + fcrmat[hole][cry][1][2] * z + fcrmat[hole][cry][1][3];
  double zt = fcrmat[hole][cry][2][0] * x + fcrmat[hole][cry][2][1] * y + fcrmat[hole][cry][2][2] * z + fcrmat[hole][cry][2][3];
  xt*=10.;
  yt*=10.;
  zt*=10.; //in mm
  return TVector3(xt,yt,zt) - fSett->TargetPos();
}
void Calibration::CalibrateIPoints(Crystal* cry){
  double sum =0;
  for(int j=0; j<cry->GetMult(); j++){
    IPoint* ipoint = cry->GetIPoint(j);
    ipoint->SetPosition(TransformCoordinates(cry->GetHole(),
					     cry->GetCrystal(),
					     ipoint->GetPosition()));
    sum+=ipoint->GetEnergy();
  }
  double core = cry->GetEnergy();
  for(int j=0; j<cry->GetMult(); j++){
    IPoint* ipoint = cry->GetIPoint(j);
    double en=ipoint->GetEnergy();
    en*=core/sum;
    ipoint->SetEnergy(en);
  }
}
vector<HitCalc*> Calibration::FilterOverflows(GretinaCalc* gr){
  vector<HitCalc*> output;
  for(int i=0; i<gr->GetMult(); i++){
    HitCalc* hit = gr->GetHit(i);
    if(hit->GetEnergy()<fSett->OverflowThreshold()){
      output.push_back(hit);
    }
  }
  return output;
}


vector<HitCalc*> Calibration::ExtractAllHits(Gretina* in){
  vector<HitCalc*> output;
  for(int i=0; i<in->GetMult(); i++){
    Crystal* cry = in->GetHit(i);
    if(cry->GetEnergy()<fSett->OverflowThreshold()){
      Short_t hole = cry->GetHole();
      Short_t crystal = cry->GetCrystal();
      long long int ts = cry->GetTS();
      Float_t t0 = cry->GetT0();
      Float_t chisq =  cry->GetChiSq();
      for(UShort_t j=0;j<cry->GetIPoints().size();j++){
	Float_t en = cry->GetIPoints()[j]->GetEnergy();
	TVector3 pos = cry->GetIPoints()[j]->GetPosition();
	output.push_back(new HitCalc(hole,crystal,en,ts,pos,t0,chisq));
      }
    }
  }
  return output;
}

void Calibration::AddBackGretinaCrystal(GretinaCalc* gr){
  vector<HitCalc*> filtered = FilterOverflows(gr);
  for(vector<HitCalc*>::iterator iter = filtered.begin(); iter!=filtered.end(); iter++){
    gr->AddHitAB(new HitCalc(*iter));
  }
}

void Calibration::AddBackGretinaHole(GretinaCalc* gr){
  //All hits within a cluster = holenumber are summed.
  vector<HitCalc*> filtered = FilterOverflows(gr);
  for(vector<HitCalc*>::iterator iter = filtered.begin(); iter!=filtered.end(); iter++){
    HitCalc* hit = *iter;
    bool addbacked = false;
    for (int j=0; j<gr->GetMultAB(); j++){
      if (gr->GetHitAB(j)->GetHole() == hit->GetHole()){
	gr->GetHitAB(j)->AddBackHitCalc(hit);
	addbacked = true;
	break;
      }
    }
    if (!addbacked){
      gr->AddHitAB(new HitCalc(*hit));
    }
  }
}

void Calibration::AddBackGretinaNeighbors(GretinaCalc* gr){
  //Depth-first search through crystals to find all neighbors.
  //All hits within a region of adjacent hits should be summed.
  vector<HitCalc*> unusedHits = FilterOverflows(gr);
  map<int,int> locations;
  int currIndex = -1;
  while (unusedHits.size() > 0){
    //First hit is automatically good.
    HitCalc* currentHit = new HitCalc(*unusedHits.back());
    currIndex++;
    locations[currentHit->GetID()] = currIndex;
    unusedHits.pop_back();
    //Give the new hit to the GretinaCalc to hold.
    gr->AddHitAB(currentHit);
    vector<int> IDsToCheck;
    IDsToCheck.push_back(currentHit->GetID());
    //Check each hit within the block for more hits adjacent to it.
    while(IDsToCheck.size() > 0){
      int IDchecking = IDsToCheck.back();
      int neighbor = -1;
      for (uint i=0; i<unusedHits.size(); i++){
	if (IsNeighbor(IDchecking, unusedHits[i]->GetID())){
	  neighbor = i;
	  break;
	}
      }
      if (neighbor == -1){
	//The current ID has no unfound neighbors.
	IDsToCheck.pop_back();
      } else {
	//A new neighbor to be added into the current hit.
	currentHit->AddBackHitCalc(unusedHits[neighbor]);
	IDsToCheck.push_back(unusedHits[neighbor]->GetID());
	locations[unusedHits[neighbor]->GetID()] = currIndex;
	unusedHits.erase(unusedHits.begin() + neighbor);
      }
    }
  }

}

//A simpler form of neighbor add-back.
//Instead of doing a search for all clusters of neighbors,
//it just checks if there are any neighbors already there.
//Will have stranger behavior for hits in 1-2-3 pattern,
//but won't add arbitrarily many hits together.
void Calibration::AddBackGretinaNeighborsNoSearch(GretinaCalc* gr){
  vector<HitCalc*> filtered = FilterOverflows(gr);
  map<int,int> locations;
  for(vector<HitCalc*>::iterator iter = filtered.begin(); iter!=filtered.end(); iter++){
    HitCalc* currentHit = *iter;
    bool found = false;
    for(int j=0; j<gr->GetMultAB(); j++){
      if(IsNeighbor(currentHit->GetID(),gr->GetHitAB(j)->GetID())){
	locations[currentHit->GetID()] = j;
	gr->GetHitAB(j)->AddBackHitCalc(currentHit);
	found = true;
	break;
      }
    }
    if(!found){
      locations[currentHit->GetID()] = gr->GetMultAB();
      gr->AddHitAB(new HitCalc(*currentHit));
    }
  }

}

 void Calibration::AddBackGretinaEverything(GretinaCalc* gr){
  //All hits within Gretina are summed.
  vector<HitCalc*> filtered = FilterOverflows(gr);
  if (filtered.size() < 1)
    return;
  for(vector<HitCalc*>::iterator iter = filtered.begin(); iter!=filtered.end(); iter++){
    HitCalc* hit = *iter;
    if(iter==filtered.begin()){
      gr->AddHitAB(new HitCalc(*hit));
    } else {
      gr->GetHitAB(0)->AddBackHitCalc(hit);
    }
  }

}

void Calibration::AllGretinaHits(GretinaCalc* gr, Gretina* in){
  //monstercluster for Ragnar
  vector<HitCalc*> cluster = ExtractAllHits(in);
  if (cluster.size() > 0 )
    gr->AddHitCL(cluster,0);
}
void Calibration::GammaTrack(GretinaCalc* gr, GretinaEvent* gt){
  ftracking->SetGretina(gr);
  ftracking->SortInClusters();    
  gt = ftracking->GetEvent();
}

void Calibration::ClusterGretina(GretinaCalc* gr, Gretina* in){
  vector<HitCalc*> unusedHits = ExtractAllHits(in);
  int c =0; //clustercounter
  //All hits within a cone should be summed.
  //cone adapts for each hit.
  while (unusedHits.size() > 0){
     //First hit is automatically good.
    vector<HitCalc*> curCluster;
    curCluster.push_back(unusedHits.back());

    //also add-back clusters
    HitCalc* curHit = new HitCalc(*unusedHits.back());
    gr->AddHitAB(curHit);

    unusedHits.pop_back();

    int added =0;
    do{
      added =0;
      int incone = -1;
      bool found = false;
      for(UShort_t k=0;k<curCluster.size();k++){
	for(UShort_t j=0;j<unusedHits.size();j++){
	  if(curCluster[k]->GetPosition().Angle(unusedHits[j]->GetPosition())<fSett->ClusterAngle()*TMath::Pi()/180.){
	    incone = j;
	    found = true;
	    break;
	  }
	}
	if(found)
	  break;
      }
      if(found&&incone>-1){
	added++;
	curCluster.push_back(unusedHits[incone]);
	curHit->AddBackHitCalc(unusedHits[incone]);
	unusedHits.erase(unusedHits.begin() + incone);
      }
      if(unusedHits.size()==0)
	break;
    }while(added>0);
    gr->AddHitCL(curCluster,c);
    c++;
  }
}
void Calibration::ResetCtrs(){
  fichctr = 0;
  fscintctr = 0;
  fobjctr = 0;
  fxfpctr = 0;
  ftofctr = 0;
  fpadctr[0] = 0;
  fpadctr[1] = 0;
  ftrackctr = 0;
  fppacctr[0] = 0;
  fppacctr[1] = 0;
  fiitrackctr = 0;
  fcard29ctr = 0;
  fgretactr = 0;
}
void Calibration::PrintCtrs(){
  cout << "fichctr    \t" << fichctr    << endl;
  cout << "fscintctr  \t" << fscintctr  << endl;
  cout << "fobjctr    \t" << fobjctr    << endl;
  cout << "fxfpctr    \t" << fxfpctr    << endl;
  cout << "ftofctr    \t" << ftofctr    << endl;
  cout << "fpadctr[0] \t" << fpadctr[0] << endl;
  cout << "fpadctr[1] \t" << fpadctr[1] << endl;
  cout << "ftrackctr  \t" << ftrackctr  << endl;
  cout << "fppacctr[0]\t" << fppacctr[0]<< endl;
  cout << "fppacctr[1]\t" << fppacctr[1]<< endl;
  cout << "fiitrackctr\t" << fiitrackctr<< endl;
  cout << "fcard29ctr \t" << fcard29ctr << endl;
  cout << "fgretactr  \t" << fgretactr  << endl;
}
