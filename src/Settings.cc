#include "Settings.hh"

using namespace std;
Settings::Settings(){
}

Settings::Settings(const char* filename){
  fInputFiles.push_back(filename);
  TEnv set(filename);
  ReadSettings(&set);
  if(fVerboseLevel>1)
    PrintSettings();
}

Settings::Settings(vector<char*> files){
  TEnv set;
  //Reverse order, because duplicate entries are ignored, instead of overwriting.
  for(vector<char*>::reverse_iterator it = files.rbegin(); it!=files.rend(); it++){
    set.ReadFile(*it,kEnvLocal);
    fInputFiles.push_back(*it);
  }
  ReadSettings(&set);
  if(fVerboseLevel>1)
    PrintSettings();
}

Settings::~Settings(){ }

void Settings::ReadSettings(TEnv* set){
  char* defaultfile = (char*)"~/analysis/settings/nocal.dat";

  fVerboseLevel = set->GetValue("VerboseLevel",0);
  fLastEvent = set->GetValue("LastEvent",-1);
  fEventTimeDiff = set->GetValue("EventTimeDiff", 500);
  fCrdcChannels = set->GetValue("Crdc.Channels",224);
  fSampleWidth = set->GetValue("Crdc.Sample.Width",12);
  fPpacChannels = set->GetValue("Ppac.Channels",256);
  fPpacWidth = set->GetValue("Ppac.Width",32);
  fPpacPitch = set->GetValue("Ppac.Pitch",1.27);
  fPpacGap = set->GetValue("Ppac.Gap",482.6);
  fPpacZ = set->GetValue("Ppac.Z",0.0);
  fIonChamberChannels = set->GetValue("IonChamber.Channels",16);
  if(fVerboseLevel>2)
    cout << "reading calfiles" << endl;
  fCalFile = set->GetValue("Crdc.File",defaultfile);
  fPedestalFile = set->GetValue("Crdc.Ped.File",defaultfile);
  for(int i=0;i<2;i++){
    fSaturation[i] = set->GetValue(Form("Saturation.%d",i),1000);
    fMSaturation[i] = set->GetValue(Form("Saturation.Multiplier.%d",i),1);
    fGravityWidth[i] = set->GetValue(Form("Crdc.Gravity.Width.%d",i),12);
    fFitWidth[i] = set->GetValue(Form("Crdc.Fit.Width.%d",i),20);
    fMethod[i] = set->GetValue(Form("Crdc.Method.%d",i),0);
    fxOffset[i] = set->GetValue(Form("Crdc.X.Offset.%d",i),0.0);
    fxSlope[i] = set->GetValue(Form("Crdc.X.Slope.%d",i),1.0);
    fyOffset[i] = set->GetValue(Form("Crdc.Y.Offset.%d",i),0.0);
    fyCoincHeadstart[i] = set->GetValue(Form("Crdc.Y.CoincHeadstart.%d",i),0.0);
    fySlope[i] = set->GetValue(Form("Crdc.Y.Slope.%d",i),1.0);
    fyCorOffset[i] = set->GetValue(Form("Crdc.Y.Correction.Offset.%d",i),0.0);
    fyCorSlope[i] = set->GetValue(Form("Crdc.Y.Correction.Slope.%d",i),0.0);
  }
  fBadFile = set->GetValue("BadPad.File",defaultfile);
  fMapFile = set->GetValue("Map.File",defaultfile);
  fTrackParams = set->GetValue("Track.Params",6);
  fTrackCoefs = set->GetValue("Track.Coefs",200);
  fGap = set->GetValue("Crdc.Gap",1073.);
  fFPShift = set->GetValue("Focal.Plane.Shift",0.);
  fShiftY = set->GetValue("Shift.y",0.0);
  fAngleA = set->GetValue("Angle.a",0.0);
  fAngleB = set->GetValue("Angle.b",0.0);
  fGretAngleA = set->GetValue("Gretina.Angle.a",0.0);
  fGretAngleB = set->GetValue("Gretina.Angle.b",0.0);


  //Read the timing corrections.
  //Use the modified names of parameters if they are found.
  //Otherwise, use the older, less legible names.
  frfE1C = set->GetValue("RF.afp.Corr",0.0);
  if (frfE1C==0.0)
    frfE1C = set->GetValue("RFe1.Corr",0.0);
  frfC = set->GetValue("RF.xfp.Corr",0.0);
  if (frfC==0.0)
    frfC = set->GetValue("RF.Corr",0.0);
  fobjE1C = set->GetValue("OBJ.afp.Corr",0.0);
  if (fobjE1C==0.0)
    fobjE1C = set->GetValue("OBJe1.Corr",0.0);
  fobjC = set->GetValue("OBJ.xfp.Corr",0.0);
  if (fobjC==0.0)
    fobjC = set->GetValue("OBJ.Corr",0.0);
  fxfpE1C = set->GetValue("XFP.afp.Corr",0.0);
  if (fxfpE1C==0.0)
    fxfpE1C = set->GetValue("XFPe1.Corr",0.0);
  fxfpC = set->GetValue("XFP.xfp.Corr",0.0);
  if (fxfpC==0.0)
    fxfpC = set->GetValue("XFP.Corr",0.0);
  ftac_objE1C = set->GetValue("TACOBJ.afp.Corr",0.0);
  if (ftac_objE1C==0.0)
    ftac_objE1C = set->GetValue("TACOBJe1.Corr",0.0);
  ftac_objC = set->GetValue("TACOBJ.xfp.Corr",0.0);
  if (ftac_objC==0.0)
    ftac_objC = set->GetValue("TACOBJ.Corr",0.0);
  ftac_xfpE1C = set->GetValue("TACXFP.afp.Corr",0.0);
  if (ftac_xfpE1C==0.0)
    ftac_xfpE1C = set->GetValue("TACXFPe1.Corr",0.0);
  ftac_xfpC = set->GetValue("TACXFP.xfp.Corr",0.0);
  if (ftac_xfpC==0.0)
    ftac_xfpC = set->GetValue("TACXFP.Corr",0.0);


  fdE_ytilt = set->GetValue("IC.Y.Corr",0.0);
  fdE_xtilt = set->GetValue("IC.X.Corr",0.0);
  fdE_x0tilt = set->GetValue("IC.X0.Corr",0.0);
  fCalFileIC = set->GetValue("IC.Cal.File",defaultfile);
  fICthresh = set->GetValue("IC.Thresh",0.0);
  fTimeCorFile = set->GetValue("Time.Corrections.File",defaultfile);
  fEvtNrFile = set->GetValue("Event.Number.File",defaultfile);

  ftargetX = set->GetValue("Target.X",0.0);
  ftargetY = set->GetValue("Target.Y",0.0);
  ftargetZ = set->GetValue("Target.Z",0.0);
  ftargetBeta = set->GetValue("Target.Beta",0.0);
  fAccepFile = set->GetValue("Acceptance.File",defaultfile);
  fAccepCutOff[0] = set->GetValue("Acceptance.Extrapol",-2.0);
  fAccepCutOff[1] = set->GetValue("Acceptance.CutOff",20.0);
  fDTACorr = set->GetValue("DTA.Corr",0.0);

  fAddBackType = set->GetValue("AddBackType",0);
  fClusterAngle = set->GetValue("ClusterAngle",20);
  fStoreAllIPoints = set->GetValue("StoreAllIPoints",0);

  fOverflowThreshold = set->GetValue("OverflowThreshold",6000);

  fMatrixFile = set->GetValue("Gretina.Matrix.File",defaultfile);
  fNeighborFile = set->GetValue("Gretina.Neighbor.File",defaultfile);
  fGretinaCalFile = set->GetValue("Gretina.Cal.File",defaultfile);
  fGretinaRecalFile = set->GetValue("Gretina.Recalibrate.File",defaultfile);

  fdet2hole.clear();
  fhole2det.clear();
  for(int det=0; det<7; det++){
    int hole = set->GetValue(Form("Detector.%d",det),-1);
    if (hole!=-1){
      fdet2hole[det] = hole;
      fhole2det[hole] = det;
    }
  }

  fHodoCalFile = set->GetValue("Hodo.Cal.File",defaultfile);
  fHodoThresh = set->GetValue("Hodo.Thresh",0.0);
  fScintThresh = set->GetValue("Scintillator.Thresh",0.0);

  fTracking = set->GetValue("DoTracking",false);

  fobj_range[0] = set->GetValue("OBJ.NBins",500);
  fobj_range[1] = set->GetValue("OBJ.Low",0);
  fobj_range[2] = set->GetValue("OBJ.High",4000);
  fxfp_range[0] = set->GetValue("XFP.NBins",500);
  fxfp_range[1] = set->GetValue("XFP.Low",0);
  fxfp_range[2] = set->GetValue("XFP.High",4000);
  ftacobj_range[0] = set->GetValue("TACOBJ.NBins",500);
  ftacobj_range[1] = set->GetValue("TACOBJ.Low",0);
  ftacobj_range[2] = set->GetValue("TACOBJ.High",4000);
  ftacxfp_range[0] = set->GetValue("TACXFP.NBins",500);
  ftacxfp_range[1] = set->GetValue("TACXFP.Low",0);
  ftacxfp_range[2] = set->GetValue("TACXFP.High",4000);
  fIC_range[0] = set->GetValue("IC.NBins",1000);
  fIC_range[1] = set->GetValue("IC.Low",0);
  fIC_range[2] = set->GetValue("IC.High",1000);
  fobjC_range[0] = set->GetValue("OBJ.Corr.NBins",500);
  fobjC_range[1] = set->GetValue("OBJ.Corr.Low",0);
  fobjC_range[2] = set->GetValue("OBJ.Corr.High",4000);
  fxfpC_range[0] = set->GetValue("XFP.Corr.NBins",500);
  fxfpC_range[1] = set->GetValue("XFP.Corr.Low",0);
  fxfpC_range[2] = set->GetValue("XFP.Corr.High",4000);
  ftacobjC_range[0] = set->GetValue("TACOBJ.Corr.NBins",500);
  ftacobjC_range[1] = set->GetValue("TACOBJ.Corr.Low",0);
  ftacobjC_range[2] = set->GetValue("TACOBJ.Corr.High",4000);
  ftacxfpC_range[0] = set->GetValue("TACXFP.Corr.NBins",500);
  ftacxfpC_range[1] = set->GetValue("TACXFP.Corr.Low",0);
  ftacxfpC_range[2] = set->GetValue("TACXFP.Corr.High",4000);
  fPP_range[0] = set->GetValue("PP.NBins",100.);
  fPP_range[1] = set->GetValue("PP.Low",5.0);
  fPP_range[2] = set->GetValue("PP.High",15.0);

}

int Settings::Det2Hole(int det){
  if (fdet2hole.count(det)==1){
    return fdet2hole[det];
  } else {
    return -1;
  }
}

int Settings::Hole2Det(int hole){
  if (fhole2det.count(hole)==1){
    return fhole2det[hole];
  } else {
    return -1;
  }
}

void Settings::PrintSettings(){
  cout << "LastEvent\t" << fLastEvent << endl;
  cout << "Crdc.Channels\t" << fCrdcChannels << endl;
  cout << "Crdc.Sample.Width\t" << fSampleWidth << endl;
  cout << "IonChamber.Channels\t" << fIonChamberChannels << endl;
  cout << "Crdc.File\t"<< fCalFile << endl;
  cout << "Crdc.Ped.File\t"<< fPedestalFile << endl;
  for(int i=0;i<2;i++){
    //cout << Form("Crdc.File.%d\t",i) << fCalFile[i] << endl;
    cout << Form("Saturation.%d\t",i) << fSaturation[i] << endl;
    cout << Form("Saturation.Multiplier.%d\t",i) << fMSaturation[i] << endl;
    cout << Form("Crdc.Gravity.Width.%d\t",i) << fGravityWidth[i] << endl;
    cout << Form("Crdc.Fit.Width.%d\t",i) << fFitWidth[i] << endl;
    cout << Form("Crdc.Method.%d\t",i) << fMethod[i] << endl;
    cout << Form("Crdc.X.Offset.%d\t",i) << fxOffset[i] << endl;
    cout << Form("Crdc.X.Slope.%d\t",i) << fxSlope[i] << endl;
    cout << Form("Crdc.Y.Offset.%d\t",i) << fyOffset[i] << endl;
    cout << Form("Crdc.Y.Slope.%d\t",i) << fySlope[i] << endl;
    cout << Form("Crdc.Y.Correction.Offset.%d\t",i) << fyCorOffset[i] << endl;
    cout << Form("Crdc.Y.Correction.Slope.%d\t",i) << fyCorSlope[i] << endl;
  }
  cout << "Crdc.Gap\t" << fGap << endl;
  cout << "Focal.Plane.Shift\t" << fFPShift << endl;
  cout << "BadPad.File\t" << fBadFile << endl;
  cout << "Ppac.Channels\t" << fPpacChannels << endl;
  cout << "Ppac.Width\t" << fPpacWidth << endl;
  cout << "Ppac.Pitch\t" << fPpacPitch << endl;
  cout << "Ppac.Gap\t" << fPpacGap << endl;
  cout << "Ppac.Z\t" << fPpacZ << endl;
  cout << "Map.File\t" << fMapFile << endl;
  cout << "Track.Params\t" << fTrackParams << endl;
  cout << "Track.Coefs\t" << fTrackCoefs << endl;
  cout << "Shift.y\t" << fShiftY << endl;
  cout << "Angle.a\t" << fAngleA << endl;
  cout << "Angle.b\t" << fAngleB << endl;
  cout << "Gretina.Angle.a\t" << fGretAngleA << endl;
  cout << "Gretina.Angle.b\t" << fGretAngleB << endl;
  cout << "RFe1.Corr\t" << frfE1C << endl;
  cout << "RF.Corr\t" << frfC << endl;
  cout << "OBJe1.Corr\t" << fobjE1C << endl;
  cout << "OBJ.Corr\t" << fobjC << endl;
  cout << "XFPe1.Corr\t" << fxfpE1C << endl;
  cout << "XFP.Corr\t" << fxfpC << endl;
  cout << "TACOBJe1.Corr\t" << ftac_objE1C << endl;
  cout << "TACOBJ.Corr\t" << ftac_objC << endl;
  cout << "TACXFPe1.Corr\t" << ftac_xfpE1C << endl;
  cout << "TACXFP.Corr\t" << ftac_xfpC << endl;
  cout << "IC.Y.Corr\t" << fdE_ytilt << endl;
  cout << "IC.X.Corr\t" << fdE_xtilt << endl;
  cout << "IC.X0.Corr\t" << fdE_x0tilt << endl;
  cout << "IC.Cal.File\t" << fCalFileIC << endl;
  cout << "IC.Thresh\t" << fICthresh << endl;
  cout << "Time.Corrections.File\t" << fTimeCorFile << endl;
  cout << "Event.Number.File\t" << fEvtNrFile << endl;

  cout << "Target.X\t" << ftargetX << endl;
  cout << "Target.Y\t" << ftargetY << endl;
  cout << "Target.Z\t" << ftargetZ << endl;
  cout << "Target.Beta\t" << ftargetBeta << endl;

  cout << "Gretina.Matrix.File\t" << fMatrixFile << endl;
  cout << "Gretina.Neighbor.File\t" << fNeighborFile << endl;
  cout << "Gretina.Cal.File\t" << fGretinaCalFile << endl;
  cout << "Gretina.Recalibrate.File\t" << fGretinaRecalFile << endl;

  cout << "Hodo.Cal.File\t" << fHodoCalFile << endl;
  cout << "Hodo.Thresh\t" << fHodoThresh << endl;

  cout << "AddBackType\t" << fAddBackType << endl;
  cout << "ClusterAngle\t" << fClusterAngle << endl;
  cout << "StoreAllIPoints\t" << fStoreAllIPoints << endl;
  cout << "OverflowThreshold\t"<< fOverflowThreshold << endl;

  cout << "DoTracking\t"<< fTracking << endl;

  //cout << "\t" << f << endl;
  //cout << "\t" << f << endl;
  //cout << "\t" << f << endl;

}
