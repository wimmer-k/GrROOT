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
#include "UnpackedEvent.hh"

using namespace std;
using namespace TMath;
double rad2deg = 180./TMath::Pi();
double deg2rad = TMath::Pi()/180.;

void swapbytes(char* a, char *b){
  char tmp=*a;
  *a=*b;
  *b=tmp;
}

// Mode 3 data is high endian format
void  HEtoLE(char* cBuf, int bytes){
  for(int i=0; i<bytes; i+=2)
    swapbytes((cBuf+i), (cBuf+i+1));
}

void bits_uint(unsigned int value){
  unsigned int bit;
  for( bit = (~0U >> 1) + 1; bit > 0; bit >>= 1 ){
    putchar(value & bit ? '1' : '0');
  }
  putchar('\n');
}

UnpackedEvent::UnpackedEvent(Settings* settings = NULL){
  if(settings!=NULL){
    fvl = settings->VLevel();
  } else {
    fvl = 0;
  }

  fSett = settings;
}

void UnpackedEvent::Init(){
  frhist = new RawHistograms(fSett);
  fchist = new CalHistograms(fSett);

  TEnv *cal = new TEnv(fSett->GretinaCalFile());
  for(int d=0;d<7;d++){
    for(int c=0;c<4;c++){
      fslope[d][c] = cal->GetValue(Form("Slope.d%d.c%d",d,c),1.0);
      foffset[d][c] = cal->GetValue(Form("Offset.d%d.c%d",d,c),1.0);
    }
  }
  fEventTimeDiff = fSett->EventTimeDiff();

  fGretina = new Gretina;
  fS800 = new S800;
  fMode3Event = new Mode3Event;
  fScaler = new Scaler;

  fS800Calc = new S800Calc;
  fMode3Calc = new Mode3Calc;
  fGretinaCalc = new GretinaCalc;
  if(fwtree){
    //setting up tree
    cout << "UnpackedEvent: " << "setting up raw tree " << endl;
    ftr = new TTree("gtr","Gretina/S800 built events");
    ftr->Branch("gretina",&fGretina, 320000);
    ftr->Branch("s800",&fS800, 320000);
    ftr->Branch("mode3Event",&fMode3Event, 320000);
    ftr->BranchRef();
    cout << "UnpackedEvent: " << "done setting up raw tree" << endl;
  }
  if(fwcaltree){
    //setting up tree
    cout << "UnpackedEvent: " << "setting up calibrated tree " << endl;
    fcaltr = new TTree("ctr","Gretina/S800 calibrated and builtevents");
    fcaltr->Branch("s800calc",&fS800Calc, 320000);
    fcaltr->Branch("gretinacalc",&fGretinaCalc, 320000);
    fcaltr->Branch("mode3calc",&fMode3Calc, 320000);
    fcaltr->BranchRef();
    cout << "UnpackedEvent: " << "done setting up calibrated tree" << endl;
  }
  fnentries = 0;
  fhits = 0;
  fstrangehits = 0;
  fGretina->Clear();
  fS800->Clear();
  fScaler->Clear();
  fMode3Event->Clear();
  fhasdata = false;
  fhasscaler = false;
  fcurrent_ts = 0;
  ffirst_ts = 0;
  fctr = 0;
  frecalibrate = 0;
  for(int i=0;i<NSCALER;i++){
    fprevscaler[i] = 0;
    foverflows[i] = 0;
  }

  fncalentries = 0;
  fS800Calc->Clear();
  fMode3Calc->Clear();

}

int UnpackedEvent::DecodeGretina(Crystal* cryst, long long int gts){

  if(ffirst_ts<1){
    if(fvl>1)
      cout << "UnpackedEvent: " << "setting first timestamp " << gts << endl;
    ffirst_ts = gts;
  }
  if(fvl>1){
    cout << "UnpackedEvent: " << "-----------------------------next hit: "<< fnentries<< endl;
    cryst->PrintEvent();
  }
  int det = cryst->GetHole();
  int cry = cryst->GetCrystal();

  if(det<0 || det > MAXDETPOS){ // i don't know whether det runs from 0 or 1
    cout << "UnpackedEvent: " << "invalid detector number " << det << endl;
    return 11;
  }
  if(cry<0 || cry > MAXCRYSTALNO-1){
    cout << "UnpackedEvent: " << "invalid crystal number " << cry << endl;
    return 12;
  }



  if(fvl>1){
    cout << "UnpackedEvent: " << "mult " << cryst->GetMult() <<"\ten " <<  cryst->GetEnergy() <<"\tts " <<  cryst->GetTS() <<"\ten max " <<  cryst->GetMaxEn()<<"\tip max " <<  cryst->GetMaxIPNr()<< endl;
    cout << "UnpackedEvent: " <<"-----------------------------"<< endl;
  }

  fhits++;

  //che events which have no good interaction points
  if(cryst->GetMaxIPNr()<0){
    fstrangehits++;
  }
  if(fvl>0 && cryst->GetMaxIPNr()<0){
    cryst->PrintEvent();;
  }
  // now check time stamps
  long long int deltaEvent = gts - fcurrent_ts;
  if(fcurrent_ts>0 && deltaEvent < 0 )
    cout << "UnpackedEvent: " << "Inconsistent Timestamp last time was " << fcurrent_ts << " this (gretina) " << gts << " difference " << deltaEvent<< endl;

  if(fvl>1){
    cout << "UnpackedEvent: " <<fnentries<< "this ts " << gts <<" current ts " << fcurrent_ts <<" difference " << deltaEvent <<endl;
  }

  if(deltaEvent  < fEventTimeDiff){
    if(fvl>1)
      cout << "UnpackedEvent: " <<fnentries<< " coincidence difference " << deltaEvent << endl;
  }
  else{
    if(fvl>1)
      cout << "UnpackedEvent: " <<fnentries << " gretina single time difference " << deltaEvent << endl;
    if(fcurrent_ts>0){
      if(fvl>2)
	cout << "UnpackedEvent: " << "Closing event due to timestamp in Gretina." << endl;
      if(fvl>1&&fMode3Event->GetMult()<1)
	cout << "UnpackedEvent: " <<__PRETTY_FUNCTION__<< " entry " << fnentries << " mode3 is empty " << endl;
      if(fvl>1&&fhasdata==false)
	cout << "UnpackedEvent: " <<__PRETTY_FUNCTION__<< " entry " << fnentries << " s800 is empty " << endl;
      fMode3Event->SetCounter(fctr);
      fctr = 0;
      this->CloseEvent();
    }
    this->ClearEvent();
  }




  fGretina->AddHit(cryst);
  //set the current timestamp
  fcurrent_ts = gts;

  if(fvl>2){
    cout << "UnpackedEvent: " << "Gretina event found with timestamp " << gts << endl;
  }
  return 0;
}

int UnpackedEvent::DecodeS800(unsigned short *pevent, long long int ts, unsigned short twords){
  // now check time stamps
  long long int deltaEvent = ts - fcurrent_ts;
  if(fcurrent_ts>0 && deltaEvent < 0 )
    cout << "UnpackedEvent: " << "Inconsistent Timestamp last time was " << fcurrent_ts << " this (s800) " << ts << " difference " << deltaEvent<< endl;



  if(fvl>1){
    cout << "UnpackedEvent: " <<fnentries<< " this ts " << ts <<" current ts " << fcurrent_ts <<" difference " << deltaEvent <<endl;
  }

  if(deltaEvent  < fEventTimeDiff){
    if(fvl>1)
      cout << "UnpackedEvent: " <<fnentries<< " coincidence difference " << deltaEvent << endl;
    if(fhasdata==true){
      if(fvl>2)
	cout << "UnpackedEvent: " << "Closing event due to two S800 entries." << endl;
      cout << "UnpackedEvent: " << " coincidence with another S800! deltaT = " << deltaEvent << " writing and clearing last event" << endl;
      if(fvl>1&&fMode3Event->GetMult()<1)
	cout << "UnpackedEvent: " <<__PRETTY_FUNCTION__<< " entry " << fnentries << " mode3 is empty " << endl;
      if(fvl>1&&fhasdata==false)
	cout << "UnpackedEvent: " <<__PRETTY_FUNCTION__<< " entry " << fnentries << " s800 is empty " << endl;
      fMode3Event->SetCounter(fctr);
      fctr = 0;
      this->CloseEvent();
    }
  } else {
    if(fvl>1)
      cout << "UnpackedEvent: " <<fnentries << " s800 single event difference " << deltaEvent << endl;
    if(fcurrent_ts>0&&fhasdata==true){
      if(fvl>2)
	cout << "UnpackedEvent: " << "Closing event due to timestamp in S800." << endl;
      if(fvl>1&&fMode3Event->GetMult()<1)
	cout << "UnpackedEvent: " <<__PRETTY_FUNCTION__<< " entry " << fnentries << " mode3 is empty " << endl;
      if(fvl>1&&fhasdata==false)
	cout << "UnpackedEvent: " <<__PRETTY_FUNCTION__<< " entry " << fnentries << " s800 is empty " << endl;
      fMode3Event->SetCounter(fctr);
      fctr = 0;
      this->CloseEvent();
    }
    this->ClearEvent();
  }
  //set the current timestamp
  fcurrent_ts = ts;


  long long int n;
  unsigned short *p=pevent, sublength, plength, ptag, ID;
  unsigned short nwords, words;
  ++p;
  bool found=false;

  twords--; // total event length is self-inclusive
  if(fvl>1){
    cout << "UnpackedEvent: " << "S800 twords " << twords << endl;;
    cout << "UnpackedEvent: " << "next entry " << endl;
  }
  // Loop over sub-events until the right one is found
  while(twords > 0){
    nwords = *p; ++p;
    sublength = nwords;
    nwords--;
    if(*p++ != S800_PACKET){
      twords -= sublength;
      p += sublength - 2;
    }
    else{
      nwords--;
      found = true;
      fhasdata = true;
      break;
    }
  }
  if(!found)
    return 0;


  // Unpack S800 data in tree
  if(*p++ != S800_VERSION){
    cout << "UnpackedEvent: " << "Wrong version of S800 sub-event. Aborting ..." << endl;
    return 1;
  }
  nwords--;
  while(nwords > 0){
    plength = *p; ++p;
    ptag = *p; ++p;
    switch(ptag){
    case S800_TRIGGER_PACKET:
      if(fvl>1)
	cout << "UnpackedEvent: " << "trigger " << endl;
      p = DecodeS800Trigger(p);
      break;

    case S800_TOF_PACKET:
      if(fvl>1)
	cout << "UnpackedEvent: " << "tof " << endl;
      p = DecodeS800TimeOfFlight(p);
      break;

    case S800_FP_SCINT_PACKET:
      if(fvl>1)
	cout << "UnpackedEvent: " << "FP scint " << endl;
      words = plength-2;
      while (words > 0) {
	ID = ((*p)&0xF000)>>12;
	if(fvl>1)
	  cout << "UnpackedEvent: " << "ID " << ID << " ID/2 " << ID/2 << endl;
	p = DecodeS800Scintillator(p, ID, ID/2);
	words -= 2;
      }
      break;

    case S800_FP_IC_PACKET:
      if(fvl>1)
      cout << "UnpackedEvent: " << "FP IC " << endl;
      p = DecodeS800IonChamber(p);
      break;


    case S800_FP_CRDC_PACKET:
      ID = *p; ++p;
      if(fvl>1)
	cout << "UnpackedEvent: " << "FP CRDC " << ID<< endl;
      p = DecodeS800Crdc(p, ID, true);
      break;

    case S800_II_CRDC_PACKET:
      ID = *p; ++p;
      //if(fvl>1)
	cout << "UnpackedEvent: " << "II CRDC " << ID<< endl;
      //p = DecodeS800Crdc(p, ID, false);
      p += plength - 2;
      break;

    case S800_II_TRACK_PACKET:
      p = DecodeS800Ppac(p);
      break;

    case S800_TIMESTAMP_PACKET:
      if(fvl>1)
	cout << "UnpackedEvent: " << "Timestamp " << endl;
      n = *p++;
      n = (*p++<<16|n);
      n = (*p++<<16|n);
      n = (*p++<<16|n);
      fS800->SetInternalTS(n);
      //p += plength - 2;
      break;

    case S800_EVENT_NUMBER_PACKET:
      if(fvl>1)
	cout << "UnpackedEvent: " << "Event nr " << endl;
      n = *p++;
      n = (*p++<<16|n);
      n = (*p++<<16|n);
      fS800->SetEvtNr(n);
      //p += plength - 2;
      break;

    case S800_FP_HODO_PACKET:
      if(fvl>1)
	cout << "UnpackedEvent: " << "Hodoscope " << endl;
      if(plength-2>0)
	p = DecodeS800Hodoscope(p,plength-2);
      break;

    case S800_LABR_PACKET:
      if(fvl>1)
	cout << "UnpackedEvent: " << "LaBr" << endl;
      p = DecodeS800LaBr(p);
      break;

    case S800_OB_PIN_PACKET:
      if(fvl>1)
	cout << "UnpackedEvent: " << "OB PIN " << endl;
      if(plength-2 >0){
	p = DecodeS800ObPin(p);
      }
      break;

    case S800_GALOTTE_PACKET:
      if(fvl>1)
	cout << "UnpackedEvent: " << "Galotte " << endl;
      if(plength-2>0)
	p = DecodeS800Galotte(p,plength-2);
      break;

    default:
      //if(fvl>0)
      cout << "UnpackedEvent: " << "unidentified s800 packet: " << (hex) << ptag <<(dec)<< endl;
      p += plength - 2;
      break;
    }
    nwords -= plength;
  }
  fS800->SetTS(ts);
  if(ffirst_ts<1){
    if(fvl>1)
      cout << "UnpackedEvent: " << "setting first timestamp " << ts << endl;
    ffirst_ts = ts;
  }


  if(fvl>2){
    cout << "UnpackedEvent: " << "S800 event found with timestamp " << ts <<  endl;
  }

  return 0;
}

unsigned short* UnpackedEvent::DecodeS800LaBr(unsigned short *p){
  UShort_t words = *(p-2)-2;
  while(words>0){
    fS800->GetLaBr()->SetADC(*p>>12,        //First 4 bits are the detector number.
			     *p & 0x0FFF);  //Remaining 12 bits are the ADC value.
    p++;
    fS800->GetLaBr()->SetTDC(*p>>12,        //First 4 bits are the detector number (again).
			     *p & 0x0FFF);  //Remaining 12 bites are the TDC value.
    p++;
    words -= 2;
  }
  return p;
}

unsigned short* UnpackedEvent::DecodeS800ObPin(unsigned short *p){
  fS800->SetObPin(*p);
  p++;
  return p;
}

unsigned short* UnpackedEvent::DecodeS800TimeOfFlight(unsigned short *p){
  UShort_t words = (*(p-2))-2, ch, dum;
  Short_t rf = -1;
  Short_t obj = -1;
  Short_t xfp = -1;
  Short_t tar = -1;
  Short_t tac_obj = -1;
  Short_t tac_xfp = -1;
  while (words > 0) {
    ch = ((*p)&0xf000)>>12;
    int tmp = *p; ++p;
    if(ch == 12)
      rf  = (tmp)&0xfff;
    else if(ch == 13)
      obj = (tmp)&0xfff;
    else if(ch == 14)
      xfp = (tmp)&0xfff;
    else if(ch == 15)
      tar = (tmp)&0xfff;
#ifdef CS800_LINK_TOFTAC
    else if(ch==5)
      tac_obj = (tmp)&0xfff;
    else if(ch==4)
      tac_xfp = (tmp)&0xfff;
#endif
    else if(ch>0 && ch<8){
      dum = (tmp)&0xfff;
      if(fvl>1){
	cout << "timeofflight dummy variable found, ch " << ch << " value " << dum << endl;
      }
    }
    /*
    if (ch ==  6) dum = (tmp)&0xfff;
    if (ch ==  7) dum = (tmp)&0xfff;
    if (ch ==  1) dum = (tmp)&0xfff;
    if (ch ==  2) dum = (tmp)&0xfff;
    if (ch ==  3) dum = (tmp)&0xfff;
    if (ch ==  4) dum = (tmp)&0xfff;
    if (ch ==  5) dum = (tmp)&0xfff;
    */
    words--;
  }
  //cout << "UnpackedEvent: " << tac_obj <<"\t" << tac_xfp << endl;
  //cout << "UnpackedEvent: " << "("<<rf<<","<<obj<<","<<xfp<<","<<tar<<") " << endl;
  fS800->GetTimeOfFlight()->Set(rf, obj, xfp, tar);
  fS800->GetTimeOfFlight()->SetTAC(tac_obj, tac_xfp);
  //cout << "UnpackedEvent: " << fS800->GetTimeOfFlight()->GetRF() << endl;
  return p;
}
unsigned short* UnpackedEvent::DecodeS800Hodoscope(unsigned short *p, unsigned short length){
  //cout << __PRETTY_FUNCTION__ << endl;
  //cout << "length " << length;
  int id = *p++;
  //cout << "\tid " << id << endl;
  length--;
  // id =0,1 energies  next numbers are0xceee with c channel and e energy
  // id 0 is channel 0-15, id =1 channel 16-31
  if(id<2){
    for(int k=0;k<length;k++){
      int ch = ((*p & 0xf000) >> 12);
      int data = (*p & 0x0fff);
      //cout << (hex) << *p << "\t" << (dec) << *p  << "\t" << ch << "\t" << data << endl;
      fS800->GetHodoscope()->Set(data,ch+(16*id));
      *p++;
    }
  }
  // id = 2 hitpattern and time
  // 0x000i 0xaaaa 0xbbbb 0xtttt
  // where i is the id, a is the coincidence register A for the first 16 channels, b the register for the second 16 channels and t the TAC time.
  else{
    if(length!=3){
      cout << "inconsistent packet length for in " <<  __PRETTY_FUNCTION__ << endl;
      cout << "length is " << length << " should be 3 " << endl;
    }
    fS800->GetHodoscope()->SetHitpattern(*p,*(p+1));
    fS800->GetHodoscope()->SetTime(*(p+2));
    //cout << *p << "\t" << *(p+1) << "\t" << *(p+2) << endl;
    p+=3;
  }
  return p;
}
unsigned short* UnpackedEvent::DecodeS800Galotte(unsigned short *p, unsigned short length){
  //cout << "length " << length << endl;
  if(length<1 || length >3){
    cout << "inconsistent packet length for in " <<  __PRETTY_FUNCTION__ << endl;
    cout << "length is " << length << " should be < 4 " << endl;
  }
  for(int k=0;k<length;k++){
    int ch = ((*p & 0xf000) >> 12);
    int data = (*p & 0x0fff);
    //cout << (hex) << *p << "\t" << (dec) << *p  << "\t" << ch << "\t" << data << endl;
    fS800->GetGalotte()->Set(data,ch);
    *p++;
  }
  return p;
}
unsigned short* UnpackedEvent::DecodeS800Trigger(unsigned short *p){
  UShort_t words = (*(p-2))-2, ch;
  int registr = -1;
  int s800 = -1;
  int external1 = -1;
  int external2 = -1;
  int secondary = -1;
  registr = *p++;
  words--;
  while(words > 0){
    ch = ((*p)&0xf000)>>12;
    if (ch == 8) s800 = (*p++)&0xfff;
    if (ch == 9) external1 = (*p++)&0xfff;
    if (ch == 10) external2 = (*p++)&0xfff;
    if (ch == 11) secondary = (*p++)&0xfff;
    words--;
  }
  fS800->GetTrigger()->Set(registr, s800, external1, external2, secondary);
  return p;
}

unsigned short* UnpackedEvent::DecodeS800Scintillator(unsigned short *p, unsigned short updown, int id){
  int de_up = -1;
  int time_up = -1;
  int de_down = -1;
  int time_down = -1;

// Even updown: UP channels.  Odd ID: DOWN channels
  if(updown%2==0){
    de_up = (*p++)&0xfff;
    time_up = (*p++)&0xfff;
  }
  else{
    de_down = (*p++)&0xfff;
    time_down = (*p++)&0xfff;
  }
  //cout << "UnpackedEvent: " << de_up <<" "<< time_up <<" "<< de_down <<" "<< time_down << endl;
  fS800->GetScintillator(id)->SetID(id);
  fS800->GetScintillator(id)->Set(de_up, time_up, de_down, time_down);
  return p;
}

unsigned short* UnpackedEvent::DecodeS800IonChamber(unsigned short *p){
  UShort_t ch=-1;
  UShort_t raw=-1;
  if(*(p+1) == S800_FP_IC_ENERGY_PACKET){
    // IC packet with times
    UShort_t length;
    length = *p++;
    p++;
    length -= 2;
    while (length > 0) {
      //cout << "UnpackedEvent: " << "length " << length << " (*p) " << hex << (*p) <<dec;
      ch = ((*p)&0xf000)>>12;
      raw = (*p++)&0xfff;
      //cout << "UnpackedEvent: " << " ch " << ch << " raw " << raw << endl;
      length--;
      fS800->GetIonChamber()->Set(ch,raw);
    }
  }
  else{
    // Old style IC packet
    UShort_t words = (*(p-2))-2;
    while(words > 0){
      ch = ((*p)&0xf000)>>12;
      raw = (*p++)&0xfff;
      words--;
      fS800->GetIonChamber()->Set(ch,raw);
    }
  }
  //cout << "UnpackedEvent: " << "ch " << ch << " raw " << raw << endl;

  return p;
}

unsigned short* UnpackedEvent::DecodeS800Crdc(unsigned short *p, int id, bool fp){
  UShort_t anode=-1;
  UShort_t tac=-1;
  if(fp)
    fS800->GetCrdc(id)->SetID(id);
  else
    fS800->GetTCrdc(id)->SetID(id);

  Int_t tag;
  if(fp)
    tag = S800_FP_CRDC_PACKET;
  else
    tag = S800_II_CRDC_PACKET;
  if(*(p+1) == tag+1){
    p = DecodeS800CrdcRaw(p,id, fp);
  }
  if(*(p+1) == tag+5){
    anode = *(p+2);
    tac = *(p+3);
    p += 4;
  }
  if(fp)
    fS800->GetCrdc(id)->SetAnodeTAC(anode, tac);
  else
    fS800->GetTCrdc(id)->SetAnodeTAC(anode, tac);

  return p;
}

unsigned short* UnpackedEvent::DecodeS800CrdcRaw(unsigned short *p, int id, bool fp){
  static ULong_t total=0, failed=0;
  Short_t sampleBegin = 0, sampleWidth, isample, ichannel, cdata[4], connector, previous_sample = 0, ch;
  Short_t maxwidth = S800_CRDC_MAXWIDTH;
  Short_t channels;
  if(fp)
    channels = S800_FP_CRDC_CHANNELS;
  else
    channels = S800_II_CRDC_CHANNELS;

  unsigned short *pStore = p;
  bool mes1=true, mes2=true, mes3=true, mes4=true;
  bool debug = S800_DEBUG;
  UShort_t length = *p++;
  short i = length-3;
  p++;	// skip packet id
  UShort_t threshold = *p++;
  while(i > 0){
    if ((*p)>>15 != 1) {
      cout << "UnpackedEvent: " << "CRDC data is corrupted!" << endl;
      p++; i--;
      continue;
    }
    else{
      isample = ((*p)&0x7FC0)>>6;
      ichannel = (*p)&0x003F;
      if(i == length-3){
	sampleBegin = isample;
	previous_sample = isample;
      }
    }
    p++; i--;
    memset(cdata, 0, sizeof(cdata));
    while ((*p)>>15 == 0) {
      connector = ((*p)&0x0C00)>>10;
      cdata[connector] = (*p)&0x3FF;
      p++; i--;
      if (i == 0) break;
    }
    if(isample < sampleBegin || isample > sampleBegin+maxwidth){
      if(debug)
	printf("Warning in Crdc Unpack: inconsistent sample number: %d (first: %d)\n", isample, sampleBegin);
      mes1 = false;
      continue;
    }
    if(isample < previous_sample){
      if(debug)
	printf("Warning in Crdc Unpack: sample number lower than previous: %d (previous: %d)\n", isample, previous_sample);
      mes2 = false;
      continue;
    }
    previous_sample = isample;
    for(int j=0; j<4; j++){
      ch = ichannel+j*64;
      if(cdata[j] != 0 && ch < channels){
	if(cdata[j] < threshold){
	  if(debug)
	    printf("Warning in Crdc Unpack: data lower than threshold: %d (threshold: %d)\n", cdata[j], threshold);
	  mes3 = false;
	}
	else{
	  if(fp)
	    fS800->GetCrdc(id)->Set(cdata[j], isample, ch);
	  else
	    fS800->GetTCrdc(id)->Set(cdata[j], isample, ch);
	}
      } //cdata[j] != 0 && ch < channels
      else if(cdata[j] != 0 && ch >= channels){
	if(debug)
	  printf("Warning in Crdc Unpack: channel greater than limit: %d (limit: %d)\n", ch, channels);
	mes4 = false;
      }
    }
    sampleWidth = isample - sampleBegin + 1;
    fS800->GetCrdc(id)->SetSampleWidth(sampleWidth);

  }// i>0
  if (!mes1 || !mes2 || !mes3 || !mes4) failed++;
  total++;
  if (failed == 1000) {
    if (debug)
      printf ("Errors in Crdc Unpackings: %g%%\n", 1.0*failed/total*100);
    total = 0;
    failed = 0;
  }
  //cout << "UnpackedEvent: " << "sampleWidth = "<<isample<<" - "<<sampleBegin<<" + 1 "<< sampleWidth  << endl;
  return (pStore+length);
}

unsigned short* UnpackedEvent::DecodeS800Ppac(unsigned short *p){
  if(*(p+1) == S800_II_TRACK_RAW_PACKET)
    p = DecodeS800PpacRaw(p);
  else
    cout << "UnpackedEvent: " << " Tracking PPAC Wrong packet " << endl;
  return p;
}

unsigned short* UnpackedEvent::DecodeS800PpacRaw(unsigned short *p){
  static ULong_t total=0, failed=0;
  Short_t sampleBegin = 0, sampleWidth, isample, ichannel, cdata[4], connector, previous_sample = 0, ch;
  Short_t maxwidth = S800_II_TRACK_MAX_WIDTH;
  Short_t channels = S800_II_TRACK_CHANNELS;
  unsigned short *pStore = p;
  bool mes1=true, mes2=true, mes3=true, mes4=true;
  bool debug = S800_DEBUG;
  UShort_t length = *p;p++;
  short i = length-3;
  p++;	// skip packet id
  UShort_t threshold = *p;p++;
  while(i > 0){
    if ((*p)>>15 != 1) {
      cout << "UnpackedEvent: " << "TPPAC data is corrupted!" << endl;
      p++; i--;
      continue;
    }
    else{
      isample = ((*p)&0x7FC0)>>6;
      ichannel = (*p)&0x003F;
      if(i == length-3){
	sampleBegin = isample;
	previous_sample = isample;
      }
    }
    p++; i--;
    memset(cdata, 0, sizeof(cdata));
    while ((*p)>>15 == 0) {
      connector = ((*p)&0x0C00)>>10;
      //cout << "UnpackedEvent: " << "id " << id << " connector " << connector << " isample " << isample << " ichannel " << ichannel<< endl;
      cdata[connector] = (*p)&0x03FF;
      p++; i--;
      if (i == 0) break;
    }
    if(isample < sampleBegin || isample > sampleBegin+maxwidth){
      if(debug)
	printf("Warning in Tppac Unpack: inconsistent sample number: %d (first: %d)\n", isample, sampleBegin);
      mes1 = false;
      continue;
    }
    if(isample < previous_sample){
      if(debug)
	printf("Warning in Tppac Unpack: sample number lower than previous: %d (previous: %d)\n", isample, previous_sample);
      mes2 = false;
      continue;
    }
    previous_sample = isample;
    for(int j=0; j<4; j++){
      ch = ichannel+j*64;
      if(cdata[j] != 0 && ch < channels){
	if(cdata[j] < threshold){
	  if(debug)
	    printf("Warning in Crdc Unpack: data lower than threshold: %d (threshold: %d)\n", cdata[j], threshold);
	  mes3 = false;
	}
	else{
	  //cout << "UnpackedEvent: " << "id " << id << " j " << j << " cdata " << cdata[j] << " isample " << isample << " ch " << ch << endl;
	  fS800->GetTPpac()->Set(cdata[j], isample, ch);
	}
      } //cdata[j] != 0 && ch < channels
      else if(cdata[j] != 0 && ch >= channels){
	if(debug)
	  printf("Warning in Crdc Unpack: channel greater than limit: %d (limit: %d)\n", ch, channels);
	mes4 = false;
      }
    }
    sampleWidth = isample - sampleBegin + 1;
    fS800->GetTPpac()->SetSampleWidth(sampleWidth);

  }// i>0
  if (!mes1 || !mes2 || !mes3 || !mes4) failed++;
  total++;
  if (failed == 1000) {
    if (debug)
      printf ("Errors in TPpac Unpackings: %g%%\n", 1.0*failed/total*100);
    total = 0;
    failed = 0;
  }
  return (pStore+length);
}

int UnpackedEvent::DecodeScaler(unsigned short *pevent, long long int ts){
  UInt_t nScaler;
  ULong64_t upper,lower;

  //Reading the header of the scaler item.
  //Note that this is different from the GEB header.
  //Scaler packets have their own internal headers.
  lower   = (ULong64_t)*pevent++;
  upper   = (ULong64_t)*pevent++;
  ULong64_t size = (ULong64_t) (upper<<16|lower);
  lower   = (ULong64_t)*pevent++;
  upper   = (ULong64_t)*pevent++;
  ULong64_t type = (ULong64_t) (upper<<16|lower);

  ULong64_t timestamp = ts;

  if (type == NONEVENT_ITEM_TYPE){
    lower   = (ULong64_t)*pevent++;
    upper   = (ULong64_t)*pevent++;
    ULong64_t timestamp1 = (ULong64_t) (upper<<16|lower);
    lower   = (ULong64_t)*pevent++;
    upper   = (ULong64_t)*pevent++;
    ULong64_t timestamp2 = (ULong64_t) (upper<<16|lower);

    timestamp =  timestamp2 << 32;
    timestamp += timestamp1;
  }
  if(fvl>3)
    cout << "size " << size << "\ttimestamp " << timestamp << "\ttype" << type  <<endl;


  if(type == PHYSICS_EVENT_COUNT_TYPE){
    return 1;
  }
  else if(type == SCALER_ITEM_TYPE || type == NONEVENT_ITEM_TYPE){

    fhasscaler = true;
    fScaler->SetTS(timestamp); //for old data this is the header timestamp, for the newer data (>jan 2013) this would be the internal TS
    fScaler->SetEventNr(max(fnentries,fncalentries));

    //First 8 bytes are the start in little-endian form.
    lower   = (ULong64_t)*pevent++;
    upper   = (ULong64_t)*pevent++;
    fScaler->SetStart((ULong64_t)(upper<<16|lower));

    //Next 8 bytes are the end
    lower   = (ULong64_t)*pevent++;
    upper   = (ULong64_t)*pevent++;
    fScaler->SetEnd((ULong64_t)(upper<<16|lower));

    if(fvl>1)
      cout << "UnpackedEvent, DecodeScaler: start " << fScaler->GetStart() << " end " << fScaler->GetEnd() << endl;

    if(type == NONEVENT_ITEM_TYPE){
      //intervalDivisor.
      lower   = (ULong64_t)*pevent++;
      upper   = (ULong64_t)*pevent++;
      if(fvl>1)
	cout << "interval divisor " << (ULong64_t)(upper<<16|lower) << endl;
    }

    //Next 8 bytes are the timestamp.
    lower   = (ULong64_t)*pevent++;
    upper   = (ULong64_t)*pevent++;
    fScaler->SetInternalTS((ULong64_t)(upper<<16|lower));

    if(fvl>1)
      cout << "UnpackedEvent, DecodeScaler: timestamp " << fScaler->GetInternalTS() << endl;

    //Next 8 bytes are the number of scalers.
    lower   = (ULong64_t)*pevent++;
    upper   = (ULong64_t)*pevent++;
    nScaler = (ULong64_t)(upper<<16|lower);


    if(fvl>1)
      cout << "UnpackedEvent, DecodeScaler: " << nScaler << " scalers found " << endl;
    if(nScaler > NSCALER){
      cout << "UnpackedEvent, DecodeScaler:  error, too many scalers found " << endl;
      cout << "UnpackedEvent, DecodeScaler: " << nScaler << " scalers found " << endl;
      cout << "UnpackedEvent, DecodeScaler: limit (from Scalerdefs.h) "<< NSCALER << endl;
      return 4;
    }

    //Scaler values are now stored every 8 bytes.
    for(UInt_t i=0; i<nScaler; i++){
      lower = (ULong64_t)*pevent++;
      upper = (ULong64_t)*pevent++;
      ULong64_t value = (ULong64_t)(upper<<16|lower);
      if(fvl>1 && i==1)
	cout << i << "\t" << fnentries<< "\t" << (ULong64_t)(upper)<< "\t" << (ULong64_t)(lower)<< "\t" << value << "\t";

      // the scalers are 24 bit, we need to count the overflows and correct for this
      value+=pow(2,24)*foverflows[i];
      // calc the difference between the previous value and the current one
      double diff =  (double)value-(double)fprevscaler[i];
      if(fvl>1 && i==1)
	cout << "this diff " << diff << "\t";

      // check whether we are at a jump
      if(diff < -1e7){
	foverflows[i]++; // one more overflow
	value+=pow(2,24); // fix value
      }

      fScaler->SetValue(i,value);
      fScaler->SetOverflows(i,foverflows[i]);
      fprevscaler[i] = value;
      if(fvl>1 && i==1)
	cout << value << endl;
      if(fvl>1)
	cout <<i<< "UnpackedEvent, DecodeScaler:  val " << fScaler->GetValue(i) << endl;
    }
  } else {
    cout << "Unknown non-event type: " << type << endl;
    return 2;
  }

  return 0;
}

int UnpackedEvent::DecodeMode3(char* cBuf, int len, long long int gts, bool card29){
  if(ffirst_ts<1){
    if(fvl>1)
      cout << "UnpackedEvent: " << "setting first timestamp " << gts << endl;
    ffirst_ts = gts;
    fcurrent_ts = gts;
  }
  fctr++;

  unsigned short *wBuf = (unsigned short*)cBuf;
  HEtoLE(cBuf, len);
  // length now in units of 16bit words
  len/=2;
  if(fvl>0)
    cout << "UnpackedEvent: " << " length " << len <<endl;

  Trace curTrace;
  curTrace.Clear();

  int tracesFound = 0;

  //As long as we still have data in the buffer.
  while(len>0){
    tracesFound++;
    if(fvl>2){
      cout << "UnpackedEvent: " << *wBuf << " " << (hex) << *wBuf << (dec) << endl;
      cout << "UnpackedEvent: " << *(wBuf+1) << " " << (hex) << *(wBuf+1) << (dec) << endl;
    }

    // 1st & 2nd word are 0xaaaa
    if( (*wBuf != 0xaaaa) && (*(wBuf+1) != 0xaaaa) ) {
      cerr << "0xAAAA header missing" << endl;
      return 9;
    } else if(fvl>1){
      cout << "UnpackedEvent: " << "-----------------------------next 0xaaaa aaaa" << endl;
    }

    wBuf+=2;

    int length  = (*wBuf & 0x07ff) * 2 + 2;
    len -= length;

    //Read out from buffer and build current trace.
    curTrace = DecodeTrace(&wBuf,length,gts);

    //Trace is now built from input.
    //Time to add it into an event.


    static bool errorShown = false;
    if(!errorShown && (curTrace.GetLED() < fcurrent_ts)){
      cout << "UnpackedEvent: " << "Found events that are not in order" << endl;
      cout << "UnpackedEvent: " << "Have you run this through GEB_HFC first?" << endl;
      errorShown = true;
    }

    //Enough global time has passed, and so I am now in a new event.
    //Close and write the event, then clear it out.
    if(curTrace.GetLED() - fcurrent_ts > fEventTimeDiff || curTrace.GetLED() < fcurrent_ts){
      if(fvl>2)
	cout << "UnpackedEvent: " << "Closing event due to timestamp in Mode3 Entry." << endl;
      if(fvl>1&&fMode3Event->GetMult()<1)
	cout << "UnpackedEvent: " <<__PRETTY_FUNCTION__ << " entry " << fnentries << " mode3 is empty " << endl;
      if(fvl>1&&fhasdata==false)
	cout << "UnpackedEvent: " <<__PRETTY_FUNCTION__ << " entry " << fnentries << " s800 is empty " << endl;
      fMode3Event->SetCounter(fctr);
      fctr = 0;
      this->CloseEvent();
    }

    //The event is now written to the tree and is cleared if we decided to make a new event.
    //The event, whether newly made or made previously, gets passed a new trace.
    fMode3Event->AddTrace(curTrace);

    if(fvl>1)
      cout << "UnpackedEvent: " << "setting current ts " << curTrace.GetLED() << endl;

    fcurrent_ts = curTrace.GetLED();

    if(fvl>2)
      cout << "UnpackedEvent: " << "remaining " << len << endl;

  }

  if(card29 && tracesFound!=2 && fvl>1){
    cout << "UnpackedEvent: Expected 2 traces in card29 entry, but found " << tracesFound << endl;
  }

  if(fvl>2){
    cout << "UnpackedEvent: " << "Mode3 event found with timestamp " << gts << endl;
    cout << "UnpackedEvent: " << tracesFound << " traces decoded from the mode3 entry" << endl;
  }
  return 0;
}

Trace UnpackedEvent::DecodeTrace(unsigned short** wBuf_p, int length, long long int gts){
  unsigned short* wBuf = *wBuf_p;

  Trace curTrace;
  curTrace.Clear();


  curTrace.SetTS(gts);
  curTrace.SetLength(length-16);
  // this length includes aaaa aaaa until next aaaa aaaa.
  // trace length is this minus 16

  if(fvl>2)
    cout << "UnpackedEvent: " << "length " << curTrace.GetLength() << endl;

  curTrace.SetBoard(*wBuf >> 11);//called GA in mario's doc 13.dez
  wBuf++;
  if(fvl>4){
    cout << "UnpackedEvent: " << "---------------------------------------------------------------------------------"<<endl;
    cout << "UnpackedEvent: " << *wBuf << "\t" << (hex) << *wBuf << (dec) <<"\t";
    bits_uint(*wBuf);
    cout << "UnpackedEvent: " << (*wBuf & 0x000f) << "\t" << (hex) << (*wBuf & 0x000f) << (dec) <<"\t";
    bits_uint(*wBuf & 0x000f);
    cout << "UnpackedEvent: " << ((*wBuf & 0x0030)>>4) << "\t" << (hex) << ((*wBuf & 0x0030)>>4) << (dec) <<"\t";
    bits_uint((*wBuf & 0x0030)>>4);
    cout << "UnpackedEvent: " << ((*wBuf & 0x00c0)>>6) << "\t" << (hex) << ((*wBuf & 0x00c0)>>6) << (dec) <<"\t";
    bits_uint((*wBuf & 0x00c0)>>6);
    cout << "UnpackedEvent: " << ((*wBuf & 0x1f00)>>8) << "\t" << (hex) << ((*wBuf & 0x1f00)>>8) << (dec) <<"\t";
    bits_uint((*wBuf & 0x1f00)>>8);
  }

  curTrace.SetChn(*wBuf & 0x000f);
  curTrace.SetSlot((*wBuf & 0x0030)>>4);
  curTrace.SetCrystal((*wBuf & 0x00c0)>>6);
  curTrace.SetHole((*wBuf & 0x1f00)>>8);
  wBuf++;


  if(fvl>2)
    cout << "UnpackedEvent: " << "board " << curTrace.GetBoard() << " chn " << curTrace.GetChn() << " slot " << curTrace.GetSlot() << " cry " << curTrace.GetCrystal() << " det " << curTrace.GetHole() << endl;

  int id = curTrace.GetHole()*4 + curTrace.GetCrystal();

  // point 5th
  long long int ts;
  ts = (long long int) *(wBuf+1);
  ts += ((long long int) *(wBuf+0)) << 16;
  ts += ((long long int) *(wBuf+3)) << 32;
  curTrace.SetLED(ts);

  if(fvl>2)
    cout << "UnpackedEvent: " << "led " << curTrace.GetLED() << endl;

  wBuf+=2; //point 7th

  int en = (int) *(wBuf+3) & 0x00ff;
  en = en << 16;
  en += (int) *(wBuf);
  int sign = *(wBuf+3) & 0x0100;

  if(sign){
    if(curTrace.GetChn()==9) //core
      en = (int)(en - (int)0x01000000); // (2^24)
    else{
      en = (int)(en - (int)0x01000000); // (2^24)
      en = - en;
    }
  }
  else{
    if(curTrace.GetChn()!=9) // not core
      en = - en;
  }


  curTrace.SetEnergy(en);
  curTrace.SetEnSign(sign);
  curTrace.SetPileUp((*(wBuf+3) & 0x8000) >> 15);
  wBuf+=2; //point 9th

  if(fvl>2)
    cout << "UnpackedEvent: " << "energy " << curTrace.GetEnergy() << " sign " << curTrace.GetEnSign() << " pileup " << curTrace.GetPileUp() << endl;

  ts = (long long int) *(wBuf+0);
  ts += ((long long int) *(wBuf+3)) << 16;
  ts += ((long long int) *(wBuf+2)) << 32;
  curTrace.SetCFD(ts);

  if(fvl>2)
    cout << "UnpackedEvent: " << "cfd " << curTrace.GetCFD() << endl;

  wBuf+=4; //point 13th

  int cfd = (int) *(wBuf+1) << 16;
  cfd += (int) *wBuf;
  curTrace.SetCFD(0,cfd);
  wBuf+=2; //point 15th

  cfd = (int) *(wBuf+1) << 16;
  cfd += (int) *wBuf;
  curTrace.SetCFD(1,cfd);
  wBuf+=2; //point 17th

  if(fvl>2)
    cout << "UnpackedEvent: " << "cfd points " << curTrace.GetCFD(0) << " and " << curTrace.GetCFD(1) << endl;

  //cout << "UnpackedEvent: " << curTrace.GetTrace().size() << " size" << endl;
  int bg =0;
  int sig =0;
  for(int i=0; i<(length-16); i+=2){
    //cout << "UnpackedEvent: " << i <<"\t"<< *(wBuf) << endl;
    curTrace.SetTrace(i ,-(short) *(wBuf+1)+512 );
    curTrace.SetTrace(i+1 ,-(short) *(wBuf)+512 );
    if(id==124){
      curTrace.SetLaBrTrace(i ,-(short) *(wBuf+1)+512 );
      curTrace.SetLaBrTrace(i+1 ,-(short) *(wBuf)+512 );
      if(i>19&&i<71){

	bg+=-(short) *(wBuf)+512;
	bg+=-(short) *(wBuf+1)+512;
	if(fvl>2)
	  cout << "UnpackedEvent: " << "labr bg " << bg << endl;
      }
      if(i>69&&i<121){

	sig+=-(short) *(wBuf)+512;
	sig+=-(short) *(wBuf+1)+512;
	if(fvl>2)
	  cout << "UnpackedEvent: " << "labr sig " << sig << endl;
      }
    }
    wBuf++;
    wBuf++;


    if(fvl>2){
      cout << "UnpackedEvent: " << i <<"\t"<< curTrace.GetTrace()[i] << endl;
      cout << "UnpackedEvent: " << i+1 <<"\t"<< curTrace.GetTrace()[i+1] << endl;
    }
  }
  if(id==124){
    if(fvl>1){
      cout << "UnpackedEvent: " << "end fof trace bg " << bg << " sig " << sig << endl;
      cout << "UnpackedEvent: " << " energy " << sig-bg << endl;
    }
    curTrace.SetLaBr(sig-bg);
  }
  if(fvl>2){
    cout << "UnpackedEvent: " << *wBuf << " " << (hex) << *wBuf << (dec) << endl;
    cout << "UnpackedEvent: " << *(wBuf+1) << " " << (hex) << *(wBuf+1) << (dec) << endl;
  }


  if(curTrace.GetBoard() == 6 && curTrace.GetChn() == 9){ //board 3 = 10 , 6 = 5 MeV range
    if(fvl>1){
      cout << "UnpackedEvent: " << fctr << " core hit id: "<<id<<"\tts: " << curTrace.GetLED() << endl;
    }
    curTrace.SetCore(true);
  }
  else{
    if(fvl>1){
      cout << "UnpackedEvent: " << fctr << " segm hit id: "<<id<<"\tts: " << curTrace.GetLED() << endl;
    }
  }

  *wBuf_p = wBuf;

  return curTrace;
}


void UnpackedEvent::ClearEvent(){

  fS800->Clear();
  fhasdata = false;
  fGretina->Clear();
  fMode3Event->Clear();
  return;
}

void UnpackedEvent::CloseEvent(){
  if(fwtree || fwhist){

    //Perform mode 2 recalibrations, addback.


    if (fwhist){
      frhist->FillHistograms(fGretina,fS800,fMode3Event);
    }
    //Write the raw tree.
    if (fwtree){
      ftr->Fill();
    }
    fnentries++;
  }
  if(fwcaltree||fwcalhist){

    //Build all of the calibrated objects, using the calibration in cal.
    //Use the data from the first three parameters, output into the last three parameters.
    //cout << "calculation BuildAllCalc called" << endl;
    fcal->BuildAllCalc(fS800,fGretina,fMode3Event,
		       fS800Calc,fGretinaCalc,fMode3Calc);

    if(fwcaltree){
      fcaltr->Fill();
      fncalentries++;
    }
    if(fwcalhist){
      fchist->FillHistograms(fGretinaCalc,fS800Calc,fMode3Calc);
    }
  }

  this->ClearEvent();
}

void UnpackedEvent::WriteLastEvent(){
  if(fvl>2)
    cout << "UnpackedEvent: " << "last event " << endl;
  fhits++;
  if(fvl>1&&fMode3Event->GetMult()<1)
    cout << "UnpackedEvent: " <<__PRETTY_FUNCTION__<< " entry " << fnentries << " mode3 is empty " << endl;
  if(fvl>1&&fhasdata==false)
    cout << "UnpackedEvent: " <<__PRETTY_FUNCTION__<< " entry " << fnentries << " s800 is empty " << endl;
  fMode3Event->SetCounter(fctr);
  fctr = 0;
  this->CloseEvent();
  fhits--;

  if(fwhist){
    frhist->Write();
  }
  if(fwcalhist){
    fchist->Write();
  }
}
double UnpackedEvent::GetRunTime(){
  //cout << __PRETTY_FUNCTION__ << endl;
  int n = -1;
  int status = -1;
  long long int first = -1;
  long long int ts0 = -1;
  long long int tg0 = -1;
  long long int tm0 = -1;
  long long int last = -1;
  long long int tsn = -1;
  long long int tgn = -1;
  long long int tmn = -1;
  if (fwtree){
    n = ftr->GetEntries()-1;
    //cout << "number of events ftr->GetEntries()-1 = " << n << endl;
    status = ftr->GetEvent(0);
    //cout << "status = ftr->GetEvent(0) = " << status << endl;
    if(fS800)
      ts0 = fS800->GetTS();
    if(fGretina->GetMult()>0)
      tg0 =  fGretina->GetHit(0)->GetTS();
    if(fMode3Event->GetMult()>0)
      tm0 = fMode3Event->GetHit(0)->GetTS();
    if(ts0>0)
      first = ts0;
    if(tg0>0 && tg0<first)
      first = tg0;
    if(tm0>0 && tm0<first)
      first = tm0;
    //cout << 0 <<"\t"<< status <<"\t"<< ts0<<"\t"<< tg0<<"\t"<< tm0 <<"\t" << first<< endl;
    while(last<0){
      //cout << "getting event n = " << n << endl;
      status = ftr->GetEvent(n);
      if(fS800)
	tsn = fS800->GetTS();
      if(fGretina->GetMult()>0)
	tgn =  fGretina->GetHit(0)->GetTS();
      if(fMode3Event->GetMult()>0)
	tmn = fMode3Event->GetHit(0)->GetTS();
      if(tsn>0)
	last = tsn;
      if(tgn>0 && tgn>last)
	last = tgn;
      if(tmn>0 && tmn>last)
	last = tmn;
      n--;
    }
    //cout << n <<"\t"<< status <<"\t"<< tsn<<"\t"<< tgn<<"\t"<< tmn <<"\t" << last<< endl;
    //if(first>0 && last>0)
    //cout << "diff " << last-first << "\t" << 1e-8*(last-first) << " sec "<<  1e-8*(last-first)/60 << " min "<< endl;
    return 1e-8*(last-first);


  }
  if(fwcaltree){
    n = fcaltr->GetEntries()-1;
    status = fcaltr->GetEvent(0);
    if(fS800Calc)
      ts0 = fS800Calc->GetTimeStamp();
    if(fGretinaCalc->GetMult()>0)
      tg0 =  fGretinaCalc->GetHit(0)->GetTS();
    if(fMode3Calc->GetMult()>0)
      tm0 = fMode3Calc->GetHit(0)->GetTime();
    if(ts0>0)
      first = ts0;
    if(tg0>0 && tg0<first)
      first = tg0;
    if(tm0>0 && tm0<first)
      first = tm0;
    //cout << 0 <<"\t"<< status <<"\t"<< ts0<<"\t"<< tg0<<"\t"<< tm0 <<"\t" << first<< endl;
    while(last<0 && n>0){
      status = fcaltr->GetEvent(n);
      if(fS800Calc)
	tsn = fS800Calc->GetTimeStamp();
      if(fGretinaCalc->GetMult()>0)
	tgn =  fGretinaCalc->GetHit(0)->GetTS();
      if(fMode3Calc->GetMult()>0)
	tmn = fMode3Calc->GetHit(0)->GetTime();
      if(tsn>0)
	last = tsn;
      if(tgn>0 && tgn>last)
	last = tgn;
      if(tmn>0 && tmn>last)
	last = tmn;
      n--;
    }
    //cout << n <<"\t"<< status <<"\t"<< tsn<<"\t"<< tgn<<"\t"<< tmn <<"\t" << last<< endl;
    //if(first>0 && last>0)
    //  cout << "diff " << last-first << "\t" << 1e-8*(last-first) << " sec "<<  1e-8*(last-first)/60 << " min "<< endl;
    return 1e-8*(last-first);




  }
  return 0;
}
