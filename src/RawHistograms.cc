#include "RawHistograms.hh"

#include <iostream>
#include <iomanip>
#include <string.h>
#include <sys/time.h>

#include <sstream>


#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH2S.h"
#include "TH1S.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TMath.h"
#include "TString.h"

#include "S800.hh"
#include "Gretina.hh"
#include "Trace.hh"

using namespace TMath;
using namespace std;

void RawHistograms::Write(){
  fhlist->Sort();
  fhlist->Write();

  for(int i=0; i<16; i++){
    cout << Form("counter[%2d]: %10d (%5.1f%) (%s,%s,%s,%s)",
		 i, fcounter[i],100.*fcounter[i]/fentry,
		 (i/1)%2 ? "G2" : "  ",
		 (i/2)%2 ? "S800" : "    ",
		 (i/4)%2 ? "G3" : "  ",
		 (i/8)%2 ? "C29" : "   ") << endl;
  }
}

void RawHistograms::FillHistograms(Gretina* gr, S800* s800, Mode3Event* m3r){
  fentry++;
  //Determine which of the systems are present in the data.
  bool hascard29 = false;
  bool hasmode3 = false;
  bool hasmode2 = gr->GetMult()!=0;
  bool hass800 = s800->GetTS()>0;
  for(int j=0; j<m3r->GetMult();j++){
    if(m3r->GetHit(j)->GetTrace(0)->GetHole()==31){
      hascard29 = true;
    } else {
      hasmode3 = true;
    }
  }

  //Time and event number differences.
  //Should ensure that all events are properly writtne to file.
  if(hass800){
    static bool firstEvent = true;
    static int s800_evtnr = -1;
    static long long int s800_ts = -1;
    static long long int gret_ts = -1;
    static long long int card29_ts = -1;
    static bool prev_gret_present = false;
    static bool prev_card29_present = false;
    if(!firstEvent){
      Fill("s800_event_diff",
	   20,-9.5,10.5,s800->GetEvtNr()-s800_evtnr);
      Fill("s800_ts_diff",
	   10000,-1e5,1e5,s800->GetTS()-s800_ts);
      if(hasmode2 && prev_gret_present){
	Fill("gret_ts_diff",
	     10000,-1e5,1e5,gr->GetHit(0)->GetTS() - gret_ts);
      }
      if(hascard29 && prev_card29_present){
	Fill("card29_ts_diff",
	     10000,-1e5,1e5,m3r->GetHit(0)->GetTrace(0)->GetLED() - card29_ts);
      }
    }
    s800_evtnr = s800->GetEvtNr();
    s800_ts = s800->GetTS();
    if(hasmode2){
      gret_ts = gr->GetHit(0)->GetTS();
    }
    if(hascard29){
      card29_ts = m3r->GetHit(0)->GetTrace(0)->GetLED();
    }
    prev_gret_present = hasmode2;
    prev_card29_present = hascard29;
    firstEvent = false;
  }


  if(hasmode2){
    FillMode2Histograms(gr);
  }
  if(hasmode3){
    FillMode3Histograms(m3r);
  }
  if(hass800){
    FillS800Histograms(s800);
  }


  //---------------------------gretina------------------------------------------

  Fill("hregistr",
       16,-0.5,15.5,s800->GetTrigger()->GetRegistr());
  bool hascoinc = false;
  bool hascoincWide = false;
  double coinc_center = 140;
  double closest_time = sqrt(-1);
  Trace* closest_core = NULL;
  for (int j=0; j<m3r->GetMult(); j++){
    Trace* core = m3r->GetHit(j)->GetCoreTrace();
    if(core==NULL){
      continue;
    }
    double tdiff = s800->GetTS() - core->GetLED();
    if (!(fabs(tdiff-coinc_center) > fabs(closest_time))){ //Double-negative to help with NaN.
      closest_time = tdiff-coinc_center;
      closest_core = core;
    }
    if((tdiff>116) && (tdiff<170)){
      hascoinc = true;
    }
    if((tdiff>93) && (tdiff<186)){
      hascoincWide = true;
    }
  }
  Fill("closest",
       200,-200,200,closest_time);
  if (hascoinc){
    Fill("hregistr_hascoinc",
	 16,-0.5,15.5,s800->GetTrigger()->GetRegistr());
  }
  if(hascoincWide){
    Fill("hregistr_hascoincWide",
	 16,-0.5,15.5,s800->GetTrigger()->GetRegistr());
  }

  if(hass800 && (closest_core!=NULL)){
    bool singlesbit = s800->GetTrigger()->GetRegistr() & (1<<0);
    bool coincbit = s800->GetTrigger()->GetRegistr() & (1<<1);
    Fill("closest_energy",
	 200,-200,200,closest_time,
	 500,0,1500,closest_core->GetEnergy()/600.0);
    if (singlesbit){
      Fill("closest_energy_singles",
	   200,-200,200,closest_time,
	   500,0,1500,closest_core->GetEnergy()/600.0);
    }
    if (coincbit){
      Fill("closest_coinc",
	   200,-200,200,closest_time);
      Fill("closest_energy_coinc",
	   200,-200,200,closest_time,
	   500,0,1500,closest_core->GetEnergy()/600.0);
    }
    if (singlesbit && (!coincbit)){
      Fill("closest_energy_singles_nocoinc",
	   200,-200,200,closest_time,
	   500,0,1500,closest_core->GetEnergy()/600.0);
    }
  }





  if(hascard29){
    Fill("htdiff_card29cfd_card29led",
	 2000,-1000,1000,
	 m3r->GetHit(0)->GetTrace(0)->GetCFD() - m3r->GetHit(0)->GetTrace(0)->GetLED());
    Fill("hhas_card29",
	 100000,0,fnentries,fentry,
	 3,-0.25,1.25,1.0);
  } else{
    Fill("hhas_card29",
	 100000,0,fnentries,fentry,
	 3,-0.25,1.25,0.0);
  }

  if(hass800 && hascard29){
    Fill("htdiff_s800_card29",
	 2000,-1000,1000,
	 s800->GetTS()-m3r->GetHit(0)->GetTrace(0)->GetTS());
    Fill("htdiff_s800_card29cfd",
	 2000,-1000,1000,s800->GetTS()-m3r->GetHit(0)->GetTrace(0)->GetCFD());
    Fill("htdiff_s800_card29led",
	 2000,-1000,1000,s800->GetTS()-m3r->GetHit(0)->GetTrace(0)->GetLED());
  }
  if(hass800 && hasmode2){
    for(int j=0;j<gr->GetMult();j++){
      Fill("htdiff_s800_gret",
	   2000,-1000,1000,s800->GetTS()-gr->GetHits()[j]->GetITS());
      Fill("htdiff_s800_grethead",
	   2000,-1000,1000,s800->GetTS()-gr->GetHits()[j]->GetTS());
      Fill("htdiff_s800_grettrig",
	   2000,-1000,1000,s800->GetTS()-gr->GetHits()[j]->GetTrigTime());
    }
  }
  if(hasmode2){
    for(int j=0;j<gr->GetMult();j++){
      Fill("htdiff_grethead_gret",
	   2000,-1000,1000,gr->GetHits()[j]->GetTS()-gr->GetHits()[j]->GetITS());
      Fill("htdiff_grethead_grettrig",
	   2000,-1000,1000,gr->GetHits()[j]->GetTS()-gr->GetHits()[j]->GetTrigTime());
    }
  }

  if(hass800 && hasmode3){
    for(int j=0; j<m3r->GetMult(); j++){
      Mode3Hit* hit = m3r->GetHit(j);
      Trace* core = hit->GetCoreTrace();
      if (core != NULL){
	bool singlesbit = s800->GetTrigger()->GetRegistr() & (1<<0);
	bool coincbit = s800->GetTrigger()->GetRegistr() & (1<<1);

	Fill("tdiff_energy",
	     1000,-1000,1000,s800->GetTS()-core->GetLED(),
	     500,0,1500,core->GetEnergy()/600.0);
	if (singlesbit){
	  Fill("tdiff_energy_singles",
	       1000,-1000,1000,s800->GetTS()-core->GetLED(),
	       500,0,1500,core->GetEnergy()/600.0);
	}
	if (coincbit){
	  Fill("tdiff_energy_coinc",
	       1000,-1000,1000,s800->GetTS()-core->GetLED(),
	       500,0,1500,core->GetEnergy()/600.0);
	}
	if (singlesbit && (!coincbit)){
	  Fill("tdiff_energy_singles_nocoinc",
	       1000,-1000,1000,s800->GetTS()-core->GetLED(),
	       500,0,1500,core->GetEnergy()/600.0);
	}
	if (!hascoinc && coincbit){
	  Fill("tdiff_energy_coincbit_nocoinc",
	       1000,-1000,1000,s800->GetTS()-core->GetLED(),
	       500,0,1500,core->GetEnergy()/600.0);
	}
	if (!hascoincWide && coincbit){
	  Fill("tdiff_energy_coincbit_nocoincWide",
	       1000,-1000,1000,s800->GetTS()-core->GetLED(),
	       500,0,1500,core->GetEnergy()/600.0);
	}
      }
    }
  }


  if(hasmode3){
    if(hasmode2){
      Fill("hm2mult_m3mult",
	   30,0,30,m3r->GetMult(),
	   30,0,30,gr->GetMult());
      for(int j=0; j<m3r->GetMult(); j++){
	Fill("hpresent_segmult",
	     41,-0.5,40.5,m3r->GetHit(j)->GetMult());
	Fill("hpresent_ID",
	     150,-0.5,149.5,m3r->GetHit(j)->GetTrace(0)->GetID());
      }
    } else {
      for(int j=0; j<m3r->GetMult(); j++){
	Mode3Hit* hit = m3r->GetHit(j);
	Fill("hmissing_segmult",
	     41,-0.5,40.5,hit->GetMult());
	Fill("hmissing_ID",
	     150,-0.5,149.5,hit->GetTrace(0)->GetID());

	for(int k=0;k<m3r->GetMult();k++){
	  Mode3Hit* crys = m3r->GetHit(k);
	  for(UShort_t j=0;j<crys->GetTrace()->size();j++){
	    Trace* trace = crys->GetTrace(j);
	    int b = trace->GetBoard();
	    int ch = trace->GetChn();
	    if(ch==9&&b==6){
	      Fill("hmissing_gamma_scaled",
		   10000,0,10000,trace->GetEnergy()/600.0);
	    }
	  }
	}

      }
    }
  }

  int lastevent = 0;
  if(hasmode2){
    lastevent += 1;
  }
  if(hass800){
    lastevent += 2;
  }
  if(hasmode3){
    lastevent += 4;
  }
  if(hascard29){
    lastevent += 8;
  }

  Fill("hdetectors_hit",
       16,-0.5,15.5,lastevent);
  Fill("hdetectors_hit_event",
       1e5,0,fnentries,fentry,
       16,-0.5,15.5,lastevent);

  fcounter[lastevent]++;
}

void RawHistograms::FillMode2Histograms(Gretina* gr){
  //Mode 2 histograms

  Fill("hmult",30,0,30,gr->GetMult());

  for(int i=0; i<gr->GetMult(); i++){
    Crystal* cr = gr->GetHit(i);
    Fill("hgamma",
	 10000,0,10000,cr->GetEnergy());
    Fill("hsegsum",
	 10000,0,10000,cr->GetSegmentSum());
    Fill("hgamma_segsum",
	 10000,0,10000,cr->GetEnergy(),
	 10000,0,10000,cr->GetSegmentSum());
    Fill("hgamma_segsum_diff",
	 20000,-10000,10000,cr->GetEnergy()-cr->GetSegmentSum());
    Fill("hgamma_IPsum",
	 10000,0,10000,cr->GetEnergy(),
	 10000,0,10000,cr->GetIPSum());
    Fill("hgamma_IPsum_diff",
	 20000,-10000,10000,cr->GetEnergy()-cr->GetIPSum());
    if (cr->GetEnergy()>300 && cr->GetEnergy()<10000){
      Fill("hgamma_segsum_diff_gated",
	   20000,-10000,10000,cr->GetEnergy()-cr->GetSegmentSum());
    }
    Fill("hm2_segmult",
	 41,-0.5,40.5,cr->GetMult());
    Fill("herror",
	 10,-0.5,9.5,cr->GetError());
    Fill("herror_ID",
	 150,-0.5,149.5,cr->GetID(),
	 10,-0.5,9.5,cr->GetError());
    Fill("hgamma_error",
	 10,-0.5,9.5,cr->GetError(),
	 500,0,10000,cr->GetEnergy());
    Fill(Form("hgamma_d%d_c%d",fSett->Hole2Det(cr->GetHole()),cr->GetCrystal()),
	 3500,0,3500,cr->GetEnergy());
    Fill("hhitpattern",
	 150,-0.5,149.5,cr->GetID());
    Fill("hhitpattern_en",
	 150,-0.5,149.5,cr->GetID(),
	 2000,0,2000,cr->GetEnergy());

    for(int i=0; i<4; i++){
      Fill(Form("hcore_e_%d",i),
	   10000,0,3e4,cr->GetCoreE(i));
    }
    Fill("hprestep",
	 1000,0,1000,cr->GetPreStep());
    Fill("hpoststep",
	 1000,0,1000,cr->GetPostStep());

    for(int j=0; j<cr->GetMult(); j++){
      IPoint* IP = cr->GetIPoint(j);
      int segID = cr->GetID()*36 + IP->GetSeg();
      Fill("hseghitpattern",
	   30*4*36,-0.5,30*4*36-0.5,segID);
      Fill("hseghitpattern_en",
	   30*4*36,-0.5,30*4*36-0.5,segID,
	   2000,0,2000,IP->GetSegEnergy());
    }
  }
}

void RawHistograms::FillMode3Histograms(Mode3Event* m3r){
  //Mode 3 histograms
  Fill("hm3ID",
       150,-0.5,149.5,m3r->GetHit(0)->GetTrace(0)->GetID());
  Fill("hm3mult",
       30,-0.5,29.5,m3r->GetMult());
  for(int k=0;k<m3r->GetMult();k++){
    Mode3Hit* crys = m3r->GetHit(k);
    Fill("hm3_segmult",
	 41,-0.5,40.5,crys->GetMult());
    for(UShort_t j=0;j<crys->GetTrace()->size();j++){
      Trace* trace = crys->GetTrace(j);
      int d = trace->GetHole();
      int cr = trace->GetCrystal();
      int b = trace->GetBoard();
      int ch = trace->GetChn();
      Fill("hm3board",
	   10,-0.5,9.5,b);
      Fill("hm3channel",
	   17,-1,16,ch);
      Fill("hm3crystal",
	   5,-1,4,cr);
      Fill("hm3hole",
	   33,-1,32,d);
      Fill("hm3tdiff",
	   1000,0,1000,trace->GetTDiff());
      Fill(Form("hm3energy_b%d_c%d_cr%d_d%d",b-3,ch,cr,d),
	   10000,0,1e7,trace->GetEnergy());
      if(ch==9&&b==6){
	Fill("hm3gamma",
	     2e5,0,2e6,trace->GetEnergy());
	int d2 = fSett->Hole2Det(d);
	if(d2>-1){
	  Fill(Form("hm3cc_d%d_c%d",d,cr),
	       2e5,0,2e6,trace->GetEnergy());
	  Fill(Form("hm3cc_cal_d%d_c%d",d,cr),
	       10000,0,10000,trace->GetEnergy()/625.0);
	}
      }
    }
  }
}

void RawHistograms::FillS800Histograms(S800* s800){
  bool hass800 = true;
  bool hasrf = s800->GetTimeOfFlight()->GetRF() != 0 &&
    s800->GetTimeOfFlight()->GetRF()!=-1;
  bool hasobj = s800->GetTimeOfFlight()->GetOBJ() != 0 &&
    s800->GetTimeOfFlight()->GetOBJ() != -1;
  bool hasxfp = s800->GetTimeOfFlight()->GetXFP() != 0 &&
    s800->GetTimeOfFlight()->GetXFP() != -1;
  bool hastar = s800->GetTimeOfFlight()->GetTAR() != 0 &&
    s800->GetTimeOfFlight()->GetTAR() != -1;
  bool hastacobj = s800->GetTimeOfFlight()->GetTACOBJ() != 0;
  bool hastacxfp = s800->GetTimeOfFlight()->GetTACXFP() != 0;
  bool hasIC = s800->GetIonChamber()->GetChannels()->size() != 0;
  bool hasCrdc0 = s800->GetCrdc(0)->GetData()->size() != 0;
  bool hasCrdc1 = s800->GetCrdc(1)->GetData()->size() != 0;
  bool hasTCrdc0 = s800->GetTCrdc(0)->GetData()->size() != 0;
  bool hasTCrdc1 = s800->GetTCrdc(1)->GetData()->size() != 0;
  bool hasTPpac = s800->GetTPpac()->GetData()->size() != 0;
  if(hass800){Fill("s800Efficiency",  13,-0.5,12.5,0);}
  if(hasrf){Fill("s800Efficiency",    13,-0.5,12.5,1);}
  if(hasobj){Fill("s800Efficiency",   13,-0.5,12.5,2);}
  if(hasxfp){Fill("s800Efficiency",   13,-0.5,12.5,3);}
  if(hastar){Fill("s800Efficiency",   13,-0.5,12.5,4);}
  if(hastacobj){Fill("s800Efficiency",13,-0.5,12.5,5);}
  if(hastacxfp){Fill("s800Efficiency",13,-0.5,12.5,6);}
  if(hasIC){Fill("s800Efficiency",    13,-0.5,12.5,7);}
  if(hasCrdc0){Fill("s800Efficiency", 13,-0.5,12.5,8);}
  if(hasCrdc1){Fill("s800Efficiency", 13,-0.5,12.5,9);}
  if(hasTCrdc0){Fill("s800Efficiency",13,-0.5,12.5,10);}
  if(hasTCrdc1){Fill("s800Efficiency",13,-0.5,12.5,11);}
  if(hasTPpac){Fill("s800Efficiency", 13,-0.5,12.5,12);}

  int hp0 = s800->GetHodoscope()->GetHitPattern(0);
  int hp1 = s800->GetHodoscope()->GetHitPattern(1);
  for(int i=0;i<16;i++){
    if(hp0 & (1<<i)){
       Fill("hodo_hitpat",32,0,32,i);
    }
    if(hp1 & (1<<i)){
       Fill("hodo_hitpat",32,0,32,i+16);
    }
 
  }
  Fill("hodo_mult",32,0,32,s800->GetHodoscope()->GetData()->size());
  for(UShort_t i=0;i<s800->GetHodoscope()->GetData()->size();i++){
    int ch = s800->GetHodoscope()->GetChannels()->at(i);
    int en = s800->GetHodoscope()->GetData()->at(i);
    Fill(Form("hodo_en_%d",ch),4000,0,4000,en);
    Fill("hodo_vs_ch",32,0,32,ch,2000,0,4000,en);
    Fill("hodo_hitchannels",32,0,32,ch);
    int ix = ch/4; //integer division
    int iy = ch%4;
    Fill("hodo_hit2d",8,0,8,ix,4,0,4,iy);
  }
  

}
