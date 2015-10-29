#include "Peaks.hh"
Int_t npeaks = 2;

void Peaks(TFile* file, char* name, char* OutputFile, int min, int max, int np, double res, double Threshold, bool rejectbg){

}

void Peaks(TFile* file, char* name, char* OutputFile, int min, int max, int np, double res, double Threshold, bool rejectbg, char* hfile){
  TH1F *h;
  double thresh=Threshold;
  int binmin= min;
  int binmax= max;
  int npeaks =np;
  int p;
  Double_t par[30] = {0};
  ofstream output(OutputFile, std::ios::out | std::ios::app);
  h=(TH1F*) file->Get(name);
  h->GetXaxis()->SetRange(binmin,binmax); 
  h->Draw();
  TSpectrum *sp = new TSpectrum(np,res);
  sp->SetResolution(res);
  Int_t nfound = sp->Search(h,res,"nobackground",thresh);
  cout << "Found " << nfound << " peaks in spectrum" << endl;
  Float_t *xpeaks = sp->GetPositionX();
  npeaks = 0;
  
  TF1 *bgline = new TF1("bgline","pol1",binmin,binmax);
  h->Fit("bgline","qn");
  //Loop on all found peaks. Eliminate peaks at the background level
  par[0] = bgline->GetParameter(0);
  par[1] = bgline->GetParameter(1);
  if(rejectbg){
    cout << "warning: rejecting bg!" << endl;
    par[0] = 0.;
    par[1] = 0.;
  }
  for(p=0;p<nfound;p++){
    Float_t xp = xpeaks[p];
    Int_t bin = h->GetXaxis()->FindBin(xp);
    Float_t yp = h->GetBinContent(bin);
    cout << " xpos " << xp << endl;
    if(yp-TMath::Sqrt(yp) < par[0]+xp*par[1]){
      cout << "Peak at " << xp << " is too low " << yp << endl;
      continue;
    }
    par[3*npeaks+0+2] = yp;
    par[3*npeaks+1+2] = xp;
    par[3*npeaks+2+2] = res;
    npeaks++;
  }
  TF1 *fit = new TF1("fit",fpeaksbg,binmin,binmax,2+3*npeaks);
  //TVirtualFitter::Fitter(h,20+3*npeaks);
  fit->SetParameters(par);
  for(int i=0; i<npeaks; i++)
    fit->SetParLimits(3*i+4,0,3*res);
  //for(int i=0; i<npeaks; i++)
    //fit->SetParLimits(3*i+3,par[3*i+3]-100,par[3*i+3]+100);
  //fit->Draw("same");
  if(rejectbg){
    fit->FixParameter(0,0.0);
    fit->FixParameter(1,0.0);
  }
  fit->SetNpx(10000);
  fit->SetLineColor(3);
  fit->SetLineWidth(1);
  h->Fit("fit","R0");
  //int order[npeaks];
  //double mean[npeaks];
  vector<int> order;
  order.resize(npeaks);
  vector<double> mean;
  mean.resize(npeaks);
  int SwapCount = 0;
  for(int i=0;i<npeaks;i++){
    mean[i] = 0;
    order[i] = i;
    mean[i] = fit->GetParameter(3+3*i);
    //cout << " before mean " << mean[i] << "  position " << i << endl;
  }
  for(int i=0;i<npeaks;i++){
    for(int j=0;j<(npeaks-1);j++){
      if( mean[j] > mean[j+1] ){
	double temp = mean[j];
	int itemp = order[j];
	mean[j] = mean[j+1];
	order[j] = order[j+1];
	mean[j+1] = temp;
	order[j+1] = itemp;
	SwapCount++;	
      }      
    }

    if(SwapCount == 0)
      break;
    else
      SwapCount = 0;
  }
  /*
  for(int i=0;i<npeaks;i++){
    cout << " after mean " << mean[i] << "  position " << order[i] << endl;
  }
  */
  for(int i=0;i<npeaks;i++)
    output << fit->GetParameter(2+3*order[i]) << "\t" << fit->GetParameter(3+3*order[i]) << "\t" << fit->GetParameter(4+3*order[i])*2.35 << "\t";
  output << fit->GetChisquare() << "\t";
  output << name << endl;
    
  TFile* HFile = new TFile(hfile,"update");
  h->Write("",TObject::kOverwrite);
  fit->Write(Form("fit%s",name),TObject::kOverwrite);
  bgline->Write("",TObject::kOverwrite);
  HFile->Close();
  
}
void indivPeaks(TFile* file, char* name, char* OutputFile, int min, int max, int np, double res, double Threshold, char* hfile){
  TH1F *h;
  double thresh=Threshold;
  int binmin= min;
  int binmax= max;
  int npeaks =np;
  int p;
  Double_t par[3] = {0,0,0};
  ofstream output(OutputFile, std::ios::out | std::ios::app);
  h=(TH1F*) file->Get(name);
  h->GetXaxis()->SetRange(binmin,binmax); 
  h->Draw();
  TSpectrum *sp = new TSpectrum(npeaks,res);
  sp->SetResolution(res);
  Int_t nfound = sp->Search(h,res,"nobackground",thresh);
  cout << "Found " << nfound << " peaks in spectrum" << endl;
  Float_t *xpeaks = sp->GetPositionX();
  npeaks = 0;
  
  int peaksused=0;
  //TF1 *fit[nfound];
  vector<TF1*> fit;
  fit.resize(nfound);
  for (p=0;p<nfound;p++){
    Float_t xp = xpeaks[p];
    Int_t bin = h->GetXaxis()->FindBin(xp);
    Float_t yp = h->GetBinContent(bin);
    cout << " xpos " << xp << endl;
    fit[peaksused] = new TF1(Form("fit_%d",p),fonegaus,xp-100,xp+100,3);
    par[0] = yp;
    par[1] = xp;
    par[2] = res;
    npeaks =1;
    fit[peaksused]->SetParameters(par);
    fit[peaksused]->SetParLimits(2,0.1,100);
    fit[peaksused]->SetParLimits(1,xp-20,xp+20);
    fit[peaksused]->SetNpx(10000);
    fit[peaksused]->SetLineColor(3);
    fit[peaksused]->SetLineWidth(1);
    h->Fit(fit[peaksused],"R0");
    peaksused++;
  }
  //int order[peaksused];
  //double mean[peaksused];
  vector<int> order;
  order.resize(peaksused);
  vector<double> mean;
  mean.resize(peaksused);
  int SwapCount = 0;
  for(int i=0;i<peaksused;i++){
    mean[i] = 0;
    order[i] = i;
    mean[i] = fit[i]->GetParameter(1);
    //cout << " before mean " << mean[i] << "  position " << i << endl;
  }
  for(int i=0;i<peaksused;i++){
    for(int j=0;j<(peaksused-1);j++){
      if( mean[j] > mean[j+1] ){
	double temp = mean[j];
	int itemp = order[j];
	mean[j] = mean[j+1];
	order[j] = order[j+1];
	mean[j+1] = temp;
	order[j+1] = itemp;
	SwapCount++;	
      }      
    }

    if(SwapCount == 0)
      break;
    else
      SwapCount = 0;
  }
  /*
  for(int i=0;i<npeaks;i++){
    cout << " after mean " << mean[i] << "  position " << order[i] << endl;
  }
  */
  for(int i=0;i<peaksused;i++)
    output << fit[order[i]]->GetParameter(0) << "\t" << fit[order[i]]->GetParameter(1) << "\t" << fit[order[i]]->GetParameter(2)*2.35 << "\t";
  if(peaksused==0)
    output << 0 << "\t" << 0 << "\t" << 0 << "\t"<< 0 << "\t" << 0 << "\t" << 0 << "\t";
  output << name << endl;
    
  TFile* HFile = new TFile(hfile,"update");
  h->Write("",TObject::kOverwrite);
  for(int i=0;i<peaksused;i++)
    fit[i]->Write(Form("fit%s_%d",name,i),TObject::kOverwrite);
  HFile->Close();
  
}
void Find(TFile* file, char* name, char* OutputFile, int min, int max, int np, double res, double Threshold, char* hfile, bool dofit, int block){
  TH1F *h;
  double thresh=Threshold;
  int binmin= min;
  int binmax= max;
  int npeaks =np;
  Double_t par[5] = {0};
  ofstream output(OutputFile, std::ios::out | std::ios::app);

  int cddet =0;
  int i=0;
  char temp[10], *ptr;
  ptr = name;       //Hilfspointer
  while(*ptr < '0' || *ptr > '9') 
    ptr++; // Bis zur ersten Zahl
  while(*ptr >= '0' && *ptr <= '9') //Alle aneinanderhängenden Zahlen
    temp[i++] = *ptr++;  // in temp speichern
  temp[i] = '\0'; //0-Terminierung 
  cddet = atoi(temp); //in int umwandeln

  h=(TH1F*) file->Get(name);
  h->GetXaxis()->SetRange(binmin,binmax); 
  h->Draw();
  TSpectrum *sp = new TSpectrum(np,res);
  sp->SetResolution(res);
  Int_t nfound = sp->Search(h,res,"nobackground",thresh);
  cout << "Found " << nfound << " peaks in spectrum" << endl;
  Float_t *xpeaks = sp->GetPositionX();
  TFile* HFile = new TFile(hfile,"update");

  string ringstrip = NULL;

  if(cddet%2 == 0)
    ringstrip = "ring";
  else
    ringstrip = "strip";
  
  if((cddet<2) && (block==0) )
    cddet=0;
  else if((cddet<2) && (block==1) )
    cddet=1;
  else if(block==0)
    cddet=2;
  else if(block==1)
    cddet=3;
  else{
    cout << "error" << endl;
    cddet =1000;
  }
      
  
  if(dofit){
    //double mean[nfound];
    //double sigma[nfound];
    //int order[nfound];
    vector<int> order;
    order.resize(nfound);
    vector<double> mean;
    mean.resize(nfound);
    vector<double> sigma;
    sigma.resize(nfound);
    
    for (int p=0;p<nfound;p++){
      h->GetXaxis()->SetRange((int)(xpeaks[p]-10),(int)(xpeaks[p]+10)); 
      TF1 *fit = new TF1("fit",fonegausbg,xpeaks[p]-10,xpeaks[p]+10,5);
      par[0] = 0.;
      par[1] = 0.;
      par[2] = h->GetBinContent(h->GetXaxis()->FindBin(xpeaks[p]));
      par[3] = xpeaks[p];
      par[4] = 1;
      fit->SetParameters(par);
      //fit->SetParLimits(0,0,0);	  
      //fit->SetParLimits(1,0,0);	  
      fit->SetLineColor(3);
      fit->SetLineWidth(1);
      h->Fit("fit");
      mean[p] = fit->GetParameter(3);
      sigma[p] = fit->GetParameter(4);
      order[p] = p;
      
      fit->Write("",TObject::kOverwrite);
      
    }
    
    int SwapCount = 0;
    for(int i=0;i<nfound;i++){
      for(int j=0;j<(nfound-1);j++){
	if( mean[j] > mean[j+1] ){
	  double temp = mean[j];
	  int itemp = order[j];
	  mean[j] = mean[j+1];
	  order[j] = order[j+1];
	  mean[j+1] = temp;
	  order[j+1] = itemp;
	  SwapCount++;	
	}      
      }
      
      if(SwapCount == 0)
	break;
      else
	SwapCount = 0;
    }
    float fwhm =0;
    for(int i=0;i<npeaks;i++){
      //CD.Quadrant0.ring.0.pos:
      output << "CD.Quadrant" << cddet << "." << ringstrip << "." <<  i <<  ".pos:\t" <<  mean[i] << endl;
      //output << mean[i]  << "\tCD.Quadrant" << cddet << "." << ringstrip << i << endl;
      fwhm += sigma[i]*2.35;
    }
    fwhm/=16.;
    cout << "CD.Quadrant" << cddet << ".FWHM:\t" << fwhm << endl;
  }  
  h->GetXaxis()->SetRange(binmin,binmax); 
  h->Write("",TObject::kOverwrite);
  
  HFile->Close();
  
}

Double_t fpeaksbg(Double_t *x, Double_t *par){
   Double_t result = par[0] + par[1]*x[0];
   for (Int_t p=0;p<npeaks+10;p++) {
      Double_t norm  = par[3*p+2];
      Double_t mean  = par[3*p+3];
      Double_t sigma = par[3*p+4];
      result += norm*TMath::Gaus(x[0],mean,sigma);
   }
   return result;
}

Double_t fpeaksbgpad(Double_t *x, Double_t *par){
  Double_t result = par[0] + par[1]*x[0];
  for (Int_t p=0;p<npeaks;p++) {
    Double_t norm  = par[3*p+2];
    Double_t mean  = par[3*p+3];
    Double_t sigma = par[3*p+4];
    result += norm*TMath::Gaus(x[0],mean,sigma);
  }
  return result;
}

Double_t fonegausbg(Double_t *x, Double_t *par){
  static Float_t sqrt2pi = TMath::Sqrt(2*TMath::Pi()), sqrt2 = TMath::Sqrt(2.);
  Double_t arg;
  /*
  par[0]   background constant
  par[1]   background slope
  par[2]   gauss content
  par[3]   gauss mean
  par[4]   gauss width
  */
  arg = (x[0]-par[3])/(sqrt2*par[4]);
  Double_t result = par[0] + par[1]*x[0];
  result += 1/(sqrt2pi*par[4]) * par[2] * exp(-arg*arg);
  return result;
}

Double_t fonegaus(Double_t *x, Double_t *par){
  static Float_t sqrt2pi = TMath::Sqrt(2*TMath::Pi()), sqrt2 = TMath::Sqrt(2.);
  Double_t arg;
  /*
  par[0]   gauss content
  par[1]   gauss mean
  par[2]   gauss width
  */
  Double_t result =0;
  arg = (x[0]-par[1])/(sqrt2*par[2]);
  result += 1/(sqrt2pi*par[2]) * par[0] * exp(-arg*arg);
  return result;
}

Double_t fmultgausbg(Double_t *x, Double_t *par){
  static Float_t sqrt2pi = TMath::Sqrt(2*TMath::Pi()), sqrt2 = TMath::Sqrt(2.);
  Double_t arg;
  /*
  par[0]   background constant
  par[1]   background slope
  par[2]   gauss0 content
  par[3]   gauss0 mean
  par[4]   gauss0 width
  par[5]   gauss1 content
  par[6]   gauss1 mean
  par[7]   gauss1 width
  ......
  */
  Double_t result = par[0] + par[1]*x[0];

  for (Int_t p=0;p<npeaks;p++) {
    Double_t norm  = par[3*p+2];
    Double_t mean  = par[3*p+3];
    Double_t sigma = par[3*p+4];
    arg = (x[0]-mean)/(sqrt2*sigma);
    result += 1/(sqrt2pi*sigma) * norm * exp(-arg*arg);
  }

  return result;
}
Double_t fmultgaus(Double_t *x, Double_t *par){
  static Float_t sqrt2pi = TMath::Sqrt(2*TMath::Pi()), sqrt2 = TMath::Sqrt(2.);
  Double_t arg;
  /*
  par[0]   gauss0 content
  par[1]   gauss0 mean
  par[2]   gauss0 width
  par[3]   gauss1 content
  par[4]   gauss1 mean
  par[5]   gauss1 width
  ......
  */
  Double_t result = 0;

  for (Int_t p=0;p<npeaks;p++) {
    Double_t norm  = par[3*p+0];
    Double_t mean  = par[3*p+1];
    Double_t sigma = par[3*p+2];
    arg = (x[0]-mean)/(sqrt2*sigma);
    result += 1/(sqrt2pi*sigma) * norm * exp(-arg*arg);
  }

  return result;
}
Double_t fgammagausbg(Double_t *x, Double_t *par){
  static Float_t sqrt2pi = TMath::Sqrt(2*TMath::Pi()), sqrt2 = TMath::Sqrt(2.);
  Double_t arg;
  Double_t result = par[0] + par[1]*x[0];
  
  

  Double_t norm  = par[2];
  Double_t mean  = par[3];
  Double_t sigma = par[4];

  Double_t step = par[5];

  arg = (x[0]-mean)/(sqrt2*sigma);
  result += 1/(sqrt2pi*sigma) * norm * exp(-arg*arg);

  result += step/pow(1+exp(sqrt2*arg),2);

  return result;

}
Double_t fgammastep(Double_t *x, Double_t *par){
  static Float_t sqrt2 = TMath::Sqrt(2.);
  Double_t arg;
  Double_t result = par[0] + par[1]*x[0];
  
  //Double_t norm  = par[2];
  Double_t mean  = par[3];
  Double_t sigma = par[4];

  Double_t step = par[5];
  arg = (x[0]-mean)/(sqrt2*sigma);
  result += step/pow(1+exp(sqrt2*arg),2);

  return result;

}
Double_t fgammagaus(Double_t *x, Double_t *par){
  static Float_t sqrt2pi = TMath::Sqrt(2*TMath::Pi()), sqrt2 = TMath::Sqrt(2.);
  Double_t arg;

  Double_t norm  = par[2];
  Double_t mean  = par[3];
  Double_t sigma = par[4];


  arg = (x[0]-mean)/(sqrt2*sigma);
  Double_t result = 1/(sqrt2pi*sigma) * norm * exp(-arg*arg);

  return result;

}
Double_t fdoublegammagausbg(Double_t *x, Double_t *par){
  static Float_t sqrt2pi = TMath::Sqrt(2*TMath::Pi()), sqrt2 = TMath::Sqrt(2.);
  Double_t arg, arg2;
  Double_t result = par[0] + par[1]*x[0];
  
  

  Double_t norm  = par[2];
  Double_t mean  = par[3];
  Double_t sigma = par[4];

  Double_t step = par[5];

  Double_t norm2  = par[6];
  Double_t mean2  = par[7];

  arg = (x[0]-mean)/(sqrt2*sigma);
  result += 1/(sqrt2pi*sigma) * norm * exp(-arg*arg);

  result += step/pow(1+exp(sqrt2*arg),2);

  arg2 = (x[0]-mean2)/(sqrt2*sigma);
  result += 1/(sqrt2pi*sigma) * norm2 * exp(-arg2*arg2);


  return result;

}
