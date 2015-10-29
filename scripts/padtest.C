char *mainch = "crdcpad0";
char *ddchar = "crdcmaxpad0";
Double_t fonegaus(Double_t *x, Double_t *par);
void PressEnterToContinue(){
  int c;
  printf( "Press ENTER to continue... " );
  fflush( stdout );
  do c = getchar(); 
  while ((c != '\n') && (c != EOF));
}
void doit(int cr){
  TFile *f = new TFile("hist/padcal0.root");
  TGraph *g[224];
  TH1F *h[4][224];
  TCanvas *c = new TCanvas("c","c",0,0,800,800);
  c->Divide(3,2);
  for(int p=70;p<224;p++){
    if(p==100)
      continue;
    h[0][p] = (TH1F*)f->Get(Form("%s_%d_%d_in68Niout64Co",mainch,cr,p));
    h[1][p] = (TH1F*)f->Get(Form("%s_%d_%d_in68Niout67Ni",mainch,cr,p));
    h[2][p] = (TH1F*)f->Get(Form("%s_%d_%d_in69Cuout65Ni",mainch,cr,p));
    h[3][p] = (TH1F*)f->Get(Form("%s_%d_%d_in69Cuout67Cu",mainch,cr,p));
    g[p] = (TGraph*)f->Get(Form("line_%d_%d",cr,p));
    for(int i=0;i<4;i++){
      c->cd(1+i);
      h[i][p]->Draw();
      h[i][p]->GetXaxis()->SetRangeUser(0,300);
    }
    c->cd(5);
    g[p]->Draw("AP*");
    c->Modified();
    c->Update();
    PressEnterToContinue();

  }


}
void check(int cr, int p){
  TFile *f = new TFile("hist/padcal0.root");
  TH1F *h[4];
  TGraph *g;
  TCanvas *c = new TCanvas("c","c",0,0,800,800);
  c->Divide(3,2);
  h[0] = (TH1F*)f->Get(Form("%s_%d_%d_in68Niout64Co",mainch,cr,p));
  h[1] = (TH1F*)f->Get(Form("%s_%d_%d_in68Niout67Ni",mainch,cr,p));
  h[2] = (TH1F*)f->Get(Form("%s_%d_%d_in69Cuout65Ni",mainch,cr,p));
  h[3] = (TH1F*)f->Get(Form("%s_%d_%d_in69Cuout67Cu",mainch,cr,p));
  g = (TGraph*)f->Get(Form("line_%d_%d",cr,p));
  for(int i=0;i<4;i++){
    c->cd(1+i);
    h[i]->Draw();
    h[i]->GetXaxis()->SetRangeUser(0,300);
  }
  c->cd(5);
  g->Draw("AP*");
}
void checkmatch(){
  TFile *f = new TFile("hist/padcal0.root");
  TH1F *h[2][4];
  TCanvas *c = new TCanvas("c","c",0,0,800,800);
  c->Divide(4,2);
  for(int cr=0;cr<2;cr++){
    int p=100;
    h[cr][0] = (TH1F*)f->Get(Form("%s_%d_%d_in68Niout64Co",mainch,cr,p));
    h[cr][1] = (TH1F*)f->Get(Form("%s_%d_%d_in68Niout67Ni",mainch,cr,p));
    h[cr][2] = (TH1F*)f->Get(Form("%s_%d_%d_in69Cuout65Ni",mainch,cr,p));
    h[cr][3] = (TH1F*)f->Get(Form("%s_%d_%d_in69Cuout67Cu",mainch,cr,p));
    for(int i=0;i<4;i++){
      c->cd(1+i+cr*4);
      h[cr][i]->Draw();
      h[cr][i]->GetXaxis()->SetRangeUser(0,300);
    }
  }
}
void checkall(){
  TFile *f = new TFile("hist/hist_padcal.root");
  TH2F *h[2][4];
  TCanvas *c = new TCanvas("c","c",0,0,800,800);
  c->Divide(4,2);
  for(int cr=0;cr<2;cr++){
    h[cr][0] = (TH2F*)f->Get(Form("%s_%d_in68Niout64Co",ddchar,cr));
    h[cr][1] = (TH2F*)f->Get(Form("%s_%d_in68Niout67Ni",ddchar,cr));
    h[cr][2] = (TH2F*)f->Get(Form("%s_%d_in69Cuout65Ni",ddchar,cr));
    h[cr][3] = (TH2F*)f->Get(Form("%s_%d_in69Cuout67Cu",ddchar,cr));
    for(int i=0;i<4;i++){
      c->cd(1+i+cr*4);
      h[cr][i]->Draw("colz");
      h[cr][i]->GetYaxis()->SetRangeUser(0,300);
    }
  }
}
void refit(int cr,int p){
  TFile *f = new TFile("hist/padcal0.root");
  TF1* fit[5];
  TH1F *h[4];
  TGraph *g;
  TCanvas *c = new TCanvas("c","c",0,0,800,800);
  c->Divide(3,2);
  double mean[4];
  h[0] = (TH1F*)f->Get(Form("%s_%d_%d_in68Niout64Co",mainch,cr,p));
  h[1] = (TH1F*)f->Get(Form("%s_%d_%d_in68Niout67Ni",mainch,cr,p));
  h[2] = (TH1F*)f->Get(Form("%s_%d_%d_in69Cuout65Ni",mainch,cr,p));
  h[3] = (TH1F*)f->Get(Form("%s_%d_%d_in69Cuout67Cu",mainch,cr,p));
  for(int i=0;i<4;i++){
    c->cd(1+i);
    h[i]->Rebin(2);
    h[i]->Draw();
    h[i]->GetXaxis()->SetRangeUser(0,300);


    fit[i] = new TF1(Form("fit_%d",i),fonegaus,0,1500,3);
    fit[i]->SetLineColor(3);
    fit[i]->SetLineWidth(1);
    
    fit[i]->SetParameter(0,h[i]->Integral());

    //good for crdc1 76 - 152 failures
    //fit[i]->SetParameter(1,130);
    //fit[i]->SetRange(100,200);

    fit[i]->SetParameter(1,100);
    fit[i]->SetRange(80,150);
 

    // if(i==1)
    //   fit[i]->SetRange(70,140);
    //fit[i]->SetParLimits(1,center-25,center+25);
    
    fit[i]->SetParameter(2,30);
    h[i]->Fit(fit[i],"Rq");
    mean[i] = fit[i]->GetParameter(1);
    
 }
  //Matching channel: 140.083 152.072 143.63 152.636
  //Matching channel: 130.752 142.59 136.227 144.149
  double match[4];
  if(cr==0){
    match[0] = 140.083;
    match[1] = 152.072;
    match[2] = 143.63;
    match[3] = 152.636;
  }
  if(cr==1){
    match[0] = 130.752;
    match[1] = 142.59;
    match[2] = 136.227;
    match[3] = 144.149;
  }
  g = new TGraph(4,mean,match);
  c->cd(5);
  g->Draw("AP*");
  TF1* line = new TF1("line","[0]+[1]*x");
  line->SetParameter(0,0.0);
  line->SetParameter(1,1.0);
  line->SetParLimits(1,0.1,2.0);
  g->Fit(line,"q");
  cout << line->GetParameter(0) << "\t" << line->GetParameter(1) << endl;
  TEnv* output = new TEnv("hist/padcal0.fixed.dat");
  output->SetValue(Form("Crdc.%d.Slope.%03d",cr,p),line->GetParameter(1));
  output->SetValue(Form("Crdc.%d.Offset.%03d",cr,p),line->GetParameter(0));
  output->SaveLevel(kEnvLocal);
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
