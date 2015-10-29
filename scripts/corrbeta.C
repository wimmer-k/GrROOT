void plot(char* name){
  //TFile *f = new TFile("hist/hist_all.root");
  TFile *f = new TFile("hist/hist_all_1611.root");
  TH1F *h[2];
  h[0] = (TH1F*)f->Get(Form("egamdc_FW_%s",name));
  h[1] = (TH1F*)f->Get(Form("egamdc_BW_%s",name));
  h[0]->Rebin(4);
  h[1]->Rebin(4);
  h[0]->Draw("");
  h[0]->SetMinimum(0);
  h[1]->Draw("same");
  h[1]->SetLineColor(2);
}

void correct(double e1, double e2, double beta = 0.327){
  double theta1 = 47.0;
  double theta2 = 105.0;
  
  
  theta1 = cos(theta1*3.141592/180.);
  theta2 = cos(theta2*3.141592/180.);

  double gamma = 1./sqrt(1.-beta*beta);
  
  double EL1 = e1/(gamma*(1.-beta*theta1));
  double EL2 = e2/(gamma*(1.-beta*theta2));

  double d1 = beta*pow(gamma,3) * (1.-beta*theta1) - gamma *theta1;
  double d2 = beta*pow(gamma,3) * (1.-beta*theta2) - gamma *theta2;

  cout << "d1 " << d1 << "\t" << "d2 " << d2 << endl;
  cout << "EL1 " << EL1 << "\t" << "EL2 " << EL2 << endl;
  cout << "t1 " << theta1 << "\t" << "t2 " << theta2 << endl;
  d1*=EL1;
  d2*=EL2;

  double diff = e1-e2;
  double dbeta = diff/(d2 - d1);
  cout << "--------------------------------------------------" << endl;
  cout << "dbeta " << dbeta << " new beta " << beta+dbeta << " corrected to " << e1+d1*dbeta << "  " << e2+d2*dbeta << endl;

}
