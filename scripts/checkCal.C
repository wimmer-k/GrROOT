void plot(int d, int c){
  //TFile *f = new TFile("calout.root");
  TFile *f = new TFile("5MeVcal.root");
  TH1F* h  = (TH1F*)f->Get(Form("hcc_d%d_c%d",d,c));
  h->Draw();
  TF1* func[8];
  for(int i=0;i<8;i++){
    cout << "c " << c << " d " << d << " c*d " << d*4+c << " l " << i << endl;
    func[i] = (TF1*)f->Get(Form("f_%d_%d",d*4+c,i));
    func[i]->Draw("ssame");
  }
}
