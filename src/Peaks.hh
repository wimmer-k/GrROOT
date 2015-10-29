#ifndef PEAK
#define PEAK
#include <string>
#include <string.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>

#include "TCanvas.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TKey.h"
#include "TH1.h"
#include "TF1.h"
#include "TSpectrum.h"
#include "TVirtualFitter.h"

using namespace std;

void Peaks(TFile*, char*, char*, int, int, int, double, double, bool);
void Peaks(TFile*, char*, char*, int, int, int, double, double, bool, char*);
void indivPeaks(TFile*, char*, char*, int, int, int, double, double, char*);
void Find(TFile*, char*, char*, int, int, int, double, double, char*, bool, int);
Double_t fpeaksbg(Double_t *x, Double_t *par);
Double_t fpeaksbgpad(Double_t *x, Double_t *par);
Double_t fonegausbg(Double_t *x, Double_t *par);
Double_t fonegaus(Double_t *x, Double_t *par);
Double_t fmultgausbg(Double_t *x, Double_t *par);
Double_t fmultgaus(Double_t *x, Double_t *par);
Double_t fgammagausbg(Double_t *x, Double_t *par);
Double_t fgammastep(Double_t *x, Double_t *par);
Double_t fgammagaus(Double_t *x, Double_t *par);
Double_t fdoublegammagausbg(Double_t *x, Double_t *par);
#endif
