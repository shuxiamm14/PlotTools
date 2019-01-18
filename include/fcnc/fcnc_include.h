#ifndef fcnc_include
#define fcnc_include

#include "TROOT.h"
#include "TObject.h"
#include "TChain.h"
#include "TFile.h"
#include "TH1F.h"
#include "TMinuit.h"
#include "THStack.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TLorentzVector.h"
#include "TLine.h"
#include "TLegend.h"
#include "TMVA/Tools.h"
#include "TLatex.h"
#include "TH2D.h"
#include <sstream>
#include <vector>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include "TMVA/Tools.h"
#include "TMVA/Factory.h"
#include "TMVA/Reader.h"
#include "TKey.h"
const Double_t GeV=1000;
const Double_t PI=3.1415926536;
using namespace std;

Double_t significance(Double_t b0, Double_t s0, Double_t db=0);

void SetMax(TH1* h1, TH1* h2, Double_t scale);

void SetMax(THStack* h1, TH1* h2, Double_t scale);

Float_t AtoF(const char* str);

void Copy(TH1F* h1, TH1F* h2);
template<typename T,typename D>
TString CharAppend(T* aa, D* bb){
  stringstream ss;
  string str;
  ss<<aa<<bb;
  str = ss.str();
  return str;
}

template<typename T,typename D>
TString CharAppend(T aa, D bb){
  stringstream ss;
  string str;
  ss<<aa<<bb;
  str = ss.str();
  return str;
}


template<typename T,typename D>
TString CharAppend(T* aa, D bb){
  stringstream ss;
  string str;
  ss<<aa<<bb;
  str = ss.str();
  return str;
}


template<typename T,typename D>
TString CharAppend(T aa, D* bb){
  stringstream ss;
  string str;
  ss<<aa<<bb;
  str = ss.str();
  return str;
}

template<typename T>
void deletepointer(T *&a){
  if (a)
  {
    delete a;
    a = NULL;
  }
}

template<typename T, typename D>
int FindBin(T* ary, int nbins, D nb){
  int hibin = nbins;
  int lowbin = 0;
  int newlimit = -1;
  if (nb>ary[nbins-1] || nb<ary[0])
  {
    printf("Error: exceed range\n");
    return -1;
  }
  for(;;){
    newlimit = (hibin + lowbin)/2;
    if (ary[newlimit]>nb) hibin = newlimit; 
    else if (ary[newlimit]<nb) lowbin = newlimit;
    if(newlimit == (hibin + lowbin)/2) return newlimit;
  }
}
template<typename T>
T Max(T a, T b) {
  if(a>b) return a;
  else return b;
}

#endif
