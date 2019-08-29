#include "fcnc_include.h"

#include <chrono>
#include <unistd.h>
#include <cstdlib>
#include <thread>
#include <sys/ioctl.h>

std::vector<TString> readTovecString(TString filename){
  std::vector<TString> v;
  std::ifstream file(filename.Data());
  std::string ss;
  while(!file.eof()){
    getline(file,ss,'\n');
    v.push_back(TString(ss));
  }
  return v;
}

double integral(TH1 *hist, double xlow, double xhigh, double *error){
  return hist->IntegralAndError(hist->FindBin(xlow), hist->FindBin(xhigh)-1, *error);
}

Double_t significance(Double_t b0, Double_t s0, Double_t db) {
  if(db==0) return sqrt(2*(s0+b0)*log(1+s0/b0)-2*s0);
  else {
    Double_t tmp = b0-db*db;
    Double_t b = 0.5*(tmp+sqrt(pow(tmp,2)+4*db*db*(b0+s0)));
    return sqrt(2*(b0+s0)*log((b0+s0)/b)-(b0+s0)+b-(b-b0)*b0/db/db);
  }
}

void SetMax(TH1* h1, TH1* h2, Double_t scale=1.0) {
  h1->SetMaximum(scale*TMath::Max(h1->GetMaximum(),h2->GetMaximum()));
  h2->SetMaximum(scale*TMath::Max(h1->GetMaximum(),h2->GetMaximum()));
}

void SetMax(THStack* h1, TH1* h2, Double_t scale=1.0) {
  h1->SetMaximum(scale*TMath::Max(h1->GetMaximum(),h2->GetMaximum()));
}

void Copy(TH1F* h1, TH1F* h2) {
  if(h1->GetNbinsX()!=h2->GetNbinsX()) {
    printf("Error: TH1Fs do not have same number of bins\n");
    return;
  }
  for(Int_t i=1; i<=h1->GetNbinsX(); i++) {
    h2->SetBinContent(i,h1->GetBinContent(i));
    h2->SetBinError(i,h1->GetBinError(i));
  }
}

Float_t AtoF(const char* str) {
  Float_t num = 1.;
  // split by '*'
  std::string Str = str;
  for(size_t i=0,n; i <= Str.length(); i=n+1) {
    n = Str.find_first_of('*',i);
    if (n == std::string::npos) n = Str.length();
    std::string tmp = Str.substr(i,n-i);
    num *= atof(tmp.c_str());
  }
  return num;
}


void PrintTime(int timeInSec){

  struct winsize w; 
  ioctl(0, TIOCGWINSZ, &w);
  auto start = chrono::steady_clock::now();
  for(;;){
    sleep(timeInSec);
    auto end = chrono::steady_clock::now();
    stringstream ss;
    ss<<"Elapsed time in seconds : "<< chrono::duration_cast<chrono::seconds>(end - start).count()
    << " sec";

    printf("%*s\n" , w.ws_col, ss.str().c_str());
    
  }
}
