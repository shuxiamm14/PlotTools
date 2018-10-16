#ifndef histSaver_h
#define histSaver_h

#include "fcnc_include.h"

class histSaver{
public:
  Int_t nbin[50];
  Float_t xlo[50];
  Float_t xhi[50];
  std::string titleX[50];
  Float_t* var1[50];
  Float_t* fweight = NULL;
  Double_t* dweight = NULL;
  int weight_type = 0;
  Int_t* var2[50];
  Bool_t MeVtoGeV[50];
  Int_t nvar;
  std::string name[50];
  Double_t xbins[50][101];
  bool ifRebin[50];
  bool dataref;
  TString unit[50];
  TString current_sample;
  map<TString, map<TString, vector<TH1D*>>> plot_lib;
  vector<TString> regions;
  TString this_region = "nominal";
  histSaver();
  virtual ~histSaver();

  template<typename I, typename T,typename D>
  void add(I nbin_, T xlo_, T xhi_, const char* titleX_, const char* name_, D* var_, bool MeVtoGeV_ = false, const char* unit_ = ""){
    if(nvar>=0 && nvar<50) {
      nbin[nvar] = nbin_; 
      xlo[nvar] = xlo_; 
      xhi[nvar] = xhi_;
      titleX[nvar] = titleX_; 
      name[nvar] = name_;
      TString Dname = typeid(*var_).name();
      if (Dname.Contains("nt")) 
        var2[nvar] = (int*)var_;
      else
        var1[nvar] = (float*)var_;
      unit[nvar] = unit_;
      MeVtoGeV[nvar] = MeVtoGeV_;
      ifRebin[nvar] = 0;
      nvar++;
    }
  }
  void show();
  //void add(int nbin_, double xlo_, double xhi_, const char* titleX_, const char* name_, float* var_, bool MeVtoGeV_, char* unit_ = "");
  Float_t getVal(Int_t i);
  void add(Int_t nbin_, const Double_t* xbins_, const char* titleX_, const char* name_, Int_t* var_, const char* unit_ = "");
  void add(Int_t nbin_, const Double_t* xbins_, const char* titleX_, const char* name_, Float_t* var_, Bool_t MeVtoGeV_, const char* unit_ = "");
  float binwidth(int i);
  void plot_stack();
  void fill_hist(TString sample, TString region);
  void fill_hist(TString sample);
  void fill_hist();
  void add_region(TString region);
  void init_sample(TString samplename,TString sampleTitle, enum EColor color);
  void set_weight(Float_t* _weight){ fweight = _weight; weight_type = 1;}
  void set_weight(Double_t* _weight){ dweight = _weight; weight_type = 2;}

  ClassDef(histSaver,1)
};


#endif