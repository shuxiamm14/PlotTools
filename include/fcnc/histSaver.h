#ifndef histSaver_h
#define histSaver_h

#include "fcnc_include.h"
#include "AtlasStyle.h"
#include "AtlasLabels.h"
class histSaver{
public:
  Int_t nbin[50];
  Float_t xlo[50];
  Float_t xhi[50];
  TString titleX[50];
  TString histfilename = "hists";
  int rebin[50];
  int nregion = 0;
  float blinding = 0;
  Float_t* var1[50];
  Double_t* var3[50];
  Float_t* fweight = NULL;
  Double_t* dweight = NULL;
  int weight_type = 0;
  TString overlaysample = "";
  Int_t* var2[50];
  Bool_t MeVtoGeV[50];
  Int_t nvar;
  TString name[50];
  Double_t xbins[50][101];
  bool ifRebin[50];
  bool dataref;
  TString unit[50];
  TString current_sample;
  map<TString, map<TString, vector<TH1D*>>> plot_lib;
  map<TString, map<TString, vector<TH1D*>>>::iterator iter;
  TFile* inputfile = 0;
  vector<TString> regions;
  vector<TString> mutedregions;
  TString lumi = "#it{#sqrt{s}} = 13TeV, 80 fb^{-1}";
  TString analysis = "FCNC tqH H#rightarrow tautau";
  TString workflow = "work in progress";
  bool fromntuple = 1;
  int histcount = 0;
  TString this_region = "nominal";
  TString read_path = "./" ;
  int debug = 1;
  histSaver();
  virtual ~histSaver();
  void clearhist();
  template<typename I, typename T,typename D>
  void add(I nbin_, T xlo_, T xhi_, const char* titleX_, const char* name_, D* var_, bool MeVtoGeV_ = false, const char* unit_ = ""){
    if(nvar>=0 && nvar<50) {
      nbin[nvar] = nbin_; 
      xlo[nvar] = xlo_; 
      xhi[nvar] = xhi_;
      titleX[nvar] = titleX_; 
      name[nvar] = name_;
      TString Dname = typeid(*var_).name();
      if(debug) printf("fill type: %s\n", Dname.Data());
      if (Dname.Contains("i")) 
        var2[nvar] = (int*)var_;
      else if(Dname.Contains("f")) 
        var1[nvar] = (float*)var_;
      else if(Dname.Contains("d"))
        var3[nvar] = (double*)var_;
      else
        printf("unknown var type: %s\n", name);
      unit[nvar] = unit_;
      MeVtoGeV[nvar] = MeVtoGeV_;
      ifRebin[nvar] = 0;
      nvar++;
    }
  }
  void show();
  void muteregion(TString region);
  void unmuteregion(TString region);
  void SetLumiAnaWorkflow(TString _lumi, TString _analysis, TString _workflow);

  void overlay(TString _overlaysample);
  TH1D* grabhist(TString sample, TString region, int ivar);
  void merge_regions(TString inputregion1, TString inputregion2, TString outputregion);
  //void add(int nbin_, double xlo_, double xhi_, const char* titleX_, const char* name_, float* var_, bool MeVtoGeV_, char* unit_ = "");
  Float_t getVal(Int_t i);
  void add(Int_t nbin_, const Double_t* xbins_, const char* titleX_, const char* name_, Int_t* var_, const char* unit_ = "");
  void add(Int_t nbin_, const Double_t* xbins_, const char* titleX_, const char* name_, Float_t* var_, Bool_t MeVtoGeV_, const char* unit_ = "");
  void add(const char* titleX_, const char* name_, Float_t* var_, Bool_t MeVtoGeV_, const char* unit_ = "");
  float binwidth(int i);
  void add(const char* titleX_, const char* name_, const char* unit_ = "", int _rebin = 1);
  void read_sample(TString samplename, TString histname, TString sampleTitle, enum EColor color, double norm);
  void plot_stack(TString outputdir);
  void fill_hist(TString sample, TString region);
  void fill_hist(TString sample);
  void fill_hist();
  void add_region(TString region);
  void init_sample(TString samplename, TString histname, TString sampleTitle, enum EColor color);
  void set_weight(Float_t* _weight){ fweight = _weight; weight_type = 1;}
  void set_weight(Double_t* _weight){ dweight = _weight; weight_type = 2;}
  void write(TFile *outputfile);
};


#endif