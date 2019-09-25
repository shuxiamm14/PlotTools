#ifndef histSaver_h
#define histSaver_h
#include <iostream>
#include <map>
#include "TH1D.h"
#include "TFile.h"
class histSaver{
public:
  TString inputfilename;
  int nregion;
  float blinding;
  Float_t* fweight;
  Double_t* dweight;
  int weight_type;
  std::vector<TString> overlaysamples;
  TFile* inputfile;
  std::map<TString,TFile*> outputfile;
  TString lumi;
  TString analysis;
  TString workflow;
  bool fromntuple;
  bool doROC;
  TString nominalfilename;
  int histcount;
  TString this_region;
  TString read_path;
  TString createdNP;
  int debug;
  Int_t nbin[50];
  Float_t xlo[50];
  Float_t xhi[50];
  TString titleX[50];
  int rebin[50];
  Float_t* var1[50];
  Double_t* var3[50];
  Int_t* var2[50];
  Bool_t MeVtoGeV[50];
  Int_t nvar;
  TString name[50];
  Double_t xbins[50][101];
  bool ifRebin[50];
  bool dataref;
  TString unit[50];
  TString trexdir;
  TString current_sample;
  std::map<TString, TString> variations; //variations[sample] = variation_name
  std::vector<TString> stackorder;
  TString outputfilename;
  TString sensitivevariable;
  std::map<TString, std::map<TString, std::map<TString, std::vector<TH1D*> > > > plot_lib; //plot_lib[sample][region][variation][var]
  std::vector<TString> regions;
  std::vector<TString> mutedregions;
  static TFile *bufferfile;
  histSaver(TString outputfilename);
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
        printf("unknown var type: %s\n", name[nvar].Data());
      unit[nvar] = unit_;
      MeVtoGeV[nvar] = MeVtoGeV_;
      ifRebin[nvar] = 0;
      nvar++;
    }
  }
  void show();
  bool find_sample(TString sample);
  TH1D* grabbkghist(TString region, int ivar);
  TH1D* grabsighist(TString region, int ivar, TString signal="");
  TH1D* grabdatahist(TString region, int ivar);
  bool add_variation(TString sample,TString variation);
  void printyield(TString region);
  double gethisterror(TH1* hist);
  double templatesample(TString fromregion, TString variation,std::string formula,TString toregion,TString newsamplename,TString newsampletitle,enum EColor color,bool scaletogap, double SF = 1);
  void muteregion(TString region);
  void unmuteregion(TString region);
  void SetLumiAnaWorkflow(TString _lumi, TString _analysis, TString _workflow);
  void write_trexinput(TString NPname = "NOMINAL", TString writeoption = "recreate");
  void overlay(TString _overlaysample);
  TH1D* grabhist(TString sample, TString region, int ivar);
  TH1D* grabhist(TString sample, TString region, TString variation, int ivar);
  TH1D* grabhist(TString sample, TString region, TString varname);
  void merge_regions(TString inputregion1, TString inputregion2, TString outputregion);
  //void add(int nbin_, double xlo_, double xhi_, const char* titleX_, const char* name_, float* var_, bool MeVtoGeV_, char* unit_ = "");
  Float_t getVal(Int_t i);
  void add(Int_t nbin_, const Double_t* xbins_, const char* titleX_, const char* name_, Int_t* var_, const char* unit_ = "");
  void add(Int_t nbin_, const Double_t* xbins_, const char* titleX_, const char* name_, Float_t* var_, Bool_t MeVtoGeV_, const char* unit_ = "");
  void add(const char* titleX_, const char* name_, Float_t* var_, Bool_t MeVtoGeV_, const char* unit_ = "");
  float binwidth(int i);
  void add(const char* titleX_, const char* name_, const char* unit_ = "", int _rebin = 1);
  void read_sample(TString samplename, TString savehistname, TString NPname, TString sampleTitle, enum EColor color, double norm, TFile *_inputfile=0);
  void plot_stack(TString NPname,TString outputdir);
  void fill_hist(TString sample, TString region, TString variation);
  void fill_hist(TString sample, TString region);
  void add_region(TString region);
  void init_sample(TString samplename, TString histname, TString sampleTitle, enum EColor color);
  void set_weight(Float_t* _weight){ fweight = _weight; weight_type = 1;}
  void set_weight(Double_t* _weight){ dweight = _weight; weight_type = 2;}
  void write();
};


#endif
