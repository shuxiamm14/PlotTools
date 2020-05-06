#include "histSaver.h"
#include "fcnc_include.h"
#include "TGaxis.h"
#include "TGraph.h"
#include "AtlasStyle.h"
#include "AtlasLabels.h"
#include "HISTFITTER.h"
#include "makechart.h"
using namespace std;
histSaver::histSaver(TString _outputfilename) {
  outputfilename = _outputfilename;
  trexdir = "trexinputs";
  nvar = 0;
  inputfilename = "hists";
  nominalfilename = "";
  nregion = 0;
  blinding = 0;
  fweight = NULL;
  dweight = NULL;
  weight_type = 0;
  doROC = 0;
  inputfile = 0;
  lumi = "#it{#sqrt{s}} = 13TeV, 80 fb^{-1}";
  analysis = "FCNC tqH H#rightarrow tautau";
  workflow = "work in progress";
  fromntuple = 1;
  histcount = 0;
  this_region = "nominal";
  read_path = "./" ;
  debug = 1;
  sensitivevariable = "";
  for(Int_t i=0; i<50; i++) {
    nbin[i] = 1; xlo[i] = 0; xhi[i] = 1; var1[i] = 0; var2[i] = 0; MeVtoGeV[i] = 0; var3[i] = 0;
  }
}

histSaver::~histSaver() {
  if(debug) printf("histSaver::~histSaver()\n");
  for(auto& samp : plot_lib){
    for(auto &reg: samp.second) {
      for(auto &variation: reg.second) {
        for (int i = 0; i < nvar; ++i){
          TH1D *target = variation.second[i];
          if(debug) cout<<"\rdeleting histogram:"<<target->GetName()<<std::endl<<std::flush;
            deletepointer(target);
          if(debug) cout<<"\rdone deleting histogram"<<std::endl<<std::flush;
        }
      }
    }
  }
  if(debug) std::cout<<"plot_lib destructed"<<std::endl;
  deletepointer(inputfile);
  if(debug) std::cout<<"inputfile destructed"<<std::endl;
  for(auto &file : outputfile)
    deletepointer(file.second);
  outputfile.clear();
  printf("histSaver::~histSaver() destructed\n");
}

void histSaver::printyield(TString region){
  printf("Print Yeild: %s\n", region.Data());
  double er;
  for(auto iter: plot_lib){
    TH1D* target = grabhist_int(iter.first,region,0,0);
    if(target){
      printf("%s: %4.3f \\pm %4.3f\n", iter.first.Data(), target->IntegralAndError(1,target->GetNbinsX(), er), er);
    }else{
      printf("Warning: histogram not found: %s, %s, %s\n", iter.first.Data(), region.Data(), name[0].Data());
    }
  }
}

int histSaver::findvar(TString varname){
  for (int i = 0; i < nvar; ++i)
  {
    if(name[i] == varname) return i;
  }
  printf("varname not found: %s\n", varname.Data());
  exit(0);
}

TH1D* histSaver::grabhist_int(TString sample, TString region, int ivar, bool vital){
  return grabhist(sample, region, "NOMINAL", ivar, vital);
}

TH1D* histSaver::grabhist(TString sample, TString region, TString variation, TString varname, bool vital){
  int ivar = -1;
  for (int i = 0; i < nvar; ++i)
  {
    if(varname == name[i]){
      ivar = i;
      break;
    }
  }
  return grabhist(sample, region, variation, ivar, vital);
}

TH1D* histSaver::grabhist(TString sample, TString region, TString variation, int ivar, bool vital){
  if(sample == "data") variation = "NOMINAL";
  if(!find_sample(sample)){
    if(debug) {
      //show();
      printf("histSaver:grabhist  Warning: sample %s not found\n", sample.Data());
    }
    return 0;
  }

  if(plot_lib[sample].find(region) == plot_lib[sample].end()){
    if(debug) printf("histSaver:grabhist  Warning: region %s for sample %s not found\n", region.Data(), sample.Data());
    if(debug) show();
    return 0;
  }
  if(plot_lib[sample][region].find(variation) == plot_lib[sample][region].end()){
    if(debug) printf("histSaver:grabhist  Warning: region %s for sample %s with variation %s not found\n", region.Data(), sample.Data(), variation.Data());
    if(debug) show();
    return 0;
  }    
  if(!plot_lib[sample][region][variation][ivar]) {
    if(debug) printf("histSaver:grabhist  WARNING: empty histogram%s\n", name[ivar].Data());
    if(vital) {
      printf("histSaver:grabhist  ERROR: empty vital histogram%s\n", name[ivar].Data());
      exit(0);
    }
  }
  return plot_lib[sample][region][variation][ivar];
}

TH1D* histSaver::grabhist(TString sample, TString region, TString varname, bool vital){
  int ivar = -1;
  for (int i = 0; i < nvar; ++i)
  {
    if(varname == name[i]){
      ivar = i;
      break;
    }
  }
  return grabhist(sample, region, ivar, vital);
}

TH1D* histSaver::grabbkghist(TString region, int ivar){
  TH1D *hist = 0;
  for(auto iter: stackorder){
    if(find(overlaysamples.begin(),overlaysamples.end(),iter)== overlaysamples.end() && iter != "data"){
      TH1D *target = grabhist(iter,region,ivar);
      if(target){
        if(hist == 0) hist = (TH1D*)target->Clone();
        else hist->Add(target);
      }
    }
  }
  return hist;
}

TH1D* histSaver::grabsighist(TString region, int ivar, TString signal){
  TH1D *hist = 0;
  if(overlaysamples.size()==0) {
    printf("signal not defined, please call overlaysample(\"signal\")\n");
    exit(0);
  }
  return grabhist(signal==""?overlaysamples[0]:signal, region, ivar);
}

TH1D* histSaver::grabdatahist(TString region, int ivar){
  return grabhist("data",region,ivar);
}

void histSaver::add(Int_t nbin_, const Double_t* xbins_, const char* titleX_, const char* name_, Int_t* var_, const char* unit_) {
  if(nvar>=0 && nvar<20) {
    nbin[nvar] = nbin_;
    for(int i=0; i<=nbin_; i++) {
      xbins[nvar][i] = xbins_[i];
    }
    titleX[nvar] = titleX_; var2[nvar] = var_;
    name[nvar] = name_;
    unit[nvar] = unit_;
    ifRebin[nvar] = 1;
    nvar++;
  }
}


void histSaver::add(const char* titleX_, const char* name_, const char* unit_, int _rebin) {
  fromntuple = 0;
  if(nvar>=0 && nvar<50) {
    titleX[nvar] = titleX_;
    name[nvar] = name_ ;
    unit[nvar] = unit_;
    rebin[nvar] = _rebin;
    nvar++;
  }
}

void histSaver::add(Int_t nbin_, const Double_t* xbins_, const char* titleX_, const char* name_, Float_t* var_, Bool_t MeVtoGeV_, const char* unit_) {
  if(nvar>=0 && nvar<20) {
    nbin[nvar] = nbin_;
    for(int i=0; i<=nbin_; i++) {
      xbins[nvar][i] = xbins_[i];
    }
    titleX[nvar] = titleX_; var1[nvar] = var_; MeVtoGeV[nvar] = MeVtoGeV_;
    name[nvar] = name_;
    unit[nvar] = unit_;
    ifRebin[nvar] = 1;
    nvar++;
  }
}

Float_t histSaver::getVal(Int_t i) {
  Float_t tmp = -999999;
  if(i>=0  && i<nvar) {
    if(var1[i])      { if(debug) printf("fill var1\n"); tmp = MeVtoGeV[i] ? *var1[i]/1000 : *var1[i]; }
    else if(var3[i]) { if(debug) printf("fill var3\n"); tmp = MeVtoGeV[i] ? *var3[i]/1000 : *var3[i]; }
    else if(var2[i]) { if(debug) printf("fill var2\n"); tmp = *var2[i]; }
    else printf("error: fill variable failed. no var available\n");
    if(debug == 1) printf("fill value: %4.2f\n", tmp);
  }
  if (!ifRebin[i]){
    if(tmp >= xhi[i]) tmp = xhi[i]*0.999999;
    if(tmp < xlo[i]) tmp = xlo[i];
  }else{
    if(tmp >= xbins[i][nbin[i]]) tmp = xbins[i][nbin[i]]*0.999999;
    if(tmp < xbins[i][nbin[i]]) tmp = xbins[i][nbin[i]];
  }
  return tmp;
}

void histSaver::show(){
  for(auto& iter:plot_lib ){
    printf("histSaver::show()\tsample: %s\n", iter.first.Data());
  }
  for(auto const& region: regions) {
    printf("histSaver::show()\tregion: %s\n", region.Data());
  }
  for (int i = 0; i < nvar; ++i)
  {
    if(var2[i]) printf("histSaver::show()\t%s = %d\n", name[i].Data(), *var2[i]);
    else if(var1[i]) printf("histSaver::show()\t%s = %4.2f\n", name[i].Data(), MeVtoGeV[i] ? *var1[i]/1000 : *var1[i]);
    else if(var3[i]) printf("histSaver::show()\t%s = %4.2f\n", name[i].Data(), MeVtoGeV[i] ? *var3[i]/1000 : *var3[i]);
  }
}

float histSaver::binwidth(int i){
  return (xhi[i]-xlo[i])/nbin[i];
}

void histSaver::merge_regions(TString inputregion1, TString inputregion2, TString outputregion){
  if(debug) printf("histSaver::merge_regions\t %s and %s into %s\n",inputregion1.Data(),inputregion2.Data(),outputregion.Data());
  bool exist = 0;
  bool input1exist = 1;
  bool input2exist = 1;

  for(auto& iter:plot_lib ){
    if(debug) printf("=====================start merging sample %s=====================\n", iter.first.Data());
    input1exist = 1;
    input2exist = 1;
    if(iter.second.find(inputregion1) == iter.second.end()){
      if(debug) printf("histSaver::merge_regions\t inputregion1: %s not found for sample %s\n",inputregion1.Data(), iter.first.Data());
      input1exist = 0;
    }
    if(iter.second.find(inputregion2) == iter.second.end()){
      if(debug) printf("histSaver::merge_regions\t inputregion2: %s not found for sample %s\n",inputregion2.Data(), iter.first.Data());
      input2exist = 0;
    }
    if(input1exist == 0 && input2exist == 0) continue;
    bool outputexist = 0;
    if(iter.second.find(outputregion) != iter.second.end()){
      if(debug) printf("histSaver::merge_regions\t outputregion %s exist, overwrite it\n",outputregion.Data());
      exist = 1;
      for(auto &variation: iter.second[outputregion]){
        for(int i = 0; i < nvar; ++i) deletepointer(variation.second[i]);
        variation.second.clear();
      }
    }
    for(auto &variation : iter.second[inputregion1])
    for (int i = 0; i < nvar; ++i)
    {
      if(input1exist == 1) iter.second[outputregion][variation.first].push_back((TH1D*)iter.second[inputregion1][variation.first][i]->Clone(iter.first + "_" + variation.first+"_"+outputregion+"_"+name[i] + "_buffer"));
      else iter.second[outputregion][variation.first].push_back((TH1D*)iter.second[inputregion2][variation.first][i]->Clone(iter.first + "_" + variation.first+"_"+outputregion+"_"+name[i] + +"_buffer"));
      if(input1exist == 1 && input2exist == 1) {
        iter.second[outputregion][variation.first][i]->Add(iter.second[inputregion2][variation.first][i]);
        if(!iter.second[outputregion][variation.first][i]->Integral()) printf("merged histogram %s is empty\n", iter.second[outputregion][variation.first][i]->GetName());
        if(debug)
          printf("add %s to %s as %s\n", iter.second[inputregion2][variation.first][i]->GetName(),iter.second[inputregion1][variation.first][i]->GetName(),iter.second[outputregion][variation.first][i]->GetName());
      }
    }
  }
  //for(auto& iter:plot_lib ){
  //  if(iter.second.find(outputregion) == iter.second.end()){
  //    printf("histSaver::merge_regions merge regions failed\n");
  //    exit(1);
  //  }
  //}
  if(!exist) regions.push_back(outputregion);
}

void histSaver::init_sample(TString samplename, TString variation, TString sampleTitle, enum EColor color){
  createdNP = variation;

  if(find_sample(samplename)) return;

  if(outputfile.find(variation) == outputfile.end()) outputfile[variation] = new TFile(outputfilename + "_" + variation + ".root","recreate");
  else outputfile[variation]->cd();
  current_sample = samplename;
  
  if(debug) printf("add new sample: %s\n", samplename.Data());
  for(auto const& region: regions) {
    for (int i = 0; i < nvar; ++i){
      TH1D *created = new TH1D(samplename + "_" + variation  + "_" +  region + "_" + name[i] + "_buffer",sampleTitle,nbin[i],xlo[i],xhi[i]);
      created->SetDirectory(0);
      plot_lib[samplename][region][variation].push_back(created);
      if (samplename != "data")
      {
        plot_lib[samplename][region][variation][i]->Sumw2();
        plot_lib[samplename][region][variation][i]->SetFillColor(color);
        plot_lib[samplename][region][variation][i]->SetLineWidth(1);
        plot_lib[samplename][region][variation][i]->SetLineColor(kBlack);
        plot_lib[samplename][region][variation][i]->SetMarkerSize(0);
      }
    }
    if(debug == 1) printf("plot_lib[%s][%s][%s]\n", samplename.Data(), region.Data(), variation.Data());
  }
  if (samplename == "data") dataref = 1;

  if(debug) printf("finished initializing %s\n", samplename.Data() );
}

vector<observable> histSaver::scale_to_data(TString scaleregion, string formula, TString scaleVariable, vector<double> slices, TString variation){
  int nslice = slices.size();
  int ivar = 0;
  for (; ivar < nvar; ++ivar)
  {
    if(name[ivar] == scaleVariable) break;
  }

  if(!nslice) {
    nslice = nbin[ivar]+1;
    slices.push_back(xlo[ivar]);
    for (int i = 0; i < nbin[ivar]; ++i)
    {
      slices.push_back(slices[i]+binwidth(ivar));
    }
  }

  if(outputfile.find(variation) == outputfile.end()) outputfile[variation] = new TFile(outputfilename + "_" + variation + ".root", "recreate");
  else outputfile[variation]->cd();
  vector<TString> tokens = split(formula.c_str()," ");
  if(tokens.size()%2) printf("Error: Wrong formula format: %s\nShould be like: 1 real 1 zll ...", formula.c_str());
  vector<observable> scalefrom;
  vector<observable> scaleto;
  for(int i = 0 ; i < nslice ; i++)
    scalefrom.push_back(observable(0,0));
  scaleto = scalefrom;
  for(auto &sample: plot_lib){
    TH1D *target = grabhist(sample.first,scaleregion,variation,scaleVariable);
    int islice = 0;
    if(target){
      auto iter = find(tokens.begin(),tokens.end(), sample.first.Data());
      if(iter != tokens.end())
      {
        double numb = 0;
        try{
          numb = stof(string((*(iter-1)).Data()));
        } 
        catch(const std::invalid_argument& e){
          printf("Error: Wrong formula format: %s\nShould be like: 1 real 1 zll ...", formula.c_str());
          exit(1);
        }
        
        if(target->GetBinLowEdge(0) > slices[0]) {
          printf("WARNING: slice 1 (%4.2f, %4.2f) is lower than the low edge of the histogram %4.2f, please check variable %s\n", slices[0], slices[1], target->GetBinLowEdge(0), scaleVariable.Data());
        }
        for (int i = 1; i <= nbin[ivar]; ++i)
        {
          if(target->GetBinLowEdge(i) < slices[0]) continue;
          if(islice == nslice-1) break;
          scalefrom[islice] += observable(target->GetBinContent(i),target->GetBinError(i))*numb;
          if(target->GetBinLowEdge(i) >= slices[islice+1]) islice+=1;
        }
      }else{
        for (int i = 1; i <= nbin[ivar]; ++i)
        {
          if(target->GetBinLowEdge(i) < slices[0]) continue;
          if(islice == nslice-1) break;
          if(sample.first == "data") 
            scaleto[islice] += observable(target->GetBinContent(i),target->GetBinError(i));
          else scaleto[islice] -= observable(target->GetBinContent(i),target->GetBinError(i));
          if(target->GetBinLowEdge(i) >= slices[islice+1]) islice+=1;
        }
      }
    }
  }
  vector<observable> scalefactor;
  for (int i = 0; i < nslice; ++i)
  {
    scalefactor.push_back(scaleto[i]/scalefrom[i]);
  }
  printf("region %s, scale variable %s in %d slices:\n", scaleregion.Data(), scaleVariable.Data(), nslice);
  for (int i = 0; i < nslice-1; ++i)
    printf("(%4.2f, %4.2f): %4.2f +/- %4.2f to %4.2f +/- %4.2f, ratio: %4.2f +/- %4.2f\n",slices[i], slices[i+1],scalefrom[i].nominal,scalefrom[i].error,scaleto[i].nominal,scaleto[i].error,scalefactor[i].nominal,scalefactor[i].error);

  return scalefactor;
}

void histSaver::scale_sample(TString scaleregion, string formula, TString scaleVariable, vector<observable> scalefactor, vector<double> slices, TString variation){
  int nslice = slices.size();
  int ivar = 0;
  for (; ivar < nvar; ++ivar)
  {
    if(name[ivar] == scaleVariable) break;
  }
  if(!nslice) {
    nslice = nbin[ivar]+1;
    slices.push_back(xlo[ivar]);
    for (int i = 0; i < nbin[ivar]; ++i)
    {
      slices.push_back(slices[i]+binwidth(ivar));
    }
  }
  vector<TString> tokens = split(formula.c_str()," ");
  if(tokens.size()%2) printf("Error: Wrong formula format: %s\nShould be like: 1 real 1 zll ...", formula.c_str());
  for (int i = 0; i < tokens.size(); ++i){
    if(!(i%2)) continue;
      int islice = 0;
      TH1D *target = grabhist(tokens[i],scaleregion,variation,scaleVariable);
      if(!target) {
        printf("histSaver::scale_sample : WARNING: hist not found grabhist(%s,%s,%s,%s)\n", tokens[i].Data(),scaleregion.Data(),variation.Data(),scaleVariable.Data());
        continue;
      }
      for (int i = 1; i <= nbin[ivar]; ++i)
      {
        if(target->GetBinLowEdge(i) < slices[0]) continue;
        if(islice == nslice-1) break;
        double scaletmp = scalefactor[islice].nominal;
        if(scaletmp<0) scaletmp = 0;
        target->SetBinContent(i,target->GetBinContent(i)*scaletmp);
        if(target->GetBinLowEdge(i) >= slices[islice+1]) islice+=1;
      }
  }
}

map<TString,vector<observable>>* histSaver::fit_scale_factor(vector<TString> *fit_regions, TString *variable, vector<TString> *scalesamples, const vector<double> *slices, TString *variation, vector<TString> *postfit_regions){
  auto *_scalesamples = new map<TString,map<TString,vector<TString>>>();
  auto *_postfit_regions = new map<TString,map<TString,vector<TString>>>();
  for(auto sample: *scalesamples){
    (*_scalesamples)[sample];
    (*_postfit_regions)[sample][sample] = *postfit_regions;
  }
  auto ret = fit_scale_factor(fit_regions, variable, _scalesamples, slices, variation, _postfit_regions);
  delete _scalesamples;
  delete _postfit_regions;
  return ret;
}

map<TString,vector<observable>>* histSaver::fit_scale_factor(vector<TString> *fit_regions, TString *variable, map<TString,map<TString,vector<TString>>> *scalesamples, const vector<double> *slices, TString *_variation, map<TString,map<TString,vector<TString>>> *postfit_regions){
  if(!postfit_regions) postfit_regions = scalesamples;
  auto *scalefactors = new map<TString,vector<observable>>();
  TString variation = _variation? *_variation:"NOMINAL";
  vector<observable> iter;
  int nbins = nbin[findvar(*variable)];
  int ihists = 0;
  vector<int> binslices;
  HISTFITTER* fitter = new HISTFITTER();
  vector<TString> params;
  for(auto sample : *scalesamples) {
    if(sample.second.size()){
      for(auto sfForReg : sample.second){
        fitter->setparam("sf_" + sample.first + "_" + sfForReg.first, 1, 0.1, 0.,2.);
        printf("fit region: ");
        for(auto regions : sfForReg.second) printf(" %s ", regions.Data());
        printf("\n");
        params.push_back("sf_" + sample.first + "_" + sfForReg.first);
      }
    }else{
      fitter->setparam("sf_" + sample.first, 1, 0.1, 0.,2.);
      printf("fit region: All\n");
      params.push_back("sf_" + sample.first);
    }
  }
  for (int i = 0; i < slices->size()-1; ++i)
  {
    auto fitsamples = stackorder;
    fitsamples.push_back("data");
    for(auto sample : fitsamples){
      for(auto reg : *fit_regions){
        TH1D *target = grabhist(sample,reg,variation,*variable);
        if(!target) continue;
        if(ihists == 0) {
          binslices = resolveslices(target,slices);
        }
        bool scale = 0;
        TString SFname = "";
        TString addsample = sample;
        for(auto ssample : *scalesamples) {
          if(ssample.first == sample) {
            if(ssample.second.size()){
              for(auto sfForReg: ssample.second){
                for(auto sfreg: sfForReg.second){
                  if (sfreg == reg)
                  {
                    addsample = sample + "_" + sfForReg.first;
                  }
                }
              }
            }
            SFname = "sf_" + addsample;
          }
        }
        fitter->addfithist(sample,target,binslices[i],binslices[i+1]-1,SFname);
        ihists++;
      }
    }
    Double_t val[10],err[10];
    //fitter->debug();
//============================ do fit here============================
    //fitter->asimovfit(100,nprong[iprong]+"ptbin"+char(ptbin+'0')+".root");
    double chi2 = fitter->fit(val,err,0);
    int ipar = 0;
    for (auto par: params)
    {
      (*scalefactors)[par].push_back(observable(val[ipar],err[ipar]));
      ipar++;
    }
    fitter->clear();
  }
  for(auto samp : *postfit_regions){
    TH1D *target;
    map<TString,vector<TString>> plotregions;
    if(samp.second.size() == 0){
      plotregions["sf_" + samp.first] = *fit_regions;
    }else if(samp.second.size() == 1){
      for(auto sfForReg : samp.second)
        plotregions["sf_" + samp.first] = sfForReg.second;
    }else{
      for(auto sfForReg : samp.second)
        plotregions["sf_" + samp.first + "_" + sfForReg.first] = sfForReg.second;
    }

    for(auto sf : plotregions){
      for(auto reg : sf.second){
        target = grabhist(samp.first,reg,variation,*variable);
        if(!target) continue;
        for (int islice = 0; islice < slices->size()-1; ++islice)
        {
          for (int i = binslices[islice]; i < binslices[islice+1]; ++i)
          {
            target->SetBinContent(i,target->GetBinContent(i) * (*scalefactors)[sf.first][islice].nominal);
            target->SetBinError(i,target->GetBinError(i) * (*scalefactors)[sf.first][islice].nominal);
          }
        }
      }
    }
  }
  printf("fit regions:");
  for(auto reg: *fit_regions){
    printf(" %s ", reg.Data());
  }
  printf("\n");
  printf("fit samples:");
  for(auto param: params){
    printf(" %s ", param.Data());
  }
  printf("\n");
  for (int i = 0; i < slices->size()-1; ++i){
    printf("(%4.2f, %4.2f): ",(*slices)[i], (*slices)[i+1]);
    for(auto par: params){
      printf("%4.2f +/- %4.2f, ", (*scalefactors)[par][i].nominal, (*scalefactors)[par][i].error);
    }
    printf("\n");
  }

  return scalefactors;
}

vector<int> histSaver::resolveslices(TH1D* target, const vector<double> *slices){
  printf("histSaver::resolveslices(): total %d bins, from %f, to %f.\n",target->GetNbinsX(),target->GetBinLowEdge(1), target->GetBinLowEdge(target->GetNbinsX())+target->GetBinWidth(target->GetNbinsX()));
  if(debug){
    printf("slices ( ");
    for (int i = 0; i < slices->size(); ++i)
    {
      printf("%f ", slices->at(i));
    }
    printf(")\n");
  }
  vector<int> ret;
  if(target->GetBinLowEdge(1) > slices->at(0) ){
    printf("WARNING: slice 1 (%4.2f, %4.2f) is lower than the low edge of the histogram %4.2f, please check histogram %s\n", slices->at(0), slices->at(1), target->GetBinLowEdge(0), target->GetName());
  }
  if(target->GetXaxis()->GetXmax() < slices->at(slices->size()-1)) {
          printf("WARNING: last slice (%4.2f, %4.2f) is lower than the low edge of the histogram %4.2f, please check histogram %s\n", slices->at(slices->size()-2), slices->at(slices->size()-1), target->GetXaxis()->GetXmax(), target->GetName());
  }
  int islice = 0;
  for (int i = 1; i <= target->GetNbinsX(); ++i)
  {
    if(target->GetBinLowEdge(i) >= slices->at(islice)) {
      ret.push_back(i);
      islice+=1;
    }
    if(islice == slices->size()) break;
  }
  if(target->GetXaxis()->GetXmax() <= slices->at(slices->size()-1)) ret.push_back(target->GetNbinsX()+1);
  if(debug){
    printf("resolved slices ( ");
    for (int i = 0; i < ret.size(); ++i)
    {
      printf("%d ", ret.at(i));
    }
    printf(")\n");
  }
  return ret;
}

void histSaver::read_sample(TString samplename, TString savehistname, TString variation, TString sampleTitle, enum EColor color, double norm, TFile *_inputfile){

  TFile *readfromfile;

  if(_inputfile) readfromfile = _inputfile;
  else{
    if(!inputfile) inputfile = new TFile(inputfilename + ".root", "read");
    if(!inputfile) inputfile = new TFile(nominalfilename + ".root", "read");
    readfromfile = inputfile;
  }
  TString filename(readfromfile->GetName());
  if (debug == 1) printf("read from file: %s\n", filename.Data());
  if (samplename == "data") dataref = 1;
  TString histnameorig(savehistname + "_");
  for(auto const& region: regions) {
    TString histname;
    if(filename.Contains("NOMINAL") && variation.Contains("Xsec")){
      histname = histnameorig + "NOMINAL_" + region + "_";
    }else{
      histname = histnameorig + variation + "_" + region + "_";
    }
    if (debug == 1)
    {
      printf("read sample %s from %s region\n", samplename.Data(), region.Data());
    }
    if(!(TH1D*)(readfromfile->Get(histname + name[0]))) {
      if(debug) printf("histogram name not found: %s\n", (savehistname + "_" + variation + "_" + region + "_" + name[0]).Data());
      continue;
    }
    if (plot_lib[samplename].find(region) != plot_lib[samplename].end())
    {
      for (int i = 0; i < nvar; ++i)
      {
        if(!(TH1D*)(readfromfile->Get(histname + name[i]))) {
          if(debug) printf("histogram name not found: %s\n", (histname+name[i]).Data());
          printf("plot_lib[%s][%s][%s][%d]\n", samplename.Data(), region.Data(),variation.Data(), i);
          show();
          exit(1);
        }

        TH1D *readhist = (TH1D*)readfromfile->Get(histname + name[i]);
        double tmp = readhist->Integral();
        if(tmp!=tmp){
          printf("Warning: %s->Integral() is nan, skip\n", (histname + name[i]).Data());
          continue;
        }
        if(tmp==0){
          printf("Warning: %s->Integral() is 0, skip\n", (histname + name[i]).Data());
          continue;
        }
        
        plot_lib[samplename][region][variation][i]->Add(readhist,norm);
        if(checkread){
          if(samplename == checkread_sample && region == checkread_region && variation == checkread_variation && i == checkread_variable){
            printf("read histogram %s, + %f\n", (histname + name[i]).Data(), readhist->GetBinContent(checkread_ibin)*norm);
          }
        }
      }
    }else{
      ++histcount;
      for (int i = 0; i < nvar; ++i)
      {
        if(!(TH1D*)(readfromfile->Get(histname + name[i]))) {
          if(debug) printf("histogram name not found: %s\n", (histname + name[i]).Data());
          printf("plot_lib[%s][%s][%s][%d]\n", samplename.Data(), region.Data(),variation.Data(), i);
          show();
          exit(1);
        }
        TH1D *readhist = (TH1D*)readfromfile->Get(histname + name[i]);
        double tmp = readhist->Integral();
        if(tmp!=tmp){
          printf("Warning: %s->Integral() is nan, skip\n", (histname + name[i]).Data());
          continue;
        }
        if(tmp==0){
          printf("Warning: %s->Integral() is 0, skip\n", (histname + name[i]).Data());
          continue;
        }
        if(checkread){
          if(samplename == checkread_sample && region == checkread_region && variation == checkread_variation && i == checkread_variable){
            printf("read histogram %s, + %f\n", (histname + name[i]).Data(), readhist->GetBinContent(checkread_ibin)*norm);
          }
        }
        plot_lib[samplename][region][variation].push_back((TH1D*)(readfromfile->Get(histname + name[i])->Clone()));
        plot_lib[samplename][region][variation][i]->SetName(samplename + "_" + variation + "_" + region + "_" + name[i] + "_buffer");
        plot_lib[samplename][region][variation][i]->Scale(norm);
        plot_lib[samplename][region][variation][i]->SetTitle(sampleTitle);
        plot_lib[samplename][region][variation][i]->SetFillColorAlpha(color,1);
        plot_lib[samplename][region][variation][i]->SetLineWidth(1);
        plot_lib[samplename][region][variation][i]->SetLineColor(kBlack);
        plot_lib[samplename][region][variation][i]->SetMarkerSize(0);
        plot_lib[samplename][region][variation][i]->SetDirectory(0);
        if(histcount == 1){
          nbin[i] = plot_lib[samplename][region][variation][i]->GetNbinsX();
          xlo[i] =  plot_lib[samplename][region][variation][i]->GetXaxis()->GetXmin();
          xhi[i] =  plot_lib[samplename][region][variation][i]->GetXaxis()->GetXmax();
        }
      }
    }
  }
}

void histSaver::add_region(TString region){
  regions.push_back(region);
  nregion += 1;
}

void histSaver::fill_hist(TString sample, TString region, TString variation){
  if (weight_type == 0)
  {
    printf("ERROR: weight not set\n");
  }
  for (int i = 0; i < nvar; ++i){
    double fillval = getVal(i);
    if(fillval!=fillval) {
      printf("Warning: fill val is nan: \n");
      printf("plot_lib[%s][%s][%d]->Fill(%4.2f,%4.2f)\n", sample.Data(), region.Data(), i, fillval, weight_type == 1? *fweight : *dweight);
    }
    if(debug == 1) printf("plot_lib[%s][%s][%s][%d]->Fill(%4.2f,%4.2f)\n", sample.Data(), region.Data(), variation.Data(), i, fillval, weight_type == 1? *fweight : *dweight);
    TH1D *target = grabhist(sample,region,variation,i);
    if(target) target->Fill(fillval,weight_type == 1? *fweight : *dweight);
    else {
      if(!add_variation(sample,variation)) printf("add variation %s failed, sample %s doesnt exist\n", variation.Data(), sample.Data());
      TH1D *target = grabhist(sample,region,variation,i);
      if(target) target->Fill(fillval,weight_type == 1? *fweight : *dweight);
      else printf("add_variation didnt work in fill_hist\n");
    }
  }
}

void histSaver::fill_hist(TString sample, TString region){
  fill_hist(sample, region, "NOMINAL");
}


bool histSaver::find_sample(TString sample){
  if(plot_lib.find(sample) == plot_lib.end()) return 0;
  return 1;
}

bool histSaver::add_variation(TString sample,TString variation){
  if(!find_sample(sample)) return 0;
  if(outputfile.find(variation) == outputfile.end()) outputfile[variation] = new TFile(outputfilename + "_" + variation + ".root", "recreate");
  else outputfile[variation]->cd();
  for (int i = 0; i < nvar; ++i){
    for(auto reg : regions){
      if(plot_lib[sample][reg].begin() == plot_lib[sample][reg].end()) {
        printf("histSaver::add_variation() ERROR: No variation defined yet, cant add new variation\n");
        exit(0);
      }
      TH1D *created = (TH1D*) plot_lib[sample][reg][createdNP].at(i);
      if(!created){
        printf("histSaver::add_variation() ERROR: hist doesn't exist: plot_lib[%s][%s][%s][%d]\n",sample.Data(), reg.Data(), plot_lib[sample][reg].begin()->first.Data(),i);
        exit(0);
      }
      created = (TH1D*) created->Clone(sample + "_" + variation + "_" + reg + "_" + name[i] + "_buffer");
      created->Reset();
      created->SetDirectory(0);
      plot_lib[sample][reg][variation].push_back(created);
    }
  }
  return 1;
}

void histSaver::write(){
  for(auto& iter: outputfile){
    for(auto& sample : plot_lib){
      for(auto& region: sample.second) {
        for(auto& variation : region.second){
          if(variation.first != iter.first) {
            continue;
          }
          double tmp = variation.second[0]->Integral();
          if(tmp == 0) continue;
          if(tmp != tmp) {
            printf("Warning: hist integral is nan, skip writing for %s\n", variation.second[0]->GetName());
            continue;
          }
          outputfile[variation.first]->cd();
          for (int i = 0; i < nvar; ++i){
            //if(grabhist(iter.first,region,i)->Integral() == 0) {
            //  printf("Warning: histogram is empty: %s, %s, %d\n", iter.first.Data(),region.Data(),i);
            //}
            TString writename = variation.second[i]->GetName();
            writename.Remove(writename.Sizeof()-8,7); //remove "_buffer"
            printf("write histogram: %s\n", writename.Data());
            variation.second[i]->Write(writename,TObject::kWriteDelete);
          }
        }
      }
    }
    iter.second->Close();
    printf("histSaver::write() Written to file %s\n", iter.second->GetName());
  }
}

void histSaver::write_trexinput(TString NPname, TString writename, TString writeoption){
  if (writename == "")
  {
    writename = NPname;
  }
  gSystem->mkdir(trexdir);
  for (int i = 0; i < nvar; ++i){
    gSystem->mkdir(trexdir + "/" + name[i]);
    for(auto const& region: regions) {
      bool muted = 0;
      for (auto const& mutedregion: mutedregions)
      {
        if(region.Contains(mutedregion))
          muted = 1;
      }
      if(muted) continue;

      gSystem->mkdir(trexdir + "/" + name[i] + "/" + region);
      for(auto& iter : plot_lib){
        if(NPname != "NOMINAL" && iter.first.Contains("data")) continue;
        TString filename = trexdir + "/" + name[i] + "/" + region + "/" + iter.first + ".root";
        TFile outputfile(filename, writeoption);
        if(debug) printf("Writing to file: %s, histoname: %s\n", filename.Data(), NPname.Data());
        TH1D *target = grabhist(iter.first,region,NPname,i);
        if(target) {
          target->Write(writename,TObject::kWriteDelete);
          if(!target->Integral()) printf("Warinig: plot_lib[%s][%s][%d] is empty\n", iter.first.Data(),region.Data(),i);
        }
        else if(debug) printf("Warning: histogram plot_lib[%s][%s][%d] not found\n", iter.first.Data(),region.Data(),i);
        outputfile.Close();
      }
    }
  }
}
void histSaver::clearhist(){
  if(debug) printf("histSaver::clearhist()\n");
  for(auto& sample : plot_lib){
    for(auto& region: sample.second) {
      for(auto& variation : region.second){
        for (int i = 0; i < nvar; ++i){
          if(variation.second[i]){
            variation.second[i]->Reset();
          }else{
            printf("histSaver::Reset() Error: histogram not found: sample: %s, variable: %s, region: %s\n",sample.first.Data(), name[i].Data(),region.first.Data());
          }
        }
      }
    }
  }
}

void histSaver::overlay(TString _overlaysample){
  overlaysamples.push_back(_overlaysample);
}

double histSaver::templatesample(TString fromregion, TString variation,string formula,TString toregion,TString newsamplename,TString newsampletitle,enum EColor color, bool scaletogap, double SF){

  if(outputfile.find(variation) == outputfile.end()) outputfile[variation] = new TFile(outputfilename + "_" + variation + ".root", "recreate");
  else outputfile[variation]->cd();
  istringstream iss(formula);
  vector<string> tokens{istream_iterator<string>{iss},
    istream_iterator<string>{}};
  if(tokens.size()%2) printf("Error: Wrong formula format: %s\nShould be like: 1 data -1 real -1 zll ...", formula.c_str());
  vector<TH1D*> newvec;
  observable scaleto(0,0);
  for (int ivar = 0; ivar < nvar; ++ivar)
  {
    newvec.push_back((TH1D*)grabhist(tokens[1],fromregion, tokens[1] == "data" ? "NOMINAL" : variation,ivar)->Clone(newsamplename+"_"+toregion+name[ivar]));
    newvec[ivar]->Reset();
    newvec[ivar]->SetNameTitle(newsamplename,newsampletitle);
    newvec[ivar]->SetFillColor(color);
  }
  for (int i = 0; i < tokens.size()/2; ++i)
  {
    int icompon = 2*i;
    float numb = 0;
    try{
      numb = stof(tokens[icompon]);
    } 
    catch(const std::invalid_argument& e){
      printf("Error: Wrong formula format: %s\nShould be like: 1 data -1 real -1 zll ...", formula.c_str());
      exit(1);
    }
    if(grabhist(tokens[icompon+1],toregion, tokens[icompon+1] == "data" ? "NOMINAL" : variation,0)){
      if(scaletogap) {
        double error = 0;
        observable tmp(grabhist(tokens[icompon+1],toregion, tokens[icompon+1] == "data" ? "NOMINAL" : variation,0)->Integral(),gethisterror(grabhist(tokens[icompon+1],toregion, tokens[icompon+1] == "data" ? "NOMINAL" : variation,0)));
        scaleto += tmp*numb;
      }
      for (int ivar = 0; ivar < nvar; ++ivar)
      {
        newvec[ivar]->Add(grabhist(tokens[icompon+1],fromregion, tokens[icompon+1] == "data" ? "NOMINAL" : variation,ivar),numb);
      }
    }
  }
  observable scalefactor;
  if(scaletogap) {
    observable scalefrom(newvec[0]->Integral(),gethisterror(newvec[0]));
    scalefactor = scaleto/scalefrom;
    printf("scale from %s: %4.2f +/- %4.2f\nto %s: %4.2f +/- %4.2f\nratio: %4.2f +/- %4.2f\n\n",
      fromregion.Data(), scalefrom.nominal, scalefrom.error,
      toregion.Data(), scaleto.nominal, scaleto.error,
      scalefactor.nominal, scalefactor.error);
    for(auto & hists : newvec){
      hists->Scale(scalefactor.nominal);
    }
  }else{
    for(auto & hists : newvec){
      hists->Scale(SF);
    }
  }
  for(int ivar = 0; ivar < nvar; ivar++){
    plot_lib[newsamplename][toregion][variation] = newvec;
  }
  return scalefactor.nominal;
}

double histSaver::gethisterror(TH1* hist){
  double error = 0;
  for (int i = 0; i < hist->GetNbinsX(); ++i)
  {
    error+=pow(hist->GetBinError(i+1),2);
  }
  return sqrt(error);
}

void histSaver::muteregion(TString region){
  mutedregions.push_back(region);
}

void histSaver::unmuteregion(TString region){
  std::vector<TString>::iterator it = find(mutedregions.begin(), mutedregions.end(), region);
  if(it != mutedregions.end()) mutedregions.erase(it);
  else printf("histSaver::unmuteregion WARNING: region %s is not in the mute list\n",region.Data());
}

void histSaver::SetLumiAnaWorkflow(TString _lumi, TString _analysis, TString _workflow){
  lumi = _lumi;
  analysis = _analysis;
  workflow = _workflow;
}
void histSaver::plot_stack(TString NPname, TString outdir){
  SetAtlasStyle();
  TGaxis::SetMaxDigits(3);
  LatexChart* yield_chart = new LatexChart("yield");
  LatexChart* sgnf_chart = new LatexChart("significance");
  gSystem->mkdir("plots_" + outdir);
  TCanvas cv("cv","cv",600,600);
  TGraph* ROC;
  TH1D *ROC_sig = 0;
  TH1D *ROC_bkg = 0;
  TH1D* datahistorig = 0;
  vector<TH1D*> buffer;
  if(doROC){
    ROC = new TGraph();
    ROC -> SetName("ROC");
    ROC -> SetTitle("ROC");
  }
  TFile *outputrocfile = new TFile (outputfilename + "_roc.root", "recreate");
  for(auto const& region: regions) {
    bool muted = 0;
    for (auto const& mutedregion: mutedregions)
    {
      if(region.Contains(mutedregion))
        muted = 1;
    }
    if(muted) continue;
    gSystem->mkdir("plots_" + outdir + "/" + region);
    for (int i = 0; i < nvar; ++i){
      cv.SaveAs("plots_" + outdir + "/" + region + "/" + name[i] + ".pdf[");
      TPad *padlow = new TPad("lowpad","lowpad",0,0,1,0.3);
      TPad *padhi  = new TPad("hipad","hipad",0,0.3,1,1);
      TH1D hmc("hmc","hmc",nbin[i]/rebin[i],xlo[i],xhi[i]);
      TH1D hmcR("hmcR","hmcR",nbin[i]/rebin[i],xlo[i],xhi[i]);
      TH1D hdataR("hdataR","hdataR",nbin[i]/rebin[i],xlo[i],xhi[i]);
      cv.cd();
      padhi->Draw();
//===============================upper pad bkg and unblinded data===============================
      padhi->SetBottomMargin(0.017);
      padhi->SetRightMargin(0.08);
      padhi->SetLeftMargin(0.12);
      padhi->cd();
      hmc.Sumw2();
      THStack *hsk = new THStack(name[i].Data(),name[i].Data());
      TLegend* lg1 = 0;
      lg1 = new TLegend(0.45,0.7,0.90,0.9,"");
      lg1->SetNColumns(2);
      TH1D *histoverlay;
      if(debug) printf("set hists\n");
      for(auto& iter:stackorder ){
        if(iter == "data") continue;
        if(debug) {
          printf("plot_lib[%s][%s][%d]\n", iter.Data(), region.Data(), i);
        }
        if(grabhist(iter,region,NPname,i)) buffer.push_back((TH1D*)grabhist(iter,region,NPname,i)->Clone());
        else continue;
        if(doROC && sensitivevariable == name[i])
        {
          if(!ROC_bkg) ROC_bkg = (TH1D*) buffer.back()->Clone();
          else ROC_bkg->Add(buffer.back());
        }
        if(rebin[i] != 1) buffer.back()->Rebin(rebin[i]);
        hsk->Add(buffer.back());
        hmc.Add(buffer.back());
        lg1->AddEntry(buffer.back(),buffer.back()->GetTitle(),"F");
      }
      if(!hsk->GetMaximum()){
        printf("histSaver::plot_stack(): ERROR: stack has no entry for region %s, var %s, continue\n", region.Data(), name[i].Data());
        continue;
      }
      double histmax = hmc.GetMaximum() + hmc.GetBinError(hmc.GetMaximumBin());

      TH1D * datahist = 0;
      if(debug) printf("set data\n");
      if (dataref) {
        datahistorig = grabhist("data",region,"NOMINAL",i);
        if(datahistorig) datahist = (TH1D*)datahistorig->Clone("dataClone");
        if(!datahist) {
          printf("histSaver::plot_stack(): WARNING: clone data histogram failed: region %s, variable %s\n", region.Data(), name[i].Data());
          exit(0);
        } 
        if(rebin[i] != 1)
          datahist->Rebin(rebin[i]);
        if(datahist->Integral() == 0) printf("Warning: data hist is empty\n");
        lg1->AddEntry(datahist,"data","LP");
        datahist->SetMarkerStyle(20);
        datahist->SetMarkerSize(0.4);
        datahist->SetMinimum(0);
        histmax = max(histmax, datahist->GetMaximum() + datahist->GetBinError(datahist->GetMaximumBin()));
      }else{
        hsk->SetMinimum(0);
      }

      if(debug) printf("set overlay\n");
      int ratio = 0;

      if(debug) printf("set hsk\n");
      hsk->SetMaximum(1.35*histmax);

      hsk->Draw("hist");
      hsk->GetXaxis()->SetTitle(unit[i] == "" ? titleX[i].Data() : (titleX[i] + " [" + unit[i] + "]").Data());
      hsk->GetXaxis()->SetLabelColor(kWhite);
      char str[30];
      sprintf(str,"Events / %4.2f %s",binwidth(i)*rebin[i], unit[i].Data());
      hsk->GetYaxis()->SetTitle(str);
      hsk->GetYaxis()->SetTitleOffset(1.6);
      hsk->GetYaxis()->SetLabelSize(hsk->GetYaxis()->GetLabelSize()*0.7);
      hsk->GetXaxis()->SetLabelSize(hsk->GetXaxis()->GetLabelSize()*0.7);
      hsk->GetYaxis()->SetTitleSize(hsk->GetYaxis()->GetTitleSize()*0.7);
      if(debug) printf("set blinding\n");

      for(Int_t j=1; j<nbin[i]+1; j++) {
        hmcR.SetBinContent(j,1);
        hmcR.SetBinError(j,hmc.GetBinContent(j)>0 ? hmc.GetBinError(j)/hmc.GetBinContent(j) : 0);
        if(dataref) hdataR.SetBinContent(j, hmc.GetBinContent(j)>0 ? datahist->GetBinContent(j)/hmc.GetBinContent(j) : 1);
        if(dataref) hdataR.SetBinError(j, ( datahist->GetBinContent(j)>0 && hmc.GetBinContent(j)>0 )? datahist->GetBinError(j)/hmc.GetBinContent(j) : 0);
      }

      if(debug) printf("setting hmcR\n");
      hmc.SetFillColor(1);
      hmc.SetLineColor(0);
      hmc.SetMarkerStyle(1);
      hmc.SetMarkerSize(0);
      hmc.SetMarkerColor(1);
      hmc.SetFillStyle(3004);
      hmc.Draw("E2 same");
      lg1->Draw("same");

      if(debug) printf("atlas label\n");
      ATLASLabel(0.15,0.900,workflow.Data(),kBlack,lumi.Data(), analysis.Data(), region.Data());
//===============================blinded data===============================
      if(blinding && dataref){
        for(auto overlaysample: overlaysamples){
          TH1D* histoverlaytmp = (TH1D*)grabhist(overlaysample,region,NPname,i);
          if(!histoverlaytmp){
            printf("histSaver::plot_stack(): Warning: signal hist %s not found\n", overlaysample.Data());
            continue;
          }
          for(Int_t j=1; j<nbin[i]+1; j++) {
            if(histoverlaytmp->GetBinContent(j)/sqrt(datahist->GetBinContent(j)) > blinding) {
              datahist->SetBinContent(j,0);
              datahist->SetBinError(j,0);
              hdataR.SetBinContent(j,0);
              hdataR.SetBinError(j,0);
            }
          }
        }
        if(sensitivevariable == name[i]){
          for(int j = nbin[i]*3/4/rebin[i] ; j <= nbin[i] ; j++){
            datahist->SetBinContent(j,0);
            datahist->SetBinError(j,0);
            hdataR.SetBinContent(j,0);
            hdataR.SetBinError(j,0);
          }
        }
      }
      if(dataref) {
        datahist->Draw("E same");
      }

//===============================lower pad===============================
      padlow->SetFillStyle(4000);
      padlow->SetGrid(1,1);
      padlow->SetTopMargin(0.03);
      padlow->SetBottomMargin(0.35);
      padlow->SetRightMargin(0.08);
      padlow->SetLeftMargin(0.12);
      padlow->cd();

      if(debug) printf("plot data ratio\n");
      if(dataref) {
        hdataR.SetMarkerStyle(20);
        hdataR.SetMarkerSize(0.8);
        hdataR.SetMaximum(1.5);
        hdataR.SetMinimum(0.5);
        hdataR.GetYaxis()->SetNdivisions(504,false);
        hdataR.GetYaxis()->SetTitle("Data/Bkg");
        hdataR.GetYaxis()->SetTitleOffset(hdataR.GetYaxis()->GetTitleOffset()*1.08);
        hdataR.GetYaxis()->CenterTitle();
        hdataR.GetXaxis()->SetTitle(unit[i] == "" ? titleX[i].Data() : (titleX[i] + " [" + unit[i] + "]").Data());
        hdataR.GetXaxis()->SetTitleSize(hdataR.GetXaxis()->GetTitleSize()*0.7);
        hdataR.GetYaxis()->SetTitleSize(hdataR.GetYaxis()->GetTitleSize()*0.7);
      }
      hmcR.SetFillColor(1);
      hmcR.SetLineColor(0);
      hmcR.SetMarkerStyle(1);
      hmcR.SetMarkerSize(0);
      hmcR.SetMarkerColor(1);
      hmcR.SetFillStyle(3004);
      if(debug) printf("plot data ratio\n");
      hdataR.Draw("E same");
      hdataR.GetXaxis()->SetTitleOffset(3.4);
      hdataR.GetXaxis()->SetLabelSize(hdataR.GetXaxis()->GetLabelSize()*0.7); 
      hdataR.GetYaxis()->SetLabelSize(hdataR.GetYaxis()->GetLabelSize()*0.7); 
      if(debug) printf("plot mc ratio\n");
      hmcR.Draw("E2 same");
      TLine line;
      line.SetLineColor(2);
      line.DrawLine(hdataR.GetBinLowEdge(1), 1., hdataR.GetBinLowEdge(hdataR.GetNbinsX()+1), 1.);
      cv.cd();
      if(debug) printf("draw low pad\n");
      padlow->Draw();
      if(debug) printf("printing\n");

//===============================upper pad signal===============================

      padhi->cd();
      if(!overlaysamples.size()) {
        cv.SaveAs("plots_" + outdir + "/" + region + "/" + name[i] + ".pdf");
      }
      std::string regtitle = region.Data();
      findAndReplaceAll(regtitle,"reg","");
      findAndReplaceAll(regtitle,"1l1tau1b3j_","TTH $\\tlhad$ ");
      findAndReplaceAll(regtitle,"1l1tau1b2j_","STH $\\tlhad$ ");
      findAndReplaceAll(regtitle,"1l1tau2b3j_","TTH $\\tlhad$ 2b ");
      findAndReplaceAll(regtitle,"1l1tau2b2j_","STH $\\tlhad$ 2b ");
      findAndReplaceAll(regtitle,"1l2tau1bnj_","$l\\thadhad$ ");
      findAndReplaceAll(regtitle,"1l2tau2bnj_","$l\\thadhad$ 2b ");
      findAndReplaceAll(regtitle,"2lSS1tau1bnj_","$2lSS\\thad$ ");
      findAndReplaceAll(regtitle,"2lSS1tau2bnj_","$2lSS\\thad$ 2b ");

      if(sensitivevariable == name[i]) {
        if(dataref){
          yield_chart->set("data",regtitle,integral(datahistorig));
        }
        yield_chart->set("background",regtitle,integral(&hmc));
        printf("Region %s, Background yield: %f\n", region.Data(), hmc.Integral());
      }
      for(auto overlaysample: overlaysamples){
        
        TLegend *lgsig = (TLegend*) lg1->Clone();
        if(debug) { printf("overlay: %s\n", overlaysample.Data()); }
        if(grabhist(overlaysample,region,NPname,i)) histoverlay = (TH1D*)grabhist(overlaysample,region,NPname,i)->Clone();
        if(doROC && sensitivevariable == name[i]) ROC_sig = (TH1D*) histoverlay->Clone();
        if(!histoverlay) continue;
        if(rebin[i] != 1) histoverlay->Rebin(rebin[i]);
        histoverlay->SetLineStyle(9);
        histoverlay->SetLineWidth(3);
        histoverlay->SetLineColor(kRed);
        histoverlay->SetFillColor(0);
        histoverlay->SetMinimum(0);
        ratio = histmax/histoverlay->GetMaximum()/3;
        if(ratio>10) ratio -= ratio%10;
        if(ratio>100) ratio -= ratio%100;
        if(ratio>1000) ratio -= ratio%1000;
        lgsig->AddEntry(histoverlay,(histoverlay->GetTitle() + (ratio > 0? "#times" + to_string(ratio) : "")).c_str(),"LP");

        if(sensitivevariable == name[i]){
          double _significance = 0;
          for(Int_t j=1; j<nbin[i]+1; j++) {
            if(histoverlay->GetBinContent(j) && hmc.GetBinContent(j)) {
              if(hmc.GetBinContent(j) > 0 && histoverlay->GetBinContent(j) > 0)
                _significance += pow(significance(hmc.GetBinContent(j), histoverlay->GetBinContent(j)),2);
            }
          }
          if(doROC && sensitivevariable == name[i]){
            double bkgintegral = ROC_bkg->Integral();
            double sigintegral = ROC_sig->Integral();
            double sigeff = 1;
            double bkgrej = 0;
            ROC->SetPoint(0,sigeff,bkgrej);
            for (int i = 1; i < ROC_sig->GetNbinsX()+1; ++i)
            {
              sigeff -= ROC_sig->GetBinContent(i)/sigintegral;
              bkgrej += ROC_bkg->GetBinContent(i)/bkgintegral;
              ROC->SetPoint(i,sigeff,bkgrej);
            }
            outputrocfile->cd();
            ROC->Write(overlaysample + "_ROC");
            ROC_sig->Write(overlaysample + "_ROC_sig");
            ROC_bkg->Write(overlaysample + "_ROC_bkg");
            deletepointer(ROC);
            deletepointer(ROC_sig);
            deletepointer(ROC_bkg);
          }
          std::string samptitle = histoverlay->GetTitle();
          findAndReplaceAll(samptitle," ","~");
          if(samptitle.find("#") != string::npos) samptitle = "$"+samptitle+"$";
          findAndReplaceAll(samptitle,"#","\\");
          findAndReplaceAll(samptitle,"%","\\%");
          findAndReplaceAll(samptitle,"rightarrow","to ");
          yield_chart->set(samptitle,regtitle,integral(histoverlay));
          sgnf_chart->set(samptitle,regtitle,sqrt(_significance));
          printf("signal %s yield: %4.2f, significance: %4.2f\n",overlaysample.Data(), histoverlay->Integral(), sqrt(_significance));
        }

        if(ratio > 0) histoverlay->Scale(ratio);
        if(overlaysample != "") histoverlay->Draw("hist same");
        lgsig->SetBorderSize(0);
        lgsig->Draw();
        padhi->Update();
        cv.SaveAs("plots_" + outdir + "/" + region + "/" + name[i] + ".pdf");
        deletepointer(histoverlay);
        deletepointer(lgsig);
      }
      deletepointer(hsk);
      deletepointer(lg1);
      deletepointer(padlow );
      deletepointer(padhi  );
      deletepointer(datahist);
      for(auto &iter : buffer) deletepointer(iter);
      if(debug) printf("end region %s\n",region.Data());
      cv.SaveAs("plots_" + outdir + "/" + region + "/" + name[i] + ".pdf]");
      cv.Clear();
    }
    if(debug) printf("end loop region\n");
  }
  if(yield_chart->rows.size()){
    yield_chart->caption = "The sample and data yield before the fit.";
    yield_chart->print("yield_chart");
    sgnf_chart->caption = "The stat. only significance of the signal in each regions with the benchmark $\\mu$ value.";
    sgnf_chart->print("significance_chart");
  }
  deletepointer(yield_chart);
  deletepointer(sgnf_chart);
  outputrocfile->Close();
  deletepointer(outputrocfile);
}
