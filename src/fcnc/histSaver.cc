#include "histSaver.h"

histSaver::histSaver() {
  nvar = 0;
  for(Int_t i=0; i<50; i++) {
    nbin[i] = 1; xlo[i] = 0; xhi[i] = 1; var1[i] = 0; var2[i] = 0; MeVtoGeV[i] = 0;
  }
}

histSaver::~histSaver() {}

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
    if(var1[i]) tmp = MeVtoGeV[i] ? *var1[i]/1000 : *var1[i];
    else if(var2[i]) tmp = *var2[i];
  }
  if (!ifRebin[i]){
    if(tmp >= xhi[i]) tmp = xhi[i]*0.999999;
    if(tmp < xlo[i]) tmp = xlo[i];
  }else{
    if(tmp >= xbins[i][nbin[i]]) tmp = xbins[i][nbin[i]]*0.999999;
    if(tmp < xbins[i][0]) tmp = xbins[i][0];
  }
  return tmp;
}

void histSaver::show(){
  for (int i = 0; i < nvar; ++i)
  {
    printf("histSaver::show()\t%s = ", name[i].c_str());
    if(var2[i]) printf("%d\n", *var2[i]);
    else printf("%f\n", MeVtoGeV[i] ? *var1[i]/1000 : *var1[i]);
  }
}

float histSaver::binwidth(int i){
  return (xhi[i]-xlo[i])/nbin[i];
}

void histSaver::init_sample(TString samplename, TString sampleTitle, enum EColor color){
  current_sample = samplename;
  if(plot_lib.find(samplename) != plot_lib.end()) return;
  printf("add new sample: %s\n", samplename.Data());
  vector<TH1D*> plots;
  for (int i = 0; i < nvar; ++i){
    plots.push_back(new TH1D(samplename + "_" + name[i].c_str(),sampleTitle,nbin[i],xlo[i],xhi[i]));
    plots[i]->Sumw2();
    plots[i]->SetLineColor(kBlack);
    plots[i]->SetFillColor(color);
    plots[i]->SetLineWidth(0.3);
  }
  for(auto const& region: regions) {
    plot_lib[samplename][region] = plots;
  }
  if (samplename == "data") dataref = 1;

  printf("finished initializing %s\n", samplename.Data() );
}

void histSaver::add_region(TString region){
  regions.push_back(region);
}

void histSaver::fill_hist(TString sample, TString region){
  if (weight_type == 0)
  {
    printf("ERROR: weight not set\n");
  }
  for (int i = 0; i < nvar; ++i)
    plot_lib[sample][region][i]->Fill(getVal(i),weight_type == 1? *fweight : *dweight);
}

void histSaver::fill_hist(TString region){
  fill_hist(current_sample,region);
}

void histSaver::fill_hist(){
  fill_hist("nominal");
}

void histSaver::plot_stack(){
  for(auto const& region: regions) {
    gSystem->mkdir(region);
    for (int i = 0; i < nvar; ++i){
      TCanvas cv;
      THStack *hsk = new THStack(name[i].c_str(),name[i].c_str());
      map<TString, map<TString, vector<TH1D*>>>::iterator iter;
      TLegend* lg1;
      lg1 = new TLegend(0.53,0.75,0.94,0.90,"");
      lg1->SetNColumns(2);
      for(iter=plot_lib.begin(); iter!=plot_lib.end(); iter++){
        if(iter->first == "data") continue;
        hsk->Add(iter->second[region][i]);
        lg1->AddEntry(iter->second[region][i],iter->second[region][i]->GetTitle(),"F");
      }

      hsk->SetMaximum(1.4*hsk->GetMaximum());

      if (dataref) plot_lib["data"][region][i]->Draw("E");
      hsk->Draw("hist same");
      hsk->GetXaxis()->SetTitle(unit[i] == "" ? titleX[i].c_str() : (titleX[i] + " [" + unit[i] + "]").Data());
      char str[30];
      sprintf(str,"Events / %4.2f %s",binwidth(i), unit[i].Data());
      hsk->GetYaxis()->SetTitle(str);
      lg1->Draw("same");
      cv.SaveAs((CharAppend(region + "/", name[i]) + ".eps"));
      deletepointer(hsk);
      deletepointer(lg1);
    }
  }
}