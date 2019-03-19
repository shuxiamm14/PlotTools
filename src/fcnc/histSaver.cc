#include "histSaver.h"
#include "TGaxis.h"
histSaver::histSaver(TString outputfilename) {
  outputfile = new TFile (outputfilename + ".root", "recreate");
  nvar = 0;
  inputfilename = "hists";
  nregion = 0;
  blinding = 0;
  fweight = NULL;
  dweight = NULL;
  weight_type = 0;
  overlaysample = "";
  inputfile = 0;
  lumi = "#it{#sqrt{s}} = 13TeV, 80 fb^{-1}";
  analysis = "FCNC tqH H#rightarrow tautau";
  workflow = "work in progress";
  fromntuple = 1;
  histcount = 0;
  this_region = "nominal";
  read_path = "./" ;
  debug = 1;
  for(Int_t i=0; i<50; i++) {
    nbin[i] = 1; xlo[i] = 0; xhi[i] = 1; var1[i] = 0; var2[i] = 0; MeVtoGeV[i] = 0; var3[i] = 0;
  }
}

histSaver::~histSaver() {
    deletepointer(inputfile);
}

TH1D* histSaver::grabhist(TString sample, TString region, int ivar){
  if(plot_lib.find(sample) == plot_lib.end()){
    printf("histSaver:grabhist  ERROR: sample %s not found\n", sample.Data());
    show();
    exit(1);
  }

  if(plot_lib[sample].find(region) == plot_lib[sample].end()){
    printf("histSaver:grabhist  ERROR: region %s for sample %s not found\n", region.Data(), sample.Data());
    show();
    exit(1);
  }
  if(!plot_lib[sample][region][ivar]) printf("histSaver:grabhist  WARNING: empty histogram%s\n", name[ivar].Data());
  return plot_lib[sample][region][ivar];
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
    name[nvar] = name_;
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
    if(debug == 1) printf("fill value: %f\n", tmp);
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
  for(auto iter:plot_lib ){
    printf("histSaver::show()\tsample: %s\n", iter.first.Data());
  }
  for(auto const& region: regions) {
    printf("histSaver::show()\tregion: %s\n", region.Data());
  }
  for (int i = 0; i < nvar; ++i)
  {
    if(var2[i]) printf("histSaver::show()\t%s = %d\n", name[i].Data(), *var2[i]);
    else if(var1[i]) printf("histSaver::show()\t%s = %f\n", name[i].Data(), MeVtoGeV[i] ? *var1[i]/1000 : *var1[i]);
    else if(var3[i]) printf("histSaver::show()\t%s = %f\n", name[i].Data(), MeVtoGeV[i] ? *var3[i]/1000 : *var3[i]);
  }
}

float histSaver::binwidth(int i){
  return (xhi[i]-xlo[i])/nbin[i];
}

void histSaver::merge_regions(TString inputregion1, TString inputregion2, TString outputregion){
  if(debug) printf("histSaver::merge_regions\t %s and %s into %s\n",inputregion1.Data(),inputregion2.Data(),outputregion.Data());
  bool exist = 0;
  for(auto iter:plot_lib ){
    if(iter.second.find(inputregion1) == iter.second.end()){
      printf("histSaver::merge_regions\t inputregion1: %s not found",inputregion1.Data());
      show();
      exit(1);
    }
    if(iter.second.find(inputregion2) == iter.second.end()){
      printf("histSaver::merge_regions\t inputregion2: %s not found",inputregion2.Data());
      show();
      exit(1);
    }
    bool outputexist = 0;
    if(iter.second.find(outputregion) != iter.second.end()){
      printf("histSaver::merge_regions\t outputregion %s exist, overwrite it\n",outputregion.Data());
      exist = 1;
      for (int i = 0; i < nvar; ++i) deletepointer(iter.second[outputregion][i]);
      iter.second[outputregion].clear();
    }
    for (int i = 0; i < nvar; ++i)
    {
      iter.second[outputregion].push_back((TH1D*)iter.second[inputregion1][i]->Clone(iter.first+"_"+outputregion+"_"+name[i]));
      if(debug)
        printf("add %s to %s as %s\n", iter.second[inputregion2][i]->GetName(),iter.second[inputregion1][i]->GetName(),(iter.first+"_"+outputregion+"_"+name[i]).Data());
      iter.second[outputregion][i]->Add(iter.second[inputregion2][i]);
    }
  }
  if(!exist) regions.push_back(outputregion);
}

void histSaver::init_sample(TString samplename, TString histname, TString sampleTitle, enum EColor color){

  outputfile->cd();
  current_sample = samplename;

  if(plot_lib.find(samplename) != plot_lib.end()) return;
  
  if(debug) printf("add new sample: %s\n", samplename.Data());

  for(auto const& region: regions) {
    for (int i = 0; i < nvar; ++i){
      plot_lib[samplename][region].push_back(new TH1D(histname + "_" + region + "_" + name[i],sampleTitle,nbin[i],xlo[i],xhi[i]));
      if (samplename != "data")
      {
        plot_lib[samplename][region][i]->Sumw2();
        plot_lib[samplename][region][i]->SetFillColor(color);
        plot_lib[samplename][region][i]->SetLineWidth(1);
        plot_lib[samplename][region][i]->SetLineColor(kBlack);
        plot_lib[samplename][region][i]->SetMarkerSize(0);
      }
    }
    if(debug == 1) printf("plot_lib[%s][%s]\n", samplename.Data(), region.Data());
  }
  if (samplename == "data") dataref = 1;

  if(debug) printf("finished initializing %s\n", samplename.Data() );
}

void histSaver::read_sample(TString samplename, TString histname, TString sampleTitle, enum EColor color, double norm){

  if(!inputfile) inputfile = new TFile(inputfilename + ".root", "read");

  if (samplename == "data") dataref = 1;
  bool newsample = plot_lib.find(samplename) == plot_lib.end();
  for(auto const& region: regions) {
    if (debug == 1)
    {
      printf("read sample %s from %s region\n", samplename.Data(), region.Data());
    }
    ++histcount;
    if (!newsample)
    {
      for (int i = 0; i < nvar; ++i)
      {
        if(!(TH1D*)(inputfile->Get(histname+"_"+region+"_"+name[i]))) {
          printf("histogram name not found: %s\n", (histname+"_"+region+"_"+name[i]).Data());
          printf("plot_lib[%s][%s][%d]\n", samplename.Data(), region.Data(), i);
          show();
          exit(1);
        }
        plot_lib[samplename][region][i]->Add((TH1D*)inputfile->Get(histname+"_"+region+"_"+name[i]),norm);
      }
    }else{
      for (int i = 0; i < nvar; ++i){
        if(!(TH1D*)(inputfile->Get(histname+"_"+region+"_"+name[i]))) {
          printf("histogram name not found: %s\n", (histname+"_"+region+"_"+name[i]).Data());
          printf("plot_lib[%s][%s][%d]\n", samplename.Data(), region.Data(), i);
          show();
          exit(1);
        }
        plot_lib[samplename][region].push_back((TH1D*)(inputfile->Get(histname+"_"+region+"_"+name[i])->Clone(histname+"_"+region+"_"+name[i])));
        plot_lib[samplename][region][i]->SetName(samplename+"_"+region+"_"+name[i]);
        plot_lib[samplename][region][i]->Scale(norm);
        plot_lib[samplename][region][i]->SetTitle(sampleTitle);
        plot_lib[samplename][region][i]->SetFillColorAlpha(color,1);
        plot_lib[samplename][region][i]->SetLineWidth(1);
        plot_lib[samplename][region][i]->SetLineColor(kBlack);
        plot_lib[samplename][region][i]->SetMarkerSize(0);
        if(histcount == 1){
          nbin[i] = plot_lib[samplename][region][i]->GetNbinsX();
          xlo[i] = plot_lib[samplename][region][i]->GetXaxis()->GetXmin();
          xhi[i] = plot_lib[samplename][region][i]->GetXaxis()->GetXmax();
        }
      }
    }
  }
}

void histSaver::add_region(TString region){
  regions.push_back(region);
  nregion += 1;
}

void histSaver::fill_hist(TString sample, TString region){
  if (weight_type == 0)
  {
    printf("ERROR: weight not set\n");
  }
  if(plot_lib.find(sample) == plot_lib.end()){
    printf("ERROR: sample %s not found while filling\n",sample.Data());
    show();
    exit(1);
  }
  if(plot_lib[sample].find(region) == plot_lib[sample].end()){
    printf("ERROR: region %s for sample %s not found while filling\n",region.Data(),sample.Data());
    show();
    exit(1);
  }
  for (int i = 0; i < nvar; ++i){
    if(debug == 1) printf("plot_lib[%s][%s][%d]->Fill(%f,%f)\n", sample.Data(), region.Data(), i, getVal(i), weight_type == 1? *fweight : *dweight);
    if(plot_lib[sample][region][i]) plot_lib[sample][region][i]->Fill(getVal(i),weight_type == 1? *fweight : *dweight);
    else printf("ERROR: histogram doesn't exist: plot_lib[%s][%s][%d]->Fill(%f,%f)\n", sample.Data(), region.Data(), i, getVal(i), weight_type == 1? *fweight : *dweight);
  }
}

void histSaver::fill_hist(TString region){
  fill_hist(current_sample,region);
}

void histSaver::fill_hist(){
  fill_hist("nominal");
}

void histSaver::write(){
  if(!outputfile) {
    printf("histSaver::write Error: outputfile pointer is empty\n");
    exit(1);
  }
  if(debug) printf("histSaver::write() Write to file: %s\n", outputfile->GetName());
  for(auto const& region: regions) {
    for (int i = 0; i < nvar; ++i){
      for(auto iter : plot_lib){
        outputfile->cd();
        if(debug) printf("histSaver::write() Check if region exist: %s\n", region.Data());
        if(iter.second.find(region) == iter.second.end()){
          show();
          printf("histSaver::write() Error: region not found in plot_lib,  sample: %s, variable: %s, region: %s\n",iter.first.Data(), name[i].Data(),region.Data());
          exit(1);
        }
        if(debug) printf("histSaver::write() Write Hist: plot_lib[%s][%s][%d]\n", iter.first.Data(),region.Data(),i);
        if(iter.second[region][i]) iter.second[region][i]->Write("",TObject::kWriteDelete);
        else{
          printf("histSaver::write() Error: histogram not found: sample: %s, variable: %s, region: %s\n",iter.first.Data(), name[i].Data(),region.Data());
        }
      }
    }
  }
}

void histSaver::clearhist(){
  if(debug) printf("histSaver::clearhist()\n");
  for(auto const& region: regions) {
    for (int i = 0; i < nvar; ++i){
      for(auto iter : plot_lib){
        if(iter.second[region][i]){
          iter.second[region][i]->Reset();
        }else{
          printf("histSaver::Reset() Error: histogram not found: sample: %s, variable: %s, region: %s\n",iter.first.Data(), name[i].Data(),region.Data());
        }
      }
    }
  }
}

void histSaver::overlay(TString _overlaysample){
  overlaysample = _overlaysample;
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
void histSaver::plot_stack(TString outputdir){
  SetAtlasStyle();
  TGaxis::SetMaxDigits(3);
  gSystem->mkdir(outputdir);
  int iregion = 0;
  TCanvas cv("cv","cv",600,600);
  for (int i = 0; i < nvar; ++i){
    cv.SaveAs(outputdir + "/" + name[i] + ".pdf[");
    for(auto const& region: regions) {
      bool muted = 0;
      for (auto const& mutedregion: mutedregions)
      {
        if(region.Contains(mutedregion))
          muted = 1;
      }
      if(muted) continue;
      TPad *padlow = new TPad("lowpad","lowpad",0,0,1,0.3);
      TPad *padhi  = new TPad("hipad","hipad",0,0.3,1,1);
      TH1D hmc("hmc","hmc",nbin[i]/rebin[i],xlo[i],xhi[i]);
      TH1D hmcR("hmcR","hmcR",nbin[i]/rebin[i],xlo[i],xhi[i]);
      TH1D hdataR("hdataR","hdataR",nbin[i]/rebin[i],xlo[i],xhi[i]);

      cv.cd();
      padhi->Draw();
//===============================upper pad===============================
      padhi->SetBottomMargin(0.017);
      padhi->cd();
      hmc.Sumw2();
      THStack *hsk = new THStack(name[i].Data(),name[i].Data());
      TLegend* lg1 = 0;
      lg1 = new TLegend(0.43,0.75,0.90,0.90,"");
      lg1->SetNColumns(2);
      TH1D *histoverlay;
      for(auto iter:plot_lib ){
        if(rebin[i] != 1) iter.second[region][i]->Rebin(rebin[i]);
        if(iter.first == "data") continue;
        if(iter.first == overlaysample){
          histoverlay = iter.second[region][i];
          continue;
        }
        if(debug) {
          printf("plot_lib[%s][%s][%d]\n", iter.first.Data(), region.Data(), i);
        }
        hsk->Add(iter.second[region][i]);
        hmc.Add(iter.second[region][i]);
        lg1->AddEntry(iter.second[region][i],iter.second[region][i]->GetTitle(),"F");
      }
      if(overlaysample != ""){
        if(debug) { printf("overlay: %s\n", overlaysample.Data()); }
        lg1->AddEntry(histoverlay,histoverlay->GetTitle(),"LP");
        histoverlay->SetLineStyle(9);
        histoverlay->SetLineWidth(3);
        histoverlay->SetLineColor(kRed);
        histoverlay->SetFillColor(0);
        SetMax(hsk,histoverlay,1);
      }
      if (dataref) {
        lg1->AddEntry(plot_lib["data"][region][i],"data","LP");
        plot_lib["data"][region][i]->SetMarkerStyle(20);
        plot_lib["data"][region][i]->SetMarkerSize(0.4);
        SetMax(hsk,plot_lib["data"][region][i],1.6);
        plot_lib["data"][region][i]->SetMinimum(0);
      }else{
        hsk->SetMaximum(1.6*hsk->GetMaximum());
        hsk->SetMinimum(0);
      }

      hsk->Draw("hist");
      hsk->GetXaxis()->SetTitle(unit[i] == "" ? titleX[i].Data() : (titleX[i] + " [" + unit[i] + "]").Data());
      hsk->GetXaxis()->SetLabelColor(kWhite);
      char str[30];
      sprintf(str,"Events / %4.2f %s",binwidth(i)*rebin[i], unit[i].Data());
      hsk->GetYaxis()->SetTitle(str);
      hsk->GetYaxis()->SetTitleOffset(1.4);
      hsk->GetYaxis()->SetLabelSize(hdataR.GetYaxis()->GetLabelSize()*0.7);
      hsk->GetXaxis()->SetLabelSize(hdataR.GetXaxis()->GetLabelSize()*0.7);
      hsk->GetYaxis()->SetTitleSize(hdataR.GetYaxis()->GetTitleSize()*0.7);
      if(blinding){
        for(Int_t j=1; j<nbin[i]+1; j++) {
          if(histoverlay->GetBinContent(j)/sqrt(plot_lib["data"][region][i]->GetBinContent(j)) > blinding) {
            plot_lib["data"][region][i]->SetBinContent(j,0);
            plot_lib["data"][region][i]->SetBinError(j,0);
          }
        }
      }

      for(Int_t j=1; j<nbin[i]+1; j++) {
        hmcR.SetBinContent(j,1);
        hmcR.SetBinError(j,hmc.GetBinContent(j)>0 ? hmc.GetBinError(j)/hmc.GetBinContent(j) : 0);
        hdataR.SetBinContent(j, hmc.GetBinContent(j)>0 ? plot_lib["data"][region][i]->GetBinContent(j)/hmc.GetBinContent(j) : 1);
        hdataR.SetBinError(j, ( plot_lib["data"][region][i]->GetBinContent(j)>0 && hmc.GetBinContent(j)>0 )? plot_lib["data"][region][i]->GetBinError(j)/hmc.GetBinContent(j) : 0);
      }

      hmc.SetFillColor(1);
      hmc.SetLineColor(0);
      hmc.SetMarkerStyle(1);
      hmc.SetMarkerSize(0);
      hmc.SetMarkerColor(1);
      hmc.SetFillStyle(3004);
      if(overlaysample != "") histoverlay->Draw("hist same");
      hmc.Draw("E2 same");
      if(dataref) plot_lib["data"][region][i]->Draw("E same");
      lg1->Draw("same");

      ATLASLabel(0.2,0.900,workflow.Data(),kBlack,lumi.Data(), analysis.Data(), region.Data());

//===============================lower pad===============================
      padlow->SetFillStyle(4000);
      padlow->SetGrid(1,1);
      padlow->SetTopMargin(0.03);
      padlow->SetBottomMargin(0.35);
      padlow->cd();


      hdataR.SetMarkerStyle(20);
      hdataR.SetMarkerSize(0.8);
      hdataR.SetMaximum(1.5);
      hdataR.SetMinimum(0.5);
      hdataR.GetYaxis()->SetNdivisions(504,false);
      hdataR.GetYaxis()->SetTitle("Data/Bkg");
      hdataR.GetYaxis()->SetTitleOffset(hdataR.GetYaxis()->GetTitleOffset()*1.06);
      hdataR.GetYaxis()->CenterTitle();
      hdataR.GetXaxis()->SetTitle(unit[i] == "" ? titleX[i].Data() : (titleX[i] + " [" + unit[i] + "]").Data());
      hdataR.GetXaxis()->SetTitleSize(hdataR.GetXaxis()->GetTitleSize()*0.7);
      hdataR.GetYaxis()->SetTitleSize(hdataR.GetYaxis()->GetTitleSize()*0.7);
      hmcR.SetFillColor(1);
      hmcR.SetLineColor(0);
      hmcR.SetMarkerStyle(1);
      hmcR.SetMarkerSize(0);
      hmcR.SetMarkerColor(1);
      hmcR.SetFillStyle(3004);
      hdataR.Draw("E same");
      hdataR.GetXaxis()->SetTitleOffset(3.4);
      hdataR.GetXaxis()->SetLabelSize(hdataR.GetXaxis()->GetLabelSize()*0.7); 
      hdataR.GetYaxis()->SetLabelSize(hdataR.GetYaxis()->GetLabelSize()*0.7); 
      hmcR.Draw("E2 same");
      TLine line;
      line.SetLineColor(2);
      line.DrawLine(hdataR.GetBinLowEdge(1), 1., hdataR.GetBinLowEdge(hdataR.GetNbinsX()+1), 1.);
      cv.cd();
      padlow->Draw();
      ++iregion;
      cv.SaveAs(outputdir + "/" + name[i] + ".pdf");
      deletepointer(hsk);
      deletepointer(lg1);
      deletepointer(padlow );
      deletepointer(padhi  );
    }
    cv.SaveAs(outputdir + "/" + name[i] + ".pdf]");
    cv.Clear();
  }
}
