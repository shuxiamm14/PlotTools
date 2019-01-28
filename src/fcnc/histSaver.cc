#include "histSaver.h"
#include "TGaxis.h"
histSaver::histSaver() {
  nvar = 0;
  for(Int_t i=0; i<50; i++) {
    nbin[i] = 1; xlo[i] = 0; xhi[i] = 1; var1[i] = 0; var2[i] = 0; MeVtoGeV[i] = 0; var3[i] = 0;
  }
}

histSaver::~histSaver() {
    deletepointer(inputfile);
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


void histSaver::add(const char* titleX_, const char* name_, const char* unit_) {
  fromntuple = 0;
  if(nvar>=0 && nvar<50) {
    titleX[nvar] = titleX_;
    name[nvar] = name_;
    unit[nvar] = unit_;
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
  for(auto const& region: regions) {
    printf("histSaver::show()\t region: %s\n", region.Data());
  }
  for (int i = 0; i < nvar; ++i)
  {
    printf("histSaver::show()\t%s = ", name[i].Data());
    if(var2[i]) printf("%d\n", *var2[i]);
    else if(var1[i]) printf("%f\n", MeVtoGeV[i] ? *var1[i]/1000 : *var1[i]);
    else if(var3[i]) printf("%f\n", MeVtoGeV[i] ? *var3[i]/1000 : *var3[i]);
  }
}

float histSaver::binwidth(int i){
  return (xhi[i]-xlo[i])/nbin[i];
}

void histSaver::init_sample(TString samplename, TString histname, TString sampleTitle, enum EColor color){

  current_sample = samplename;

  if(plot_lib.find(samplename) != plot_lib.end()) return;

  if(debug) printf("reading hist file: %s\n", (histfilename + ".root").Data());
  if(fromntuple){
    inputfile = new TFile(histfilename + ".root", "recreate");
  }else{
    inputfile = new TFile(histfilename + ".root", "read");
  }
  
  if(debug) printf("add new sample: %s\n", samplename.Data());

  for(auto const& region: regions) {
    for (int i = 0; i < nvar; ++i){
      if(inputfile->Get(histname + "_" + region + "_" + name[i]))
        plot_lib[samplename][region].push_back((TH1D*)inputfile->Get(histname + "_" + region + "_" + name[i]));
      else{
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
    }
    if(debug == 1) printf("plot_lib[%s][%s]\n", samplename.Data(), region.Data());
  }
  if (samplename == "data") dataref = 1;

  if(debug) printf("finished initializing %s\n", samplename.Data() );
}

void histSaver::read_sample(TString samplename, TString histname, TString sampleTitle, enum EColor color, double norm){

  if(!inputfile) inputfile = new TFile(histfilename + ".root", "read");

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
        if(debug == 1) {
          printf("histogram name: %s\n", (histname+"_"+region+"_"+name[i]).Data());
          printf("plot_lib[%s][%s][%d]\n", samplename.Data(), region.Data(), i);
        }
        plot_lib[samplename][region][i]->Add((TH1D*)inputfile->Get(histname+"_"+region+"_"+name[i]),norm);
      }
    }else{
      for (int i = 0; i < nvar; ++i){
        if(debug == 1) {
          printf("histogram name: %s\n", (histname+"_"+region+"_"+name[i]).Data());
          printf("plot_lib[%s][%s][%d]\n", samplename.Data(), region.Data(), i);
        }
        plot_lib[samplename][region].push_back((TH1D*)(inputfile->Get(histname+"_"+region+"_"+name[i])->Clone(histname+"_"+region+"_"+name[i]+ "_" + samplename)));
        plot_lib[samplename][region][i]->Scale(norm);
        plot_lib[samplename][region][i]->SetTitle(sampleTitle);
        plot_lib[samplename][region][i]->SetFillColor(color);
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
  for(auto const& region: regions) {
    for (int i = 0; i < nvar; ++i){
      map<TString, map<TString, vector<TH1D*>>>::iterator iter;
      for(iter=plot_lib.begin(); iter!=plot_lib.end(); iter++){
        inputfile->cd();
        iter->second[region][i]->Write("",TObject::kWriteDelete);
      }
    }
  }
}

void histSaver::clearhist(){
  for(auto const& region: regions) {
    for (int i = 0; i < nvar; ++i){
      map<TString, map<TString, vector<TH1D*>>>::iterator iter;
      for(iter=plot_lib.begin(); iter!=plot_lib.end(); iter++){
        iter->second[region][i]->Reset();
      }
    }
  }
}

void histSaver::overlay(TString _overlaysample){
  overlaysample = _overlaysample;
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

      TPad *padlow = new TPad("lowpad","lowpad",0,0,1,0.3);
      TPad *padhi  = new TPad("hipad","hipad",0,0.3,1,1);
      TH1D hmc("hmc","hmc",nbin[i]/irebin,xlo[i],xhi[i]);
      TH1D hmcR("hmcR","hmcR",nbin[i]/irebin,xlo[i],xhi[i]);
      TH1D hdataR("hdataR","hdataR",nbin[i]/irebin,xlo[i],xhi[i]);

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
      map<TString, map<TString, vector<TH1D*>>>::iterator iter;
      TH1D *histoverlay;
      for(iter=plot_lib.begin(); iter!=plot_lib.end(); iter++){
        if(irebin != 1) iter->second[region][i]->Rebin(irebin);
        if(iter->first == "data") continue;
        if(iter->first == overlaysample){
          histoverlay = iter->second[region][i];
          continue;
        }
        if(debug) {
          printf("plot_lib[%s][%s][%d]\n", iter->first.Data(), region.Data(), i);
        }
        hsk->Add(iter->second[region][i]);
        hmc.Add(iter->second[region][i]);
        lg1->AddEntry(iter->second[region][i],iter->second[region][i]->GetTitle(),"F");
      }
      if(overlaysample != ""){
        if(debug) { printf("overlay: %s\n", overlaysample.Data()); }
        lg1->AddEntry(histoverlay,histoverlay->GetTitle(),"LP");
        histoverlay->SetLineStyle(9);
        histoverlay->SetLineWidth(3);
        histoverlay->SetLineColor(kRed);
        histoverlay->SetFillColor(0);
      }
      //hsk->GetXaxis()->SetTitle(unit[i] == "" ? titleX[i].Data() : (titleX[i] + " [" + unit[i] + "]").Data());
      if (dataref) {
        lg1->AddEntry(plot_lib["data"][region][i],"data","LP");
        plot_lib["data"][region][i]->GetXaxis()->SetTitle(unit[i] == "" ? titleX[i].Data() : (titleX[i] + " [" + unit[i] + "]").Data());
        plot_lib["data"][region][i]->GetXaxis()->SetLabelColor(kWhite);
        //plot_lib["data"][region][i]->SetMaximum(1.8*plot_lib["data"][region][i]->GetMaximum());
        char str[30];
        sprintf(str,"Events / %4.2f %s",binwidth(i)*irebin, unit[i].Data());
        plot_lib["data"][region][i]->GetYaxis()->SetTitle(str);
        plot_lib["data"][region][i]->GetYaxis()->SetTitleOffset(1.2);
        plot_lib["data"][region][i]->SetMarkerStyle(20);
        plot_lib["data"][region][i]->SetMarkerSize(0.4);
        plot_lib["data"][region][i]->Draw("E0 same");
        SetMax(hsk,plot_lib["data"][region][i],1.6);
        plot_lib["data"][region][i]->SetMinimum(0);
      }else{
        hsk->SetMaximum(1.6*hsk->GetMaximum());
        hsk->SetMinimum(0);
      }

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
      ATLASLabel(0.2,0.900,"work in progress",kBlack, region);
      hsk->Draw("hist same");
      if(overlaysample != "") histoverlay->Draw("hist same");
      hmc.Draw("E2,same");
      if(dataref) plot_lib["data"][region][i]->Draw("E+ same");
      lg1->Draw("same");

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
      hdataR.GetYaxis()->SetTitleOffset(hdataR.GetYaxis()->GetTitleOffset()*1.05);
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
