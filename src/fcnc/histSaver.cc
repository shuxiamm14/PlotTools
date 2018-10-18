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
  for (int i = 0; i < nvar; ++i)
  {
    printf("histSaver::show()\t%s = ", name[i].Data());
    if(var2[i]) printf("%d\n", *var2[i]);
    else printf("%f\n", MeVtoGeV[i] ? *var1[i]/1000 : *var1[i]);
  }
}

float histSaver::binwidth(int i){
  return (xhi[i]-xlo[i])/nbin[i];
}

void histSaver::init_sample(TString samplename, TString histname, TString sampleTitle, enum EColor color){
  current_sample = samplename;
  if(plot_lib.find(samplename) != plot_lib.end()) return;
  if(debug) printf("add new sample: %s\n", samplename.Data());
  vector<TH1D*> plots;
  for(auto const& region: regions) {
    for (int i = 0; i < nvar; ++i){
      plots.push_back(new TH1D(histname + "_" + region + "_" + name[i].Data(),sampleTitle,nbin[i],xlo[i],xhi[i]));
      plots[i]->Sumw2();
      plots[i]->SetFillColor(color);
      plots[i]->SetLineWidth(0.9);
      plots[i]->SetLineColor(kBlack);
      plots[i]->SetMarkerSize(0);
    }
    if(debug == 1) printf("plot_lib[%s][%s]\n", samplename.Data(), region.Data());
    plot_lib[samplename][region] = plots;
    plots.clear();
  }
  if (samplename == "data") dataref = 1;

  if(debug) printf("finished initializing %s\n", samplename.Data() );
}

void histSaver::add_region(TString region){
  regions.push_back(region);
}

void histSaver::fill_hist(TString sample, TString region){
  if (weight_type == 0)
  {
    printf("ERROR: weight not set\n");
  }
  for (int i = 0; i < nvar; ++i){
    if(debug == 1) printf("plot_lib[%s][%s][%d]->Fill(%f,%f)\n", sample.Data(), region.Data(), i, getVal(i), weight_type == 1? *fweight : *dweight);
    plot_lib[sample][region][i]->Fill(getVal(i),weight_type == 1? *fweight : *dweight);
  }
}

void histSaver::fill_hist(TString region){
  fill_hist(current_sample,region);
}

void histSaver::fill_hist(){
  fill_hist("nominal");
}

void histSaver::plot_stack(){
  SetAtlasStyle();
  for(auto const& region: regions) {
    gSystem->mkdir(region);
    gSystem->mkdir(region + "/eps");
    gSystem->mkdir(region + "/root");
    for (int i = 0; i < nvar; ++i){

      TCanvas cv("cv","cv",600,600);

      TFile savehist(region + "/root/" + name[i] + ".root","recreate");

      TPad *padlow = new TPad("lowpad","lowpad",0,0,1,0.3);
      TPad *padhi  = new TPad("hipad","hipad",0,0.3,1,1);
      TH1D hmc("hmc","hmc",nbin[i],xlo[i],xhi[i]);
      TH1D *hmcR   = new TH1D("hmcR","hmcR",nbin[i],xlo[i],xhi[i]);
      TH1D *hdataR = new TH1D("hdataR","hdataR",nbin[i],xlo[i],xhi[i]);

      cv.cd();

//===============================upper pad===============================
      padhi->SetBottomMargin(0.015);
      padhi->cd();
      hmc.Sumw2();
      THStack *hsk = new THStack(name[i].Data(),name[i].Data());
      map<TString, map<TString, vector<TH1D*>>>::iterator iter;
      TLegend* lg1 = 0;
      lg1 = new TLegend(0.43,0.75,0.94,0.90,"");
      lg1->SetNColumns(2);
      for(iter=plot_lib.begin(); iter!=plot_lib.end(); iter++){
        savehist.cd();
        iter->second[region][i]->Write();
        if(iter->first == "data") continue;
        hsk->Add(iter->second[region][i]);
        hmc.Add(iter->second[region][i]);
        lg1->AddEntry(iter->second[region][i],iter->second[region][i]->GetTitle(),"F");
      }
      //hsk->GetXaxis()->SetTitle(unit[i] == "" ? titleX[i].Data() : (titleX[i] + " [" + unit[i] + "]").Data());
      if (dataref) {
        lg1->AddEntry(plot_lib["data"][region][i],"data","LP");
        plot_lib["data"][region][i]->GetXaxis()->SetTitle(unit[i] == "" ? titleX[i].Data() : (titleX[i] + " [" + unit[i] + "]").Data());
        plot_lib["data"][region][i]->GetXaxis()->SetLabelColor(kWhite);
        //plot_lib["data"][region][i]->SetMaximum(1.8*plot_lib["data"][region][i]->GetMaximum());
        char str[30];
        sprintf(str,"Events / %4.2f %s",binwidth(i), unit[i].Data());
        plot_lib["data"][region][i]->GetYaxis()->SetTitle(str);
        plot_lib["data"][region][i]->SetMarkerStyle(20);
        plot_lib["data"][region][i]->SetMarkerSize(0.8);
        plot_lib["data"][region][i]->Draw("E1 same");
      }else{
        hsk->SetMaximum(1.8*hsk->GetMaximum());
      }
      lg1->Draw("same");


      hmc.SetFillColor(1);
      hmc.SetLineColor(0);
      hmc.SetMarkerStyle(1);
      hmc.SetMarkerSize(0);
      hmc.SetMarkerColor(1);
      hmc.SetFillStyle(3004);

      hsk->Draw("hist same");
      hmc.Draw("E2,same");
      if(dataref) plot_lib["data"][region][i]->Draw("E+ same");
      cv.cd();
      padhi->Draw();
//===============================lower pad===============================
      padlow->SetFillStyle(4000);
      padlow->SetGrid(1,1);
      padlow->SetTopMargin(0.03);
      padlow->SetBottomMargin(0.35);
      padlow->cd();

      for(Int_t j=1; j<nbin[i]+1; j++) {
        hmcR->SetBinContent(j,1);
        hmcR->SetBinError(j,hmc->GetBinContent(j)>0 ? hmc->GetBinError(j)/hmc->GetBinContent(j) : 0);
        hdataR->SetBinContent(j, hmc->GetBinContent(j)>0 ? plot_lib["data"][region][i]->GetBinContent(j)/hmc->GetBinContent(j) : 1);
        hdataR->SetBinError(j, ( plot_lib["data"][region][i]->GetBinContent(j)>0 && hmc->GetBinContent(j)>0 )? plot_lib["data"][region][i]->GetBinError(j)/hmc->GetBinContent(j) : 0);
      }

      hdataR->SetMarkerStyle(20);
      hdataR->SetMarkerSize(0.8);
      hdataR->SetMaximum(1.5);
      hdataR->SetMinimum(0.5);
      hdataR->GetYaxis()->SetNdivisions(504,false);
      hdataR->GetYaxis()->SetTitle("Data/Bkg");
      hdataR->GetYaxis()->SetTitleOffset(hdataR->GetYaxis()->GetTitleOffset()*1.1);
      hdataR->GetYaxis()->CenterTitle();
      hdataR->GetXaxis()->SetTitle(unit[i] == "" ? titleX[i].Data() : (titleX[i] + " [" + unit[i] + "]").Data());
      hdataR->GetXaxis()->SetTitleSize(hdataR->GetXaxis()->GetTitleSize()*0.7);
      hdataR->GetYaxis()->SetTitleSize(hdataR->GetYaxis()->GetTitleSize()*0.7);
      hmcR->SetFillColor(1);
      hmcR->SetLineColor(0);
      hmcR->SetMarkerStyle(1);
      hmcR->SetMarkerSize(0);
      hmcR->SetMarkerColor(1);
      hmcR->SetFillStyle(3004);
      hdataR->Draw("E same");
      hdataR->GetXaxis()->SetTitleOffset(3.4);
      hdataR->GetXaxis()->SetLabelSize(hdataR->GetXaxis()->GetLabelSize()*0.7); 
      hdataR->GetYaxis()->SetLabelSize(hdataR->GetYaxis()->GetLabelSize()*0.7); 
      hmcR->Draw("E2 same");
      TLine *line = new TLine();
      line->SetLineColor(2);
      line->DrawLine(hdataR->GetBinLowEdge(1), 1., hdataR->GetBinLowEdge(hdataR->GetNbinsX()+1), 1.);
      cv.cd();
      padlow->Draw();
      cv.SaveAs((CharAppend(region + "/eps/", name[i]) + ".eps"));
      savehist.Close();
      deletepointer(hsk);
      deletepointer(lg1);
      deletepointer(padlow );
      deletepointer(padhi  );
      deletepointer(hmcR   );
      deletepointer(hdataR );
      deletepointer(line   );
    }
  }
}