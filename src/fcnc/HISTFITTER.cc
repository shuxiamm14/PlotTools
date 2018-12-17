#include "HISTFITTER.h"

TMinuit* gM = 0;

HISTFITTER::HISTFITTER(){}

void HISTFITTER::addfithist(TString component,  TH1D* inputhist, int begin, int end){
	if ( fithists.find(component) == fithists.end() )
	{
		fithists[component] = new TH1D(component, component, nregion, 0, nregion);
		iregion[component] = 0;
	}
	iregion[component]++;
	fithists[component]->SetBinContent(iregion[component], inputhist->Integral(begin,end));

	double error = 0;
	for (int ib = begin; ib < end+1; ++ib)
	{
		if (inputhist->GetBinContent(ib))
			error += pow(inputhist->GetBinError(ib),2);
	}
	error = sqrt(error);
	fithists[component]->SetBinError(iregion[component],error);
}

void HISTFITTER::fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag) {
	f = 0;
	map<TString, TH1D*>* histforfit = (map<TString, TH1D*>*) gM->GetObjectFit();
	if (!histforfit)
	{
	   printf("hist isn't found\n");
	   exit(1);
	}	
	map<TString, TH1D*>::iterator iter;
	iter = histforfit->begin();
	TH1D *htot = (TH1D*)iter->second->Clone("htot");
	TH1D *data;
	htot->Scale(par[0]);
	iter++;
	int i = 1;
	for (; iter != histforfit->end(); iter ++)
	{

		if(iter->first == "data") {
			data = iter->second;
			continue;
		}
		htot->Add(iter->second,i<4?par[i]:1);
		i++;
	}
	for (int i = 1; i <= data->GetNbinsX(); ++i)
		if(htot->GetBinContent(i))
			f += pow((data->GetBinContent(i) - htot->GetBinContent(i))/htot->GetBinError(i),2);
	deletepointer(htot);
}

double HISTFITTER::fit(double *bstvl, double *error){


	gM = new TMinuit(5);
	gM->SetFCN(fcn);
	gM->SetPrintLevel(-1);
	
	Double_t arglist[10];
	Int_t ierflg = 0;
	
	arglist[0] = 1;
	gM->mnexcm("SET ERR", arglist ,1,ierflg);

	gM->mnparm(0, "sf_b", 1, 0.1, 0.,2.,ierflg);
	gM->mnparm(1, "sf_c", 1, 0.1, 0.,2.,ierflg);
	gM->mnparm(2, "sf_g", 2, 0.1, 0.,3.,ierflg);
	gM->mnparm(3, "sf_j", 1, 0.1, 0.,2.,ierflg);
	gM->SetObjectFit((TObject*)&fithists);
   
    arglist[0] = 4;
    arglist[1] = 60.;
    Double_t val[4],err[4];
   
    gM->mnexcm("SCAN", arglist ,2,ierflg);
    for (int i = 0; i < 4; ++i) gM->GetParameter(i,val[i],err[i]);
	gM->mnparm(0, "sf_b", val[0], 0.1, 0.,2.,ierflg);
	gM->mnparm(1, "sf_c", val[1], 0.1, 0.,2.,ierflg);
	gM->mnparm(2, "sf_g", val[2], 0.1, 0.,10.,ierflg);
	gM->mnparm(3, "sf_j", val[3], 0.1, 0.,2.,ierflg);
   
    arglist[0] = 1000;
    arglist[1] = 0;
	gM->mnexcm("MIGRADE", arglist ,2,ierflg);
	for (int i = 0; i < 4; ++i) gM->GetParameter(i,bstvl[i],error[i]);

	bstvl = val;
	error = err;
	Double_t minf;
	Double_t  	fedm;
	Double_t  	errdef;
	Int_t	npari;
	Int_t	nparx;
	Int_t	istat;
	gM->mnstat(minf,fedm,errdef,npari,nparx,istat);
	return minf;

}

void HISTFITTER::debug(){
	TFile debugfile("debugfile","recreate");
	for (iter = fithists.begin(); iter != fithists.end(); iter ++){
		iter->second->Write();
	}
}

void HISTFITTER::clear(){

	for (iter = fithists.begin(); iter != fithists.end(); iter ++){
		deletepointer(iter->second);
	}
	fithists.clear();
}