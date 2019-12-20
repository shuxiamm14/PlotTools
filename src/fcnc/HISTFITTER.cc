#include "HISTFITTER.h"
using namespace std;
TMinuit* gM = 0;

HISTFITTER::HISTFITTER(){
	nparam = 0;
	htot = NULL;
	debug = 1;
	nregion = 100;
}
HISTFITTER::~HISTFITTER(){}
void HISTFITTER::addfithist(TString component,  TH1D* inputhist, int begin, int end, TString fitparam){
	TString componentwofit = component;
	bool fit = 0;
	if(fitparam != ""){
		for (int i = 0; i < nparam; ++i)
		{
			if(paramname[i] == fitparam) {
				component += ("_fit" + to_string(i)).c_str();
				fit = 1;
			}
		}
	}
	if(fit){
		int nfit = 1;
		for(auto comp : fithists){
			if(comp.first.Contains("_fit") && split(comp.first.Data(),"_fit")[0] == componentwofit){
				iregion[comp.first]++;
				nfit = iregion[comp.first];
			}
		}
		if(iregion[component]==0) iregion[component] = nfit;
	}else
		iregion[component]++;
	if ( fithists.find(component) == fithists.end() )
	{
		fithists[component] = new TH1D(component, component, nregion, 0, nregion);
	}
	if(fit){
		int nfit = 0;
		for(auto comp : fithists){
			if(comp.first.Contains("_fit") && split(comp.first.Data(),"_fit")[0] == componentwofit){
				iregion[comp.first]++;
				nfit = iregion[comp.first];
			}
		}
		if(iregion[component]==0) iregion[component] = nfit;
	}else
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
void HISTFITTER::calculateEigen(){

	double covariance_matrix[100];
	gM->mnemat(&covariance_matrix[100],nparam);
	if(debug){
		for (int i = 0; i < nparam; ++i){
			printf("covariance matrix: ");
			for (int j = 0; j < nparam; ++j)
				printf(" %f", covariance_matrix[j*nparam + i]);
			printf("\n");
		}
	}
	float **covariance_matrix2 = new float*[nparam];
	for (int i = 0; i < nparam; ++i)
	{
		covariance_matrix2[i] = new float[nparam];
		for (int j = 0; j < nparam; ++j)
		{
			covariance_matrix2[i][j] = covariance_matrix[i*nparam + j];
		}
	}
	eigenval = new float[nparam];
	eigenvector = new float*[nparam];
	for (int i = 0; i < nparam; ++i)
	{
		eigenvector[i] = new float[nparam];
	}
	EigenVectorCalc(covariance_matrix2,nparam,eigenval,eigenvector);
	if(debug){
		printf("eigen values: ");
		for (int i = 0; i < nparam; ++i)
		{
			printf(" %f", eigenval[i]);
		}
		printf("\n");
		for (int i = 0; i < nparam; ++i){
			printf("eigenvectors: ");
			for (int j = 0; j < nparam; ++j)
				printf(" %f", eigenvector[j][i]);
			printf("\n");
		}
	}
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
	
	TH1D *htot = NULL;
	TH1D *data;
	for (iter = histforfit->begin(); iter != histforfit->end(); iter ++)
	{
		if(iter->first == "metadata") continue;
		if(histforfit->find("asimovdata")!=histforfit->end()){
			if(iter->first == "data") {
				continue;
			}
			if (iter->first == "asimovdata")
			{
				data = iter->second;
				continue;
			}
		}else{
			if (iter->first == "data")
			{
				data = iter->second;
				continue;
			}
		}
		int ipar = parsecomponentname(iter->first);
		double histscale = ipar >=0? par[ipar] : 1;
		if (!htot)
		{
			htot = (TH1D*)iter->second->Clone("htotfcn");
			htot->Scale(histscale);
		}else
			htot->Add(iter->second,histscale);
	}

	if(!htot) {
		printf("HISTFITTER::fcn() : ERROR: htot is empty, please check if histforfit has histograms:\n");
		for (auto hist: *histforfit)
		{
			printf(" %s ", hist.first.Data());
		}
		exit(0);
	}

	for (int i = 1; i <= data->GetNbinsX(); ++i){
		if(htot->GetBinContent(i))
			f += pow((data->GetBinContent(i) - htot->GetBinContent(i))/htot->GetBinError(i),2);
	}
	deletepointer(htot);
}

int HISTFITTER::parsecomponentname(TString name){
	if(name.Contains("fit")){
		return int(stof(split(name.Data(),"fit")[1].Data()));
	}
	return -1;
}

void HISTFITTER::asimovfit(int fitnumber, TString outfile){
	TTree *asimovFitResults = new TTree("asimovFitResults","asimovFitResults");
	double *val, *err;
	val = (double*)malloc(nparam*sizeof(double));
	err = (double*)malloc(nparam*sizeof(double));

	for (int i = 0; i < nparam; ++i)
		asimovFitResults->Branch(paramname[i],&(val[i]));
	double chi2;
	asimovFitResults->Branch("Chi2", &chi2);
	for (int i = 0; i < fitnumber; ++i)
	{
		chi2 = fit(val,err,1);
		asimovFitResults->Fill();
	}
	TFile *fitresultfile = new TFile(outfile,"recreate");
	asimovFitResults->Write();
	fitresultfile->Close();
	deletepointer(fitresultfile);
}
void HISTFITTER::setparam(TString _paramname, double _startpoint, double _stepsize, double _lowrange, double _highrange){
	printf("set parameter: %s\n", _paramname.Data());
	paramname.push_back(_paramname);
	startpoint.push_back(_startpoint);
	stepsize.push_back(_stepsize);
	lowrange.push_back(_lowrange);
	highrange.push_back(_highrange);
	nparam++;
}
double HISTFITTER::fit(double *bstvl, double *error, bool asimov){

	gRandom->SetSeed(0);
	gM = new TMinuit(5);
	gM->SetFCN(fcn);
	gM->SetPrintLevel(-1);
	
	Double_t arglist[10];
	Int_t ierflg = 0;
	
	arglist[0] = 1;
	gM->mnexcm("SET ERR", arglist ,1,ierflg);

	for (int i = 0; i < nparam; ++i)
		gM->mnparm(i, paramname[i], startpoint[i], stepsize[i], lowrange[i], highrange[i],ierflg);

	double stat = 1;
	if (asimov)
	{
		if(htot == NULL)
			for (iter = fithists.begin(); iter != fithists.end(); iter ++){
				if(iter->first == "data" || iter->first == "asimovdata") continue;
				if(htot == NULL) {
					htot = (TH1D*)iter->second->Clone("htot");
				}
				else{
					htot->Add(iter->second);
				}
			}
		fithists["asimovdata"] = (TH1D*)(fithists["data"]->Clone("asimovdata"));
		fithists["asimovdata"] ->Reset();
		double tmpintegral = htot->Integral()*stat;
		double tmpweight = 1./stat;
		for (int i = 0; i < tmpintegral; ++i)
		{
			fithists["asimovdata"] -> Fill(htot->GetRandom(),tmpweight);
		}
	}else if(fithists.find("asimovdata")!=fithists.end()){
		deletepointer(fithists["asimovdata"]);
		fithists.erase("asimovdata");
	}
	//fithists["metadata"] = new TH1D("metadata","metadata",100,0,100);
	//savemetadata(h_metadata, "nparam",nparam);
	gM->SetObjectFit((TObject*)&fithists);
   
    arglist[0] = nparam;	//number of scan dimentions
    arglist[1] = 60.; //number of scan points ,maximum 100
    Double_t val[10],err[10];
   
    gM->mnexcm("SCAN", arglist ,2,ierflg);
    for (int i = 0; i < nparam; ++i) gM->GetParameter(i,val[i],err[i]);
	for (int i = 0; i < nparam; ++i)
		gM->mnparm(i, paramname[i], val[i], stepsize[i], lowrange[i], highrange[i],ierflg);
   
    arglist[0] = 1000; //max calls
    arglist[1] = 0.1;	//tolerance
	gM->mnexcm("MIGRADE", arglist ,2,ierflg);
	for (int i = 0; i < nparam; ++i) gM->GetParameter(i,bstvl[i],error[i]);
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

void HISTFITTER::savemetadata(TH1D *metadatahist, TString what, double value){
	TAxis *xaxis = metadatahist->GetXaxis();
	int i = 1;
	while (i <= metadatahist->GetNbinsX() && TString(xaxis->GetBinLabel(i))!="" )
	{
		if(what == xaxis->GetBinLabel(i)) {
			metadatahist->SetBinContent(i,value);
			return;
		}
		i++;
	}
	metadatahist->SetBinContent(i,value);
	metadatahist->GetXaxis()->SetBinLabel(i,what);
}

double HISTFITTER::readmetadata(TH1D *metadatahist, TString what){
	TAxis *xaxis = metadatahist->GetXaxis();
	int i = 1;
	while (i <= metadatahist->GetNbinsX() && TString(xaxis->GetBinLabel(i))!="" )
	{
		if(what == xaxis->GetBinLabel(i)) {
			return metadatahist->GetBinContent(i);
		}
		i++;
	}
	printf("HISTFITTER::readmetadata() : ERROR : metadata %s not found\n", what.Data());
	exit(0);
}

void HISTFITTER::debugfile(){
	TFile debugfile("debugfile","recreate");
	for (iter = fithists.begin(); iter != fithists.end(); iter ++){
		iter->second->Write();
	}
}

void HISTFITTER::clear(){

	for (iter = fithists.begin(); iter != fithists.end(); iter ++){
		deletepointer(iter->second);
	}
	deletepointer(htot);
	iregion.clear();
	fithists.clear();
}
