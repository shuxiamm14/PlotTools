#include "cutflow.h"
#include "TFile.h"
#include "TH1D.h"
#include "TAxis.h"
#include "fcnc_include.h"
CutFlow::CutFlow(TString _sample, TString _region) :
sample(_sample), region(_region),
n_tracked_event(0), n_cuts(0), i_cut(0), if_track(0), weight_type(0), ievt_track(0),
fweight(NULL), dweight(NULL), event_number(NULL)
{
}


CutFlow::~CutFlow(){
}

void CutFlow::clear(){
	for(auto &iter : cutflow_raw) iter = 0;
	for(auto &iter : cutflow_weighted) iter = 0;
	for(auto &iter : cutflow2) iter = 0;
	if(n_tracked_event) {
		std::vector<ULong64_t> tmp = event_track[0];
		event_track.clear();
		event_track[0] = tmp;
	}
	i_cut = 0;
}

void CutFlow::newEvent(){
	if_track = 0;
	if(n_tracked_event)
		for(auto evt: event_track[0]){
			if(*event_number == evt)
				if_track = 1;
		}
	i_cut = 0;
}

void CutFlow::trackEvent(ULong64_t evtnumber) {
	n_tracked_event++;
	if(n_tracked_event == 1)
		event_track.push_back(std::vector<ULong64_t>());
	event_track[0].push_back(evtnumber);

}


void CutFlow::fill(TString cut_name){
	double weight = weight_type == 1? *fweight : *dweight;
	if(cutflow_raw.size() <= i_cut) {
		cut_names.push_back(cut_name);
		cutflow_raw.push_back(1);
		cutflow_weighted.push_back(weight);
		cutflow2.push_back(weight*weight);
		n_cuts ++;
	}else{
		cutflow_raw[i_cut] += 1;
		cutflow_weighted[i_cut] += weight;
		cutflow2[i_cut] += weight*weight;
	}
	if(if_track){
		if(event_track.size() <= i_cut+1) event_track.push_back(std::vector<ULong64_t>());
		event_track[i_cut+1].push_back(*event_number);
	}
	i_cut ++;
}

void CutFlow::save(int total_cuts){
	if(n_cuts == 0) {
		printf("CutFlow::save() : WARNING : no Cut applied, nothing saved\n");
		return;
	}
	int nbins(n_cuts);
	if(total_cuts > n_cuts) nbins = total_cuts;
	TFile *save_file = new TFile("cutflow_" + sample + ".root", "update");
	TH1D *save_hist = (TH1D*)save_file->Get(region);
	if(!save_hist) save_hist = new TH1D(region,region,nbins,0,nbins);
	TAxis *xaxis = save_hist->GetXaxis();
	for (int i = 1; i <= nbins; ++i)
	{
		if(xaxis->GetBinLabel(i) == cut_names[0] || xaxis->GetBinLabel(i) == TString("")){
			for (int j = 0; j < n_cuts; ++j)
			{
				save_hist->SetBinContent(i+j,cutflow_weighted[j]);
				save_hist->SetBinError(i+j,sqrt(cutflow2[j]));
				xaxis->SetBinLabel(i+j,cut_names[j]);
			}
			break;
		}
	}
	save_file->cd();
	save_hist->Write(region,TObject::kWriteDelete);
	save_file->Close();
	deletepointer(save_file);
}

void CutFlow::print(){
	printf("cutflow:");
	if(n_cuts == 0) {
		printf("CutFlow::print() : no Cut applied\n");
		return;
	}

	for (int i = 0; i < n_cuts; ++i)
	{
		printf(" %s, ", cut_names[i].Data());
	}
	printf("\n");
	for (int i = 0; i < n_cuts; ++i)
	{
		printf(" %f+/-%f", cutflow_weighted[i], sqrt(cutflow2[i]));
	}
	printf("\n");
	printf("cutflow_raw:");
	for (int i = 0; i < n_cuts; ++i)
	{
		printf(" %ld", cutflow_raw[i]);
	}
	printf("\n");
	if(n_tracked_event){
		for(auto cut: event_track){
			printf("cut: ");
			for(auto evt : cut)
				printf("%llu,",evt);
			printf("\n");
		}
	}
}
