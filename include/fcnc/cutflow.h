#ifndef CUTFLOW
#define CUTFLOW

#include <vector>
#include <iostream>
#include "TString.h"
class CutFlow
{
public:
	CutFlow(TString _sample = "", TString _region = "");
	~CutFlow();

	Float_t* fweight;
	Double_t* dweight;
	ULong64_t* event_number;
	TString sample;
	TString region;
	std::vector<std::vector<ULong64_t>> event_track;
	int n_tracked_event;
	int n_cuts;
	int i_cut;
	bool if_track;
	int weight_type;
	int ievt_track;
	std::vector<long> cutflow_raw;
	std::vector<double> cutflow_weighted;
	std::vector<double> cutflow2;
	std::vector<TString> cut_names;
	void trackEvent(ULong64_t evt_number);
	void setEventNumber(ULong64_t* _event_number){ event_number = _event_number; }
	void setWeight(Float_t* _weight){ fweight = _weight; weight_type = 1;}
	void setWeight(Double_t* _weight){ dweight = _weight; weight_type = 2;}
	void clear();
	void newEvent();
	void print();
	void save(int total_cuts);
	void fill(TString cut_name = "");
};
#endif