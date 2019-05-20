#include <vector>
#include <iostream>
#include "TString.h"
class cutflow
{
public:
	cutflow();
	~cutflow();

	Float_t* fweight;
	Double_t* dweight;
	int weight_type;

	std::vector<long> cutflowraw;
	std::vector<double> cutflowweighted;
	std::vector<double> cutflow2;
	int nCuts;
	int iCut;
	std::vector<TString> cutname;
	void set_weight(Float_t* _weight){ fweight = _weight; weight_type = 1;}
	void set_weight(Double_t* _weight){ dweight = _weight; weight_type = 2;}
	void clear();
	void newEvent();
	void print();
	void fill();
};