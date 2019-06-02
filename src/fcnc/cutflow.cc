#include "cutflow.h"

cutflow::cutflow(){
	nCuts = 0;
	iCut = 0;
	fweight = NULL;
	dweight = NULL;

}


cutflow::~cutflow(){
}

void cutflow::clear(){
	for(auto &iter : cutflowraw) iter = 0;
	for(auto &iter : cutflowweighted) iter = 0;
	for(auto &iter : cutflow2) iter = 0;
	iCut = 0;
}

void cutflow::newEvent(){
	iCut = 0;
}

void cutflow::fill(){
	double weight = weight_type == 1? *fweight : *dweight;
	if(cutflowraw.size() <= iCut) {
		cutflowraw.push_back(1);
		cutflowweighted.push_back(weight);
		cutflow2.push_back(weight*weight);
		nCuts ++;
	}else{
		cutflowraw[iCut] += 1;
		cutflowweighted[iCut] += weight;
		cutflow2[iCut] += weight*weight;
	}
	iCut ++;
}

void cutflow::print(){
	printf("cutflow:");
	for (int i = 0; i < nCuts-1; ++i)
	{
		printf(" %f +/- %f,", cutflowweighted[i], sqrt(cutflow2[i]));
	}
	printf(" %f +/- %f\n", cutflowweighted[nCuts-1], sqrt(cutflow2[nCuts-1]));
	printf("cutflowraw:");
	for (int i = 0; i < nCuts-1; ++i)
	{
		printf(" %ld,", cutflowraw[i]);
	}
	printf(" %ld\n", cutflowraw[nCuts-1]);
}