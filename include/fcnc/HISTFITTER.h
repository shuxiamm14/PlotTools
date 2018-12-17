#include "fcnc_include.h"

class HISTFITTER
{
public:
	HISTFITTER();
	~HISTFITTER();

	map<TString, TH1D*> fithists;
	map<TString, int> iregion;
	map<TString, TH1D*>::iterator iter;
	int nregion = 100;
	void addfithist(TString component, TH1D* inputhist, int begin, int end);
	static void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);
	double fit(double *bstvl, double *error);
	void debug();
	void clear();

};