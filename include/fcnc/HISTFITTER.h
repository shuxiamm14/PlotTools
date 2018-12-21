#include "fcnc_include.h"

class HISTFITTER
{
public:
	HISTFITTER();
	~HISTFITTER();

	map<TString, TH1D*> fithists;
	int nparam;
	map<TString, int> iregion;
	map<TString, TH1D*>::iterator iter;
	TH1D *htot = NULL;
	int nregion = 100;
	TString paramname[100];
	double startpoint[100];
	double stepsize[100];
	double lowrange[100];
	double highrange[100];
	void asimovfit(int fitnumber, TString outfile);
	void addfithist(TString component, TH1D* inputhist, int begin, int end);
	static void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);
	double fit(double *bstvl, double *error, bool asimov);
	double setparam(TString _paramname, double _startpoint, double _stepsize, double _lowrange, double _highrange);
	void debug();
	void clear();

};