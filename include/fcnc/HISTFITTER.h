#include "fcnc_include.h"
#include "EigenVectorCalc.h"

class HISTFITTER
{
public:
	HISTFITTER();
	~HISTFITTER();

	int nparam;
	TH1D *htot;
	bool debug;
	int nregion;

	std::map<TString, TH1D*> fithists;
	std::map<TString, int> iregion;
	std::map<TString, TH1D*>::iterator iter;
	float *eigenval;
	float **eigenvector;
	TH1D *h_metadata;
	std::vector<TString> paramname;
	std::vector<double > startpoint;
	std::vector<double > stepsize;
	std::vector<double > lowrange;
	std::vector<double > highrange;
	void asimovfit(int fitnumber, TString outfile);
	void addfithist(TString component, TH1D* inputhist, int begin, int end, TString fitparam = "");
	static int parsecomponentname(TString name);
	static void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);
	static void savemetadata(TH1D *metadatahist, TString what, double value);
	static double readmetadata(TH1D *metadatahist, TString what);
	double fit(double *bstvl, double *error, bool asimov);
	void setparam(TString _paramname, double _startpoint, double _stepsize, double _lowrange, double _highrange);
	void debugfile();
	void clear();
	void calculateEigen();

};