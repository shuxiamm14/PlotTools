Author: Boyang Li

email: boyang.li@cern.ch

Project git: https://gitlab.cern.ch/boyang/plotTools

//=====================================installation=====================================

need external package: eigen3 (http://eigen.tuxfamily.org/index.php?title=Main_Page)
in CMakelist.txt: set(EIGEN3DIR <PATH TO EIGEN>)

. env.sh
mkdir build
cd build
cmake ..
make install

add ./include/fcnc ./include/atlasstyle ./include/external into your Makefile include flag
add ./lib into your Makefile link lib flag

//=====================================Usage1: histSaver=====================================

//==========Plot stacks with n-tuples=========
#include "histSaver.h"

// sequence matters : variable -> region -> samples
tau_plots = new histSaver();
tau_plots->set_weight(&weight);
tau_plots->debug = 0/1;
tau_plots->histfile = "outputhistograms.root";
//seq: nbins, xlow, xhigh, hist title, hist name , variable address, if convert from MeV to GeV, unit
tau_plots->add(10,25.,125.,"p_{T,#tau}","taupt",&tau_pt_0,true,"GeV");
tau_plots->add(10,25.,125.,"p_{T,b}","bpt",&pt_b,true,"GeV");
tau_plots->add(10,25.,125.,"p_{T,light-jet}","ljetpt",&pt_ljet,true,"GeV");

tau_plots->add_region("the regions you have 1");
tau_plots->add_region("the regions you have 2");
tau_plots->add_region("the regions you have 3");
tau_plots->add_region("the regions you have 4");

//seq: sample name, sample fill name (in case memory leak), sample title name, sample color in the stack.
tau_plots->init_sample("data","data","data",kBlack);
tau_plots->init_sample("other","other","Other samples",kYellow);
tau_plots->init_sample("ttbar_g","ttbar_g","t#bar{t}(gluon fake #tau)",(enum EColor)7);
tau_plots->init_sample("ttbar_j","ttbar_j","t#bar{t}(light-jet fake #tau)",kBlue);
tau_plots->init_sample("ttbar_b","ttbar_b","t#bar{t}(b-jets fake #tau)",kViolet);
tau_plots->init_sample("ttbar_lep","ttbar_lep","t#bar{t}(lepton fake #tau)",kGreen);
tau_plots->init_sample("ttbar_real","ttbar_real","t#bar{t}(real #tau)",kRed);
tau_plots->init_sample("ttbar_c","ttbar_c","t#bar{t}(c-jets fake #tau)",kOrange);
tau_plots->init_sample("ttbar_nomatch","ttbar_nomatch","t#bar{t}(no truth matched fake #tau)",kGray);

tau_plots->init_sample("signal1","signal1","signal1",kRed);
tau_plots->init_sample("signal2","signal2","signal2",kRed);
tau_plots->init_sample("signal2","signal2","signal2",kRed);
//plot_lib[sample][region][ivar]
if("you wanna see the list of input variable") tau_plots->show();

for (Long64_t jentry=0; jentry<nentries;jentry++) {
	fChain->GetEntry(jentry);

	if(sampleisGluon && region1cut) tau_plots->fill_hist("ttbar_g","the regions you have 1");
		//...
}
tau_plots->stackorder.push_back("ttbar_g")
tau_plots->stackorder.push_back("ttbar_j")
tau_plots->stackorder.push_back("ttbar_b")
...

tau_plots->overlay("signal1")
tau_plots->overlay("signal2")
tau_plots->overlay("signal3")

tau_plots->write(TFile* outputfile) // write the histograms into rootfile for further use
tau_plots->plot_stack();

//=============Read from a histogram============
tau_plots = new histSaver();

tau_plots->histfile = "yourhistfile.root" //histogram in the file with name: samplefillname_regionname_varname

tau_plots->add("p_{T,#tau}","taupt","GeV",nrebin1);
tau_plots->add("p_{T,b}","bpt","GeV",nrebin2);
tau_plots->add("p_{T,light-jet}","ljetpt","GeV",nrebin3);

tau_plots->add_region("the regions you have 1");
tau_plots->add_region("the regions you have 2");
tau_plots->add_region("the regions you have 3");
tau_plots->add_region("the regions you have 4");
//seq: sample name, sample fill name (in case memory leak), sample title name, sample color in the stack.
tau_plots->read_sample("data","data","data",kBlack);
tau_plots->read_sample("other","other","Other samples",kYellow);
tau_plots->read_sample("ttbar_g","ttbar_g","t#bar{t}(gluon fake #tau)",(enum EColor)7);
tau_plots->read_sample("ttbar_j","ttbar_j","t#bar{t}(light-jet fake #tau)",kBlue);
tau_plots->read_sample("ttbar_b","ttbar_b","t#bar{t}(b-jets fake #tau)",kViolet);
tau_plots->read_sample("ttbar_lep","ttbar_lep","t#bar{t}(lepton fake #tau)",kGreen);
tau_plots->read_sample("ttbar_real","ttbar_real","t#bar{t}(real #tau)",kRed);
tau_plots->read_sample("ttbar_c","ttbar_c","t#bar{t}(c-jets fake #tau)",kOrange);
tau_plots->read_sample("ttbar_nomatch","ttbar_nomatch","t#bar{t}(no truth matched fake #tau)",kGray);

tau_plots->stackorder.push_back("ttbar_g")
tau_plots->stackorder.push_back("ttbar_j")
tau_plots->stackorder.push_back("ttbar_b")
...

tau_plots->overlay("signal1")
tau_plots->overlay("signal2")
tau_plots->overlay("signal3")

tau_plots->plot_stack();

//================features===========
void muteregion(TString keyword);			
void unmuteregion(TString keyword);			//decide if the region contains the keyword is plotted
void overlay(TString _overlaysample);		//the sample _overlaysample is shown as overlay in the plots instead of stack
TH1D* grabhist(TString sample, TString region, int ivar);	//get specific histogram
void merge_regions(TString inputregion1, TString inputregion2, TString outputregion); //merge 2 regions into another region and keeps the inputs

double templatesample(TString fromregion,string formula,TString toregion,TString newsamplename,TString newsampletitle,enum EColor color,bool scaletogap, double SF = 1); //used for fake estimation methods.
//example: tau_plots->templatesample("ss_region","1 data -1 smhiggs -1 wjet -1 diboson -1 zll -1 ztautau -1 top -1 fake","os_region","fakeSS","Fake",kYellow,0,1.31597);

void write_trexinput(TString NPname = "NOMINAL", TString writeoption = "recreate"); //in case you are using TRexFitter (https://gitlab.cern.ch/TRExStats/TRExFitter) This function will generate the histogram inputs.


//=====================================Usage2: HISTFITTER=====================================
//fit the data with norm of different component (1 bin per addfithist call)
//give the fitresult and calculate the eigen vector, eigen value of covariance matrix
#include "HISTFITTER.h"
HISTFITTER* fitter = new HISTFITTER();
//setparam and addfit should follow the same order
fitter->setparam("sf_b", 1, 0.1, 0.,2.);
fitter->setparam("sf_c", 1, 0.1, 0.,2.);
fitter->setparam("sf_g", 2, 0.1, 0.,10.);
fitter->setparam("sf_j", 1, 0.1, 0.,2.);
fitter->addfithist("data",hist1,binlow,binhigh);
fitter->addfithist("data",hist1,binlow,binhigh);
fitter->addfithist("b",hist1,binlow,binhigh);
fitter->addfithist("b",hist2,binlow,binhigh);
fitter->addfithist("c",hist1,binlow,binhigh);
fitter->addfithist("c",hist2,binlow,binhigh);
fitter->addfithist("g",hist1,binlow,binhigh);
fitter->addfithist("g",hist2,binlow,binhigh);
fitter->addfithist("j",hist1,binlow,binhigh);
fitter->addfithist("j",hist2,binlow,binhigh);
double chi2 = fitter->fit(val,err,0);
printf("%s, ptbin: %d, b: %f+/-%f, c: %f+/-%f, g: %f+/-%f, j: %f+/-%f;  Chi2:%f\n",nprong[iprong].Data(), ptbin+1, val[0],err[0], val[1],err[1], val[2],err[2], val[3],err[3], chi2);
fitter->calculateEigen();
fitter->clear();

//=====================================Usage3: EigenVector=====================================
//calculator: Calculate the eigen vector and eigen value for a given matrix. Run ./bin/test_run (util/test_run) to see how to use
#include "EigenVectorCalc.h"
void EigenVectorCalc(float **matrix, int matrixsize, float *eigenval, float **eigenvectors)


//=====================================Usage4: cutflow=====================================
#include "cutflow.h"

cutflow mycut;

mycut->trackevent(event_Number1);
mycut->trackevent(event_Number2);
mycut->trackevent(event_Number3);
mycut->trackevent(event_Number4);
...

mycut.set_weight(&weight);

for (int i = 0; i < entries; ++i)
{
	tree->GetEntry(i);
	mycut->newEvent();
	if(failcut1) continue;
	mycut->fill();
	if(failcut2) continue;
	mycut->fill();
	...
}
mycut.print();