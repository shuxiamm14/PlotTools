//Usage:

//=====================================Plot stacks with n-tuples=====================================

// sequence matters : variable -> region -> samples

tau_plots = new histSaver();
tau_plots->set_weight(&weight);
tau_plots->debug = 0/1;

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

//plot_lib[sample][region][ivar]
if("you wanna see the list of input variable") tau_plots->show();

for (Long64_t jentry=0; jentry<nentries;jentry++) {
	fChain->GetEntry(jentry);

	if(sampleisGluon && region1cut) tau_plots->fill_hist("ttbar_g","the regions you have 2");
		//...
}

tau_plots->plot_stack();

//=====================================Read from a histogram=====================================
// structure: regionname/root/varname.root : samplefillname_regionname_varname

tau_plots = new histSaver();

tau_plots->add("p_{T,#tau}","taupt","GeV");
tau_plots->add("p_{T,b}","bpt","GeV");
tau_plots->add("p_{T,light-jet}","ljetpt","GeV");

tau_plots->add_region("the regions you have 1");
tau_plots->add_region("the regions you have 2");
tau_plots->add_region("the regions you have 3");
tau_plots->add_region("the regions you have 4");

tau_plots->read_sample("data","data","data",kBlack);
tau_plots->read_sample("other","other","Other samples",kYellow);
tau_plots->read_sample("ttbar_g","ttbar_g","t#bar{t}(gluon fake #tau)",(enum EColor)7);
tau_plots->read_sample("ttbar_j","ttbar_j","t#bar{t}(light-jet fake #tau)",kBlue);
tau_plots->read_sample("ttbar_b","ttbar_b","t#bar{t}(b-jets fake #tau)",kViolet);
tau_plots->read_sample("ttbar_lep","ttbar_lep","t#bar{t}(lepton fake #tau)",kGreen);
tau_plots->read_sample("ttbar_real","ttbar_real","t#bar{t}(real #tau)",kRed);
tau_plots->read_sample("ttbar_c","ttbar_c","t#bar{t}(c-jets fake #tau)",kOrange);
tau_plots->read_sample("ttbar_nomatch","ttbar_nomatch","t#bar{t}(no truth matched fake #tau)",kGray);
