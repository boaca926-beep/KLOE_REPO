#include "syst_plot.h"

//int syst_plot(const TString para_name = "Mass_omega") {
int syst_plot() {

  gErrorIgnoreLevel = kError;
  TGaxis::SetMaxDigits(3);
  //gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);  

  //cout << "Plot syst. for choice of 3pi cross section models; BW vs VMD" << endl;

  // define variables
  double YLIST[1000], YLIST_Err[1000];
  
  // loop over all input root files
  const TString syst_path = "./path_list";

  cout << "\nROOT input file list: " << syst_path << "\n"
       << para_name << "\n";
  
  string line;
  int file_indx = 0;

  ifstream filelist(syst_path); 

  if (filelist.is_open()) {
    while (!filelist.eof()) {
      if (getline(filelist, line, '\n')) {
	if (line[0] != '!') {

	  // Input file
	  TFile* f_input = new TFile(line.data());
	  TString f_name(line.data());

	  cout << file_indx + 1 << ": " << f_name << endl;

	  // Tree
	  /*
	  TTree* TCRX3PI = (TTree*)f_input -> Get("TCRX3PI");
  
	  TIter next_tree(f_input -> GetListOfKeys());
	  TString objnm_tree, classnm_tree;

	  // Loop over entries
	  int i = 0;
	  
	  while ( (key = (TKey *) next_tree() ) ) {
	    
	    objnm_tree   = key -> GetName();
	    classnm_tree = key -> GetClassName();
	    key -> GetSeekKey();

	    cout << " tree" << i << ": classnm = " << classnm_tree << ", objnm = " << objnm_tree << endl;

	  }
	  */

	  TTree *TCRX3PI = (TTree*)f_input -> Get("TCRX3PI");
	  
	  for (Int_t irow = 0; irow < TCRX3PI -> GetEntries(); irow++) {// begin for loop
	    
	    TCRX3PI -> GetEntry(irow);

	  }
	  
	  YLIST[file_indx] = TCRX3PI -> GetLeaf("Br_" + para_name + "_fit") -> GetValue(0);
	  YLIST_Err[file_indx] = TCRX3PI -> GetLeaf("Br_" + para_name + "_err_fit") -> GetValue(0);

	  //cout << file_indx + 1 << ": " << YLIST[file_indx] << " +/- " << YLIST_Err[file_indx] << endl;
	    


	  file_indx ++;

	  f_input -> Close();

	}
      }
    }
  }

  
  //cout << "\n" << endl;
  
  // calculate errors
  double y_tmp = YLIST[0], y_err_tmp = YLIST_Err[0];
  double y_diff = 0., y_diff_abs = 0.;
  double y_uncorr_err = 0.; // uncorrelated error
  double z = 0.; // significance

  for (int i = 0; i < file_indx; i ++) {

    y_diff = YLIST[i] - y_tmp;
    y_diff_abs = TMath::Abs(y_diff);
    
    y_uncorr_err = TMath::Sqrt(TMath::Abs(YLIST_Err[i] * YLIST_Err[i] - y_err_tmp * y_err_tmp));

    if (y_err_tmp != YLIST_Err[i]) z = y_diff / TMath::Sqrt(TMath::Abs(YLIST_Err[i] * YLIST_Err[i] - y_err_tmp * y_err_tmp));
    
    
    cout << i + 1 << ": y_a = " << YLIST[i] << " +/- " << YLIST_Err[i] << ", y_b = " << y_tmp << " +/- " << y_err_tmp << ", y_diff = " << y_diff << ", |y_diff| / 2 = " << y_diff_abs / 2.  << ", y_uncorr_err = " << y_uncorr_err << ", y_diff / y_uncorr_err = " << z << endl;

    y_tmp = YLIST[i];
    y_err_tmp = YLIST_Err[i];
    
  }
  
  // output root file
  TFile *f_output = new TFile("./syst_output.root", "recreate");

  TTree* TPara = new TTree("TPara_" + para_name, "recreate");
  TPara -> SetAutoSave(0);

  TString s = Form("%i", file_indx);
  
  TPara -> Branch("Br_value", &YLIST, "Br_value[" + s + "]/D");
  TPara -> Branch("Br_err", &YLIST_Err, "Br_err[" + s + "]/D");

  TPara -> Fill();
  TPara -> Write();

  f_output -> Close();
  
  return 0;
  
}
