#include "../header/sm_para.h"
#include "../header/method.h"
#include "../header/syst_summary.h"

int syst() {

  gErrorIgnoreLevel = kError;
  TGaxis::SetMaxDigits(4);
  //gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  // loop over input files
  string line;
  int f_indx = 0;

  double nega_err2_sum = 0.;
  double plus_err2_sum = 0.;
  
  ifstream filelist(input_folder + "/input.txt"); 

  cout << "input_folder: " << input_folder << endl;
  //cout << para_label << endl;
  
  if (filelist.is_open()) {
    while (!filelist.eof()) {
      if (getline(filelist, line, '\n')) {
	if (line[0] != '!') {
	  TFile* f_in = new TFile(line.data());
	  TString fname(line.data());

	  cout << "line.data = " << line.data() << endl;

	  //getObj(f_in);
	  //cout << line.data() << endl;

	  TTree *tree_in = (TTree*) f_in -> Get("TRESULT");

	  int err_type = 0;
	  double nega_err = 0.;
	  double plus_err = 0.;
	    
	  int entries = tree_in -> GetEntries();
	  for (Int_t irow = 0; irow < entries; irow++) {// tree-loop
	    
	    tree_in -> GetEntry(irow);
	    
	    nega_err = tree_in -> GetLeaf("Br_SYST_ERR") -> GetValue(0);
	    plus_err = tree_in -> GetLeaf("Br_SYST_ERR") -> GetValue(1);
	    
	  }// end tree-loop

	  err_type = tree_in -> GetLeaf("Br_err_type") -> GetValue(0);

	  nega_err2_sum += nega_err * nega_err;
	  plus_err2_sum += plus_err * plus_err;
	  
	  cout << "err_type = " << err_type << "\n"
	       << "nega_err2_sum = " << nega_err2_sum << ", plus_err2_sum = " << plus_err2_sum << "\n"
	       << "\tsyst. err = (" << nega_err << ", " << plus_err << ")\n\n";
	  
	  f_indx ++;
	  
	}
      }
    }
  }

  double syst_nega = 0., syst_plus = 0.;
  syst_nega = -TMath::Sqrt(nega_err2_sum);
  syst_plus = TMath::Sqrt(plus_err2_sum);
  
  cout << "\nTotla Syst. Err = (" << syst_nega << ", " << syst_plus << ")" << endl;

  return 0;
  
}
