#include "../header/path.h"
#include "../header/sm_para.h"
#include "../header/plot_sfw_syst.h"
#include "../header/method.h"

int plot_sfw_syst() {

  gErrorIgnoreLevel = kError;
  TGaxis::SetMaxDigits(3);
  //gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);  

  // loop over input files

  string line;
  int f_indx = 0;

  /*
  double w_y = 0.;
  double w_sum = 0.;
  double wy_sum = 0.;
  */
  
  //ifstream filelist(syst_path); 
  //ifstream filelist(syst_path + "sfw_output.txt"); 
  ifstream filelist(syst_path + pathType); 

  cout << syst_path + pathType << endl;
  
  if (filelist.is_open()) {
    while (!filelist.eof()) {
      if (getline(filelist, line, '\n')) {
	if (line[0] != '!') {

	  TFile* f_in = new TFile(line.data());
	  TString fname(line.data());
	  getObj(f_in);

	  double cut_value[1];
	  double x = 0., x_err = 0.;
	  /*
	  double isr3pi_sfw = 0., isr3pi_sfw_err = 0.;
	  double omegapi_sfw = 0., omegapi_sfw_err = 0.; 
	  double etagam_sfw = 0., etagam_sfw_err = 0.; 
	  double ksl_sfw = 0., ksl_sfw_err = 0.; 
	  double mcrest_sfw = 0., mcrest_sfw_err = 0.; 
	  */
	  
	  double SF[6]; // sfw2d MC scaling factors
	  double SF_ERR[6];
	  
	  TTree* tree_info = (TTree*)f_in -> Get("TRESULT");
	  tree_info -> SetBranchAddress("Br_SF", SF);
	  tree_info -> SetBranchAddress("Br_SF_ERR", SF_ERR);
	  tree_info -> SetBranchAddress("Br_cut_value", cut_value);

	  // Branches
	  int entries = tree_info -> GetEntries();
	  for (Int_t irow = 0; irow < entries; irow++) {// tree-loop
	    
	    tree_info -> GetEntry(irow);
	    //entries = tree_info -> GetLeaf("Br_entries") -> GetValue(0);
	    
	  }// end tree-loop

	  /*
	  isr3pi_sfw = SF[1];
	  isr3pi_sfw_err = SF_ERR[1];

	  omegapi_sfw = SF[2];
	  omegapi_sfw_err = SF_ERR[2];

	  etagam_sfw = SF[3];
	  etagam_sfw_err = SF_ERR[3];

	  ksl_sfw = SF[4];
	  ksl_sfw_err = SF_ERR[4];

	  mcrest_sfw = SF[5];
	  mcrest_sfw_err = SF_ERR[5];
	  */
	  
	  //<<"\tisr3pi_sfw = "<< isr3pi_sfw << "+/-" << isr3pi_sfw_err << "\n"
	  //<<"\tomegapi_sfw = "<< omegapi_sfw << "+/-" << omegapi_sfw_err << "\n"
	  //<<"\tetagam_sfw = "<< etagam_sfw << "+/-" << etagam_sfw_err << "\n"
	  //<<"\tksl_sfw = "<< ksl_sfw << "+/-" << ksl_sfw_err << "\n"
	  //<<"\tmcrest_sfw = "<< mcrest_sfw << "+/-" << mcrest_sfw_err << "\n";
	       
	  // Fill list
	  XLIST[f_indx] = cut_value[0];
	  XLIST_ERR[f_indx] = 0.;

	  YLIST[f_indx] = SF[para_indx];
	  YLIST_ERR[f_indx] = SF_ERR[para_indx];

	  cout << f_indx << ", " << cut_label << " = " << cut_value[0] << ", " << para_label << " = " << YLIST[f_indx] << "+/-" << YLIST_ERR[f_indx] << "\n";

	  /*
	  // normial value
	  */
	  
	  f_indx ++;
	  
	}
      }
    }
  }

  // Determine NORM
  int norm_indx = f_indx / 2;
  cout << "f_indx = " << f_indx << ", norm_indx = " << norm_indx << endl;

  getnorm(Y_NORM, Y_ERR_NORM, YLIST, YLIST_ERR, norm_indx, f_indx);
  getnorm(X_NORM, X_ERR_NORM, XLIST, YLIST_ERR, norm_indx, f_indx);
  // Type I syst. error. Weighted average +/- weighted standard deviation

  double waverage = get_waverage(YLIST, YLIST_ERR, f_indx);
  double wy_dv = get_wy_sdv(YLIST, YLIST_ERR, f_indx, waverage);
  cout << waverage << "+/-" << wy_dv << endl;
  
  //int length = sizeof(YLIST) / sizeof(YLIST[0]);
  //double syst_err = syst_typeI(YLIST, length);

  // Type II syst. error. 
  // Fill uncorrelated errors
  fill_uncorr(UNCORR_ERR, YLIST, YLIST_ERR, f_indx);
  
  // Create graphs
  // gf_syst: variation as function of cut variable, uncorr. errors
  gf_syst = get_graph_syst(XLIST, YLIST, XLIST_ERR, UNCORR_ERR, f_indx);
  
  // gf_band: norm
  gf_band = get_graph_band(XLIST, BAND, X_ERR_NORM, BAND_ERR, f_indx);
  gf_norm = get_graph_norm(X_NORM, Y_NORM, X_ERR_NORM, Y_ERR_NORM, 1); //new TGraphErrors(1, );
  
  gf_syst -> RemovePoint(norm_indx);
  TMultiGraph *gf_mg = new TMultiGraph();
  gf_mg -> SetTitle("Exclusion graphs");
  gf_mg -> SetName("gf_mg");  
  gf_mg -> Add(gf_band);

  TCanvas *cv_syst = new TCanvas("cv_syst", "syst. " + cut_label, 0, 0, 1000, 800);
  cv_syst -> SetLeftMargin(0.16);
  cv_syst -> SetBottomMargin(0.15);

  gf_syst -> Draw("AP");
  gf_norm -> Draw("P");
  gf_mg -> Draw("C3");

  // save
  //TFile *f_output = new TFile(outputPlot + "plot.root", "recreate");

  cv_syst -> SaveAs(outputPlot + para_label + "_" + cut_label + ".pdf");
  
  //f_output -> Close();

  return 0;
  
}
