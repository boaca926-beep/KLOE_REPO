const TString infile_tmp = "/home/bo/Desktop/input_efficy_norm_TDATA/omega_fit/omega_fit.root";
// efficiency correction factor:/home/bo/Desktop/KLOE_REPO/run_plot/plot_efficy.C
double R_ratio = 1.09474;
double R_ratio_err = 0.0217026;

double rela_err(double a, double b, double delta_a, double delta_b, double c) {

  double ratio1 = delta_a / a;
  double ratio2 = delta_b / b;
  double c_err = c * TMath::Sqrt(ratio1 * ratio1 + ratio2 * ratio2);
  
  return c_err;
  
}

int fit_result_efficy_norm() {

  gErrorIgnoreLevel = kError;
  TGaxis::SetMaxDigits(4);
  //gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetFitFormat("6.4g");

  // check relative error
  double L = 10., delta_L = 0.2;
  double W = 5., delta_W = 0.1;
  double A = L * W;
  double A_err = rela_err(L, W, delta_L, delta_W, A);

  cout << "A = " << A << "+/-" << A_err << endl;
  
  cout << infile_tmp << endl;

  TFile *intree = new TFile(infile_tmp);
  
  TIter next_tree(intree -> GetListOfKeys());
  
  TString objnm_tree, classnm_tree;
  
  int i = 0;
  TKey *key;
  
  while ( (key = (TKey *) next_tree() ) ) {// start tree while lop
    
    i ++;
    
    objnm_tree   =  key -> GetName();
    classnm_tree = key -> GetClassName();
    key -> GetSeekKey();
    
    cout << " tree/histo" << i << ": classnm = " << classnm_tree << ", objnm = " << objnm_tree << endl;
    
  }

  // fill lists
  const int list_size = 1000;

  double M3PI_FIT[list_size], M3PI_FIT_ERR[list_size];
  double N3PI_OBS_FIT[list_size], N3PI_OBS_FIT_ERR[list_size];
  double N3PI_FIT[list_size], N3PI_FIT_ERR[list_size];
  double N3PI_DIFF[list_size], N3PI_DIFF_ERR[list_size];
  double ISRLUMI[list_size], CRX3PI_BW[list_size];
  double EFFICY[list_size], EFFICY_ERR[list_size];
  
  double Mass_omega_fit = 0., Mass_omega_err_fit = 0.;
  double Gam_omega_fit = 0., Gam_omega_err_fit = 0.;
  double BB_fit = 0., BB_err_fit = 0.;
  double chi2_sum = 0.;
  int ndf = 0;
  int bin_indx = 0;
  
  int fit_indx = 0;

  double efficy_corr = 0., efficy_corr_err = 0.;
  
  TTree *TCRX3PI = (TTree*)intree -> Get("TCRX3PI");
  
  for (Int_t irow = 0; irow < TCRX3PI -> GetEntries(); irow++) {// start fill
    
    TCRX3PI -> GetEntry(irow); //cout << irow << endl;

    //
    bin_indx = TCRX3PI -> GetLeaf("Br_bin_indx") -> GetValue(0);

    chi2_sum = TCRX3PI -> GetLeaf("Br_chi2_sum_crx3pi") -> GetValue(0);
    ndf = TCRX3PI -> GetLeaf("Br_ndf") -> GetValue(0);

    Mass_omega_fit = TCRX3PI -> GetLeaf("Br_OMEGA_PARA") -> GetValue(0);
    Mass_omega_err_fit = TCRX3PI -> GetLeaf("Br_OMEGA_PARA_ERR") -> GetValue(0);

    Gam_omega_fit = TCRX3PI -> GetLeaf("Br_OMEGA_PARA") -> GetValue(1);
    Gam_omega_err_fit = TCRX3PI -> GetLeaf("Br_OMEGA_PARA_ERR") -> GetValue(1);
    
    BB_fit = TCRX3PI -> GetLeaf("Br_OMEGA_PARA") -> GetValue(2);
    BB_err_fit = TCRX3PI -> GetLeaf("Br_OMEGA_PARA_ERR") -> GetValue(2);
    
    //
    M3PI_FIT[irow] = TCRX3PI -> GetLeaf("Br_M3PI") -> GetValue(bin_indx);
    M3PI_FIT_ERR[irow] = 0.;

    EFFICY[irow] = TCRX3PI -> GetLeaf("Br_EFFICY") -> GetValue(bin_indx);
    EFFICY_ERR[irow] = TCRX3PI -> GetLeaf("Br_EFFICY_ERR") -> GetValue(bin_indx);
    efficy_corr = EFFICY[irow] * R_ratio;
    //R_ratio = 1;
    //R_ratio_err = 0.;
    efficy_corr_err = rela_err(EFFICY[irow], R_ratio, EFFICY_ERR[irow], R_ratio_err, efficy_corr);
    
    N3PI_OBS_FIT[irow] = TCRX3PI -> GetLeaf("Br_N3PI_OBS") -> GetValue(bin_indx);
    N3PI_OBS_FIT_ERR[irow] = TCRX3PI -> GetLeaf("Br_N3PI_OBS_ERR") -> GetValue(bin_indx);

    N3PI_FIT[irow] = TCRX3PI -> GetLeaf("Br_N3PI_SMEAR") -> GetValue(bin_indx);
    N3PI_FIT_ERR[irow] = TCRX3PI -> GetLeaf("Br_N3PI_SMEAR_ERR") -> GetValue(bin_indx);

    N3PI_DIFF[irow] = TCRX3PI -> GetLeaf("Br_N3PI_DIFF") -> GetValue(bin_indx);
    N3PI_DIFF_ERR[irow] = TCRX3PI -> GetLeaf("Br_N3PI_DIFF_ERR") -> GetValue(bin_indx);
    //N3PI_DIFF[irow] = TCRX3PI -> GetLeaf("Br_N3PI_DIFF") -> GetValue(bin_indx) / N3PI_DIFF_ERR[irow];
    
    //
    ISRLUMI[irow] = TCRX3PI -> GetLeaf("Br_ISRLUMI") -> GetValue(bin_indx);
    CRX3PI_BW[irow] = TCRX3PI -> GetLeaf("Br_CRX3PI_BW") -> GetValue(bin_indx);
    
    fit_indx ++;

    cout << irow + 1 << ": n3pi obs. = " << N3PI_OBS_FIT[irow] << "+/-" << N3PI_OBS_FIT_ERR[irow] << ", n3pi diff. = " << N3PI_DIFF[irow] << "+/-" << N3PI_DIFF_ERR[irow] << endl;

    cout << "list indx = " << irow + 1 << "\n"
	 << "\tbin_indx = " << bin_indx << "\n"
	 << "\tm3pi = " << M3PI_FIT[irow] << " +/- " << M3PI_FIT_ERR[irow] << "\n"
	 << "\tn3pi_obs = " << N3PI_OBS_FIT[irow] << " +/- " << N3PI_OBS_FIT_ERR[irow] << "\n"
	 << "\tn3pi_fit = " << N3PI_FIT[irow] << " +/- " << N3PI_FIT_ERR[irow] << "\n"
	 << "\tn3pi_diff = " << N3PI_DIFF[irow] << " +/- " << N3PI_DIFF_ERR[irow] << "\n"
	 << "\tisr_lumi = " << ISRLUMI[irow] << "\n"
	 << "\tcrx3pi_bw = " << CRX3PI_BW[irow] << "\n"
      //<< "\tefficy = " << EFFICY[irow] * R_ratio * 100. << " +/- " << EFFICY_ERR[irow] * R_ratio_err * 100. << " [%]\n"
	 << "\tefficy = " << efficy_corr * 100. << " +/- " << efficy_corr_err * 100. << " [%], " << EFFICY_ERR[irow] * 100. << "\n"
      	 << "\chi2_sum = " << chi2_sum << "\n\n";
    
  }

  cout << Mass_omega_fit << "+/-" << Mass_omega_err_fit << "\n"
       << Gam_omega_fit  << "+/-" << Gam_omega_err_fit << "\n"
       << BB_fit << "+/-" << BB_err_fit << endl;
    
  return 0;
  
  
}
