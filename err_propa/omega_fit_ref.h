TH1D *hufo;
TH1D *hdata;
TH1D *heeg, *eeg_hist;
TH1D *hksl, *ksl_hist;
TH1D *homegapi, *omegapi_hist;
TH1D *hetagam, *etagam_hist;
TH1D *hmcrest, *mcrest_hist;

TH1D *hsig;
TH1D *hsig_true;  
TH1D *hsig_gen;
TH1D *hefficy;

TH2D *hcorrmatrix;
TH2D *hsmearmatr;
TH2D *hsmearmatr_err;
    
TF1 *fun_bw;

const int list_size = 1000;
int binsize = 0;
double hmin = 0.;
double hmax = 0.;

double M3PI[list_size]; 
double EFFICY[list_size], EFFICY_ERR[list_size]; 
double ISRLUMI[list_size], ISRLUMI_ERR[list_size]; 
double CRX3PI_BW[list_size], CRX3PI_BW_ERR[list_size]; 

double N3PI_SMEAR[list_size], N3PI_SMEAR_ERR[list_size];
double N3PI_OBS[list_size], N3PI_OBS_ERR[list_size];
double N3PI_PRE[list_size], N3PI_PRE_ERR[list_size];
double N3PI_DIFF[list_size], N3PI_DIFF_ERR[list_size];
double BIN_INDX[list_size];
double M3PI_FIT[list_size], M3PI_FIT_ERR[list_size]; 

double chi2_sum_crx3pi = 0.;
int residul_size_crx3pi = 0;
bool fit_condi = false;

TFile *f_hist = new TFile(outputHist + "hist.root"); 
TFile *f_cut = new TFile(outputCut + "tree_pre.root");  
TFile *f_out = new TFile(outputOmega + "omega_fit.root", "recreate");

TTree* TRESULT = new TTree("TRESULT", "recreate");
TTree* TCRX3PI = new TTree("TCRX3PI", "recreate");
  
TList *HIM3pi_fit = (TList *) f_hist -> Get("HIM3pi_fit");
TList *HSIG = (TList *) f_hist -> Get("HSIG");

hsig_true = (TH1D *) HSIG -> FindObject("hsig_true");  
hsig_gen = (TH1D *) HSIG -> FindObject("hsig_gen");

  
  
// parameter gsf, TDATA and TUFO

//double gsf = 1;

//const double crx_scale = 0.3894e6; // 1 [GeV]^{-2}=0.3894 mb=0.3894e6nb

const double fit_min = 760;
const double fit_max = 800;

// methods
double Get_W0_full(double * m3pi) {

  double W0_full = 0.;
  double mm = 0., x = 0.;
  
  mm = m3pi[0] * 1e-3;
  x = 1.-TMath::Power(mm/sqrtS, 2.);
  
  W0_full = (alphapi/x)*(TMath::Log(TMath::Power(sqrtS/me,2))-1)*(2-2*x+x*x);

  //cout << mm << endl;
  
  return W0_full;
 
}

//
double GetISRLumi_apprx(double m3pi, double m3pi_lower, double m3pi_upper, double W0_full) {

  double isrlumi = 0.;
  double Delta_m3pi = (m3pi_upper - m3pi_lower) * 1e-3;
  
  m3pi = m3pi * 1e-3;

  		 
  isrlumi = W0_full * (2. * m3pi / TMath::Power(sqrtS, 2.)) * Lumi_int * Delta_m3pi;
  //isrlumi = W0_full * (2. * m3pi / TMath::Power(sqrtS, 2.)) * 38794.1 * Delta_m3pi;
    
  //cout << "m3pi = " << m3pi << " [GeV/c^2], width = " << Delta_m3pi << " [GeV/c^2], W0_full = " << W0_full << ", isrlumi = " << isrlumi << ", Lumi_int = " << Lumi_int << "\n";

  return isrlumi;
  
}

//
double binomial_err(double nb_true, double nb_gen) {
  double error = 0.;
  double ratio = 0.; 

  if (nb_gen != 0.) {
    ratio = nb_true / nb_gen;
    error = TMath::Sqrt(ratio * (1. - ratio) / nb_gen);
  }
   
  //cout << "true = " << nb_true << ", gen = " << nb_gen << ", ratio = " << ratio << ", error = " << error << endl;

  return error;
}

//
double GetCrx3piMax_new(double BB, double Mass_V) {

  double crx3pi = 12 * pi * BB / Mass_V / Mass_V * 1e6 * crx_scale;

  /*
  cout << "\n\tBB = " << BB << "\n"
       << "\tM_V = " << mass_omega << " [MeV]\n"
       << "\tcrx3pi = " << crx3pi << " [nb]\n";
  */
  
  return crx3pi;
  
}

void getSmearSig() {

  //getObj(f_cut);

  double m3pi_true = 0., m3pi_corr = 0.;

  hsig = new TH1D("hsig", "", binsize, hmin, hmax);
  hsig -> Sumw2();

  TTree *TISR3PI_SIG = (TTree *) f_cut -> Get("TISR3PI_SIG");

  for (Int_t irow = 0; irow < TISR3PI_SIG -> GetEntries(); irow++) {// loop chain

    TISR3PI_SIG -> GetEntry(irow);

    m3pi_true = TISR3PI_SIG -> GetLeaf("Br_IM3pi_true") -> GetValue(0);

    m3pi_corr = DetectorEvent(TMath::Abs(m3pi_true));

    //cout << m3pi_corr << endl;
    
    hsig -> Fill(m3pi_corr);
    
  }
  
}

//
void scaleGSF() {

  cout << gsf << endl;
  heeg -> Scale(gsf); // scale EEG
  homegapi -> Scale(gsf); // scale OMEGAPI
  hetagam -> Scale(gsf); // scale ETAGAM
  hksl -> Scale(gsf); // scale KSL 
  hmcrest -> Scale(gsf); // scale MCREST
  hsig -> Scale(gsf); // scale signal
  
  Lumi_int = Lumi_int * gsf; // scale the total luminosity
  
}




//
void MCNorm() {

  // sig
  hsig -> Scale(sig_sfw);
  
  // eeg
  heeg = (TH1D *) HIM3pi_fit -> FindObject("h1d_IM3pi_TEEG");
  heeg -> Scale(eeg_sfw);
  heeg -> SetName("heeg");
  //heeg -> Draw();
  cout << "nb_eeg = " << heeg -> Integral(1, heeg -> GetNbinsX()) << endl;
  
  // omegapi
  homegapi = (TH1D *) HIM3pi_fit -> FindObject("h1d_IM3pi_TOMEGAPI");
  homegapi -> Scale(omegapi_sfw);
  homegapi -> SetName("homegapi");
  //homegapi -> Draw();
  cout << "nb_omegapi = " << homegapi -> Integral(1, homegapi -> GetNbinsX()) << endl;
  
  // ksl
  hksl = (TH1D *) HIM3pi_fit -> FindObject("h1d_IM3pi_TKSL");
  hksl -> Scale(ksl_sfw);
  hksl -> SetName("hksl");
  //hksl -> Draw();
  cout << "nb_ksl = " << hksl -> Integral(1, hksl -> GetNbinsX()) << endl;
  
  // etagam
  hetagam = (TH1D *) HIM3pi_fit -> FindObject("h1d_IM3pi_TETAGAM");
  hetagam -> Scale(etagam_sfw);
  hetagam -> SetName("hetagam");
  //hetagam -> Draw();
  cout << "nb_etagam = " << hetagam -> Integral(1, hetagam -> GetNbinsX()) << endl;
  
  // mcrest
  TH1D *h1d_IM3pi_TKPM = (TH1D *) HIM3pi_fit -> FindObject("h1d_IM3pi_TKPM");
  TH1D *h1d_IM3pi_TRHOPI = (TH1D *) HIM3pi_fit -> FindObject("h1d_IM3pi_TRHOPI");
  TH1D *h1d_IM3pi_TBKGREST = (TH1D *) HIM3pi_fit -> FindObject("h1d_IM3pi_TBKGREST");
  
  hmcrest = (TH1D*) h1d_IM3pi_TBKGREST -> Clone();
  hmcrest -> Add(h1d_IM3pi_TKPM, 1.);
  hmcrest -> Add(h1d_IM3pi_TRHOPI, 1.);
  hmcrest -> Scale(mcrest_sfw);
  hmcrest -> SetName("hmcrest");
  //hmcrest -> Draw();
  cout << "nb_mcrest = " << hmcrest -> Integral(1, hmcrest -> GetNbinsX()) << endl;
  
}

//
void get_efficy() {

  int binsize = hsig_true -> GetNbinsX();
  double hmin = hsig_true -> GetXaxis() -> GetXmin();
  double hmax = hsig_true -> GetXaxis() -> GetXmax();

  //cout << "binsize = " << binsize << ", hmin = " << hmin << ", hmax = " << hmax << endl;
  
  hefficy = new TH1D("hefficy", "", binsize, hmin, hmax);
  hefficy -> Sumw2();

  double nb_gen = 0., nb_true = 0.;
  double efficy = 0., efficy_err = 0.;
  
  for (int i = 1; i <= binsize; i ++ ) {

    nb_gen = hsig_gen -> GetBinContent(i);
    
    nb_true = hsig_true -> GetBinContent(i);
    
    efficy = nb_true / nb_gen;

    efficy_err = binomial_err(nb_true, nb_gen);
    
    if (nb_true == 0. || nb_gen == 0.) {

      efficy = 0.;
      
      efficy_err = 0.;
      
    }

    hefficy -> SetBinContent(i, efficy);
    hefficy -> SetBinError(i, efficy_err);

  }
  
}

//
double Func_bw(double *x, double *par) {

  double bw = 0.;
  double m_3pi = x[0];
  double m_V = par[0]; // Vector meson mass in GeV/c^2
  double W_V = par[1]; // Vector meson mass width in GeV/c^2
  double BB = par[2]; // B_ee * B_3pi
  double sigma_max = GetCrx3piMax_new(BB, m_V);// 12. * pi * BB / (m_V * m_V); // Normalization factor, 12pi * B_ee * B_3pi / m^2_V

  // s-M^{2}_V + i * sqrt{s} * W_V
  double denomin_sqr_real = TMath::Power(m_3pi * m_3pi - m_V * m_V, 2.) + TMath::Power(m_V * W_V, 2.);
  TComplex denomin(m_3pi * m_3pi - m_V * m_V, m_V * W_V);
  TComplex denomin_sqr = denomin.Rho2();
  
  
  bw = sigma_max * (m_3pi * m_3pi) * (W_V * W_V) / denomin_sqr;
  //bw = 1. / denomin_sqr;

  /*
  cout << "\nm_3pi = " << m_3pi << " [GeV/c^2]\n"
    //<< "sqrtS = " << sqrtS << " [GeV/c^2]\n"
       << "para: m_V = " << m_V << " [GeV/c^2], W_V = " << W_V << " [GeV/c^2], BB = " << BB << "\ n"
       << "denomin = (" << denomin.Re() << ", " << denomin.Im() << "), Rho2 (check) = " << denomin.Rho2() << "(" << denomin_sqr_real << ") [GeV/c^2]\n"
       << "sigma_max = " << sigma_max << "\n" 
       << "bw = " << bw << "\n";
  */

  //cout << "sigma_max = " << sigma_max << ", m_V = " << m_V << ", BB = " << BB << endl;

  //cout << "m3pi = " << m_3pi << ", bw = " << bw << endl;

  return bw;//nb
  
}

//
TH2D *getCorrMatrix(int bins, double xmin, double xmax) {

  TTree *TISR3PI_SIG = (TTree*)f_cut -> Get("TISR3PI_SIG");

  // create corrlation matrix
  TH2D *h2d_corr = new TH2D("", "", bins, xmin, xmax, bins, xmin, xmax);
  
  //cout << "correlation matrix is created!!! " << TISR3PI_SIG -> GetName() << endl;

  double IM3pi_true = 0.;
  double IM3pi_corred = 0.;
  double PARA[5] = {frac, smallBias, smallSigma, wideBias, wideSigma};
  
  for (Int_t irow = 0; irow < TISR3PI_SIG -> GetEntries(); irow++) {// loop chain

    TISR3PI_SIG -> GetEntry(irow);
    
    IM3pi_true = TISR3PI_SIG -> GetLeaf("Br_IM3pi_true") -> GetValue(0);
    
    //IM3pi_corred = DetectorEvent(TMath::Abs(IM3pi_true));
    IM3pi_corred = DetectorEvent_fcn(TMath::Abs(IM3pi_true), PARA);
    //cout << IM3pi_corred << ", " << DetectorEvent(TMath::Abs(IM3pi_true), para) << endl;
    
    h2d_corr -> Fill(IM3pi_true, IM3pi_corred);
    
  }
  
  return h2d_corr;
  
}

//
double get_ratio(double nb_true, double nb_gen) {

  double ratio = 0.;
  
  if (nb_gen != 0.) {

    ratio = nb_true / nb_gen;
    
  }
  
  return ratio;
}

//
TH2D* getSmearMatrix(TH2D *h2d_corr, int bins, double xmin, double xmax) {

  TH2D *h2d_smear = new TH2D("", "", bins, xmin, xmax, bins, xmin, xmax);
  
  double z_value = 0.;
  double nb_true = 0.;
  double p_value = 0., p_value_err = 0.;

  int entries = 0;
  
  for (int j = 1; j <= h2d_corr -> ProjectionX() -> GetNbinsX(); j ++ ) {

    for (int i = 1; i <= h2d_corr -> ProjectionY() -> GetNbinsX(); i ++ ) {

      z_value = h2d_corr -> GetBinContent(j, i);
      nb_true = h2d_corr -> ProjectionX() -> GetBinContent(j);
      
      p_value = get_ratio(z_value, nb_true);

      if (p_value != 0.) {

	p_value_err = binomial_err(z_value, nb_true) / p_value; //Relative err

      }
      else {

	p_value_err = 0.;

      }

      h2d_smear -> SetBinContent(j, i, p_value);
      //h2d_smear_matrix_err -> SetBinContent(j, i, p_value_err);
      entries ++;

    }
	
  }
  
  return h2d_smear;
  
}

//
TH2D* getSmearMatrixErr(TH2D *h2d_corr, int bins, double xmin, double xmax) {

  TH2D *h2d = new TH2D("", "", bins, xmin, xmax, bins, xmin, xmax);
  
  double z_value = 0.;
  double nb_true = 0.;
  double p_value = 0., p_value_err = 0.;

  int entries = 0;
  
  for (int j = 1; j <= h2d_corr -> ProjectionX() -> GetNbinsX(); j ++ ) {

    for (int i = 1; i <= h2d_corr -> ProjectionY() -> GetNbinsX(); i ++ ) {

      z_value = h2d_corr -> GetBinContent(j, i);
      nb_true = h2d_corr -> ProjectionX() -> GetBinContent(j);
      
      p_value = get_ratio(z_value, nb_true);

      if (p_value != 0.) {

	//p_value_err = binomial_err(z_value, nb_true) / p_value; //Relative err
	p_value_err = binomial_err(z_value, nb_true); //Relative err

      }
      else {

	p_value_err = 0.;

      }

      //cout << z_value << ", " << nb_true << ", " << p_value << "+/-" << p_value_err << endl;
      
      h2d -> SetBinContent(j, i, p_value_err);
      entries ++;

    }
	
  }
  
  return h2d;
  
}


void fcn_crx3pi(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag){

  double Mass_V = par[0];
  double Gam_V = par[1];
  double BB = par[2];

  fun_bw -> SetParameters(Mass_V, Gam_V, BB);

  int counter = 0;
  double n3pi_pre_sigma2 = 0.;
  double term3 = 0.;

  cout << Mass_V << endl;
  
  for (Int_t irow = 0; irow < TRESULT -> GetEntries(); irow++) {// begin for loop

    TRESULT -> GetEntry(irow);

    M3PI[counter] = TRESULT -> GetLeaf("Br_m3pi") -> GetValue(0);

    EFFICY[counter] = TRESULT -> GetLeaf("Br_efficy") -> GetValue(0);
    EFFICY_ERR[counter] = TRESULT -> GetLeaf("Br_efficy_err") -> GetValue(0);

    ISRLUMI[counter] = TRESULT -> GetLeaf("Br_isrlumi_apprx") -> GetValue(0);
    ISRLUMI_ERR[counter] = 0.;
    
    CRX3PI_BW[counter] = fun_bw -> Eval(M3PI[counter]);
    CRX3PI_BW_ERR[counter] = 0.;
    
    // updated from TUFO
    N3PI_OBS[counter] = TRESULT -> GetLeaf("Br_nb_isr3pi_obs") -> GetValue(0);
    N3PI_OBS_ERR[counter] = TRESULT -> GetLeaf("Br_nb_isr3pi_obs_err") -> GetValue(0);
	       
    // predicted number of 3pi events after the acceptance correction
    // apply globle scaling factor
    N3PI_PRE[counter] = ISRLUMI[counter] * CRX3PI_BW[counter] * EFFICY[counter];
    //n3pi_pre_sigma2 = get_sigma2_y(ISRLUMI[counter], CRX3PI_BW[counter], EFFICY[counter], ISRLUMI_ERR[counter], CRX3PI_BW_ERR[counter], EFFICY_ERR[counter]);
    //N3PI_PRE_ERR[counter] = TMath::Sqrt(n3pi_pre_sigma2);
    N3PI_PRE_ERR[counter] = ISRLUMI[counter] * CRX3PI_BW[counter] * EFFICY_ERR[counter];
  

    //cout << "ISRLUMI = " << ISRLUMI[counter] << ", CRX3PI_BW = " << CRX3PI_BW[counter] << "\n";
    //cout << "N3PI_OBS = " << N3PI_OBS[counter] << "\n";
    //cout << "N3PI_PRE = " << N3PI_PRE[counter] << " +/- " << N3PI_PRE_ERR[counter] << ", " << term3 << "\n";
    

    /*
    cout << "bin indx = " << irow + 1 << ", M3PI = " << M3PI[counter] << ", EFFICY [%] = " << EFFICY[counter] * 100. << " +/- " << EFFICY_ERR[counter] * 100. << ", ISRLUMI = " << ISRLUMI[counter] << ", CRX3PI_BW = " << CRX3PI_BW[counter] << ", N3PI_OBS = " << N3PI_OBS[counter] << " +/- " << N3PI_OBS_ERR[counter] << ", N3PI_PRE = " << N3PI_PRE[counter] << "\n";
    */
    
    counter ++;
    
  }

  //cout << "\n"; 

  
  
  // smearing
  
  int total_bins = 0;

  for (int j = 1; j <= hcorrmatrix -> ProjectionY() -> GetNbinsX(); j ++ ) {

    double prob_tmp = 0., prob_tmp_err = 0.;
    
    for (int i = 1; i <= hcorrmatrix -> ProjectionX() -> GetNbinsX(); i ++ ) {
      
      prob_tmp = hsmearmatr -> GetBinContent(i, j);
      prob_tmp_err = hsmearmatr -> GetBinError(i, j);

      cout << "smearing matrix entries (i, j) = " << prob_tmp << " +/- " << prob_tmp_err << endl;
      
    }
  }

  //fit_condi = true;
  if (fit_condi) return;
 
  
  //
  double diff_sqr = 0., err_sqr = 0.;
  double chi2_sum_tmp = 0., chi2_tmp = 0., chi2_tmp_err = 0.;
  
  int fit_indx = 0;
  int range_indx = 0, sub_indx = 0;

  for (int i = 0; i < counter; i ++) {
    
    if (M3PI[i] > fit_min && M3PI[i] < fit_max) {

      diff_sqr = TMath::Abs((N3PI_SMEAR[i] - N3PI_OBS[i]) * (N3PI_SMEAR[i] - N3PI_OBS[i]));
      err_sqr = N3PI_OBS_ERR[i] * N3PI_OBS_ERR[i] + N3PI_SMEAR_ERR[i] * N3PI_SMEAR_ERR[i];

      chi2_tmp = diff_sqr / err_sqr;
      chi2_sum_tmp += chi2_tmp;
      
      BIN_INDX[fit_indx] = i;
      N3PI_DIFF[i] = N3PI_OBS[i] - N3PI_SMEAR[i];
      N3PI_DIFF_ERR[i] = TMath::Sqrt(N3PI_OBS_ERR[i] * N3PI_OBS_ERR[i] + N3PI_SMEAR_ERR[i] * N3PI_SMEAR_ERR[i]);

      fit_indx ++;
      
    }
      
  }

  //chi2_sum_crx3pi = chi2_sum_tmp,
  //residul_size_crx3pi = fit_indx;

  f = 0.; //chi2_sum_tmp;
  
}

//
double GetBB_new(double sigma_max, double Mass_V) {

  double BB = sigma_max * 1e-6 / crx_scale * Mass_V * Mass_V / 12. / pi;

  return BB;

}


//
double GetBkgErr(double NData, double nMC, double NMC, double f, double f_err) {

  double r = nMC / NMC;
  double r_err = 0.; 
  double r_ratio = 0.;
  
  if (nMC != 0.) {
    r = nMC / NMC;
    r_err = binomial_err(nMC, NMC);
    r_ratio = r_err / r;
  }
  
  double f_ratio = f_err / f;
  
  double n_scaled = NData * r * f;

  //cout << "n_scaled = " << n_scaled << endl;
  
  double n_scaled_err = n_scaled * TMath::Sqrt(r_ratio * r_ratio + f_ratio * f_ratio);
  
  return n_scaled_err;
  
}

