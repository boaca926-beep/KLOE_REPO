const TString inpath = "../../analysis/chains_qa";
const double sample_size = 8373;
const double Lumi_int = 1724470 * (sample_size / 8373);
const double sample_frac = Lumi_int / 1724470 * 100.;
const double Delta_FWHM = 1.76707;
const double Delta_BB = 0.0000101533;
//Cut conditions:
const bool cut_cond = true;
const TString cut_nm = "";
//Cut paramters:
const double chi2_cut = 43;
const double angle_cut = 138;
const double deltaE_cut = -150;
const double beta_cut = 1.98;
const int cut_indx = 0;
//omega peak correction paramters:
const double shift_fact = 0.3;
const double smear_fact = 3.8;
const double isr3pi_sfw1d = 0.0422997;
const double isr3pi_sfw1d_err = 0.00532178;
//Constants:
const double c0 = 0.11;
const double c1 = 0.8;
//Hist attributes:
const double mass_sigma_nb = 1;
const double sfw2d_sigma_nb = 1;
//Eisr:
const double Eisr_min = 50;
const double Eisr_max = 500;
const double Eisr_sigma = 2.48;
const int Eisr_bin = TMath::Nint((Eisr_max - Eisr_min) / sfw2d_sigma_nb / Eisr_sigma);
//ppIM:
const double ppIM_min = 200;
const double ppIM_max = 700;
const double ppIM_sigma = 2.30;
const int ppIM_bin = TMath::Nint((ppIM_max - ppIM_min) / sfw2d_sigma_nb / ppIM_sigma);
//IM3pi:
const double IM3pi_min = 380;
const double IM3pi_max = 1020;
const double IM3pi_sigma = 2.65;
const int IM3pi_bin = TMath::Nint((IM3pi_max - IM3pi_min) / mass_sigma_nb / IM3pi_sigma);
const double fit_min = 760;
const double fit_max = 800;
//Double guassian smearing
const double wideSigma_sigma_nb = 1;
const double smallSigma_sigma_nb = 1;
const double frac = 0.20;
const double wideBias = 0.5;
const double wideSigma = 13;
const double smallBias = 0.34;
const double smallSigma = 3.2;
//mass bins
const double mass_size = 5.30;
const int mass_steps = TMath::Nint((fit_max - fit_min) / mass_size);
//smearing matrix
const double rec_min = 600;
const double rec_max = 900;
//Cut information (nominal)
double cut_value = 0;

TTree* TSFW2D = new TTree("TSFW2D", "recreate");
TTree* TRESULT = new TTree("TRESULT", "recreate");
  
double nb_data_sum = 0., nb_eeg_sum = 0., nb_ksl_sum = 0., nb_omegapi_sum = 0., nb_etagam_sum = 0., nb_isr3pi_sum = 0., nb_mcrest_sum = 0.;

double chi2_sfw2d_sum = 0., residul_size_sfw2d = 0.;
double chi2_sum_crx3pi = 0.;
int residul_size_crx3pi = 0;

const int list_size = 1000;

double M3PI[list_size], M3PI_LOWER[list_size], M3PI_UPPER[list_size];
double EFFICY[list_size], EFFICY_ERR[list_size];
double ISRLUMI[list_size];
double CRX3PI_BW[list_size];

double N3PI_TRUE[list_size];
double N3PI_GEN[list_size];
double N3PI_CORR[list_size];
double N3PI_PRE[list_size];
double N3PI_SMEAR[list_size], N3PI_SMEAR_ERR[list_size];
double N3PI_OBS[list_size];
double N3PI_OBS_ERR[list_size];
double N3PI_DIFF[list_size], N3PI_DIFF_ERR[list_size];
double BIN_INDX[list_size];

TF1 *fun_bw;

TRandom *rnd=0;

//TH1D *h1d_IM3pi_ALLCHAIN_GEN;

TH2D *hsmearmatr = new TH2D("hsmearmatr", "", IM3pi_bin, IM3pi_min, IM3pi_max, IM3pi_bin, IM3pi_min, IM3pi_max);

TH2D *hsmearmatr_err = new TH2D("hsmearmatr_err", "", IM3pi_bin, IM3pi_min, IM3pi_max, IM3pi_bin, IM3pi_min, IM3pi_max);
 
  
TH2D *h2d_sfw_TDATA;
TH2D *h2d_sfw_TEEG; //done
TH2D *h2d_sfw_TOMEGAPI;
TH2D *h2d_sfw_TETAGAM;
TH2D *h2d_sfw_TKSL;
TH2D *h2d_sfw_TKPM;
TH2D *h2d_sfw_TISR3PI_SIG;
TH2D *h2d_sfw_TBKGREST;
TH2D *h2d_sfw_TRHOPI;
TH2D *h2d_sfw_MCREST;
TH2D *h2d_sfw_MCSUM;

TH1D *h1d_IM3pi_TDATA;
TH1D *h1d_IM3pi_TEEG; //done
TH1D *h1d_IM3pi_TOMEGAPI;
TH1D *h1d_IM3pi_TETAGAM;
TH1D *h1d_IM3pi_TKSL;
TH1D *h1d_IM3pi_TKPM;
TH1D *h1d_IM3pi_TISR3PI_SIG;
TH1D *h1d_IM3pi_TBKGREST;
TH1D *h1d_IM3pi_TRHOPI;
TH1D *h1d_IM3pi_TISR3PI_SIG_TRUE;
TH1D *h1d_IM3pi_TISR3PI_SIG_GEN;

double GetFBeta(double a1_temp, double b1_temp, double c1_temp, double m2pi_temp) {
  m2pi_temp = m2pi_temp / 1000.;
  double fbeta = a1_temp + 1. / (exp((m2pi_temp - c1_temp) / b1_temp) - 1.);
  /*cout << "a1 = " << a1 << ", a2 = " << a2 << "\n"
    << "b1 = " << b1 << ", b2 = " << b2 << "\n"
    << "c1 = " << c1 << ", c2 = " << c2 << "\n\n";*/
  //cout << "fbeta = " << fbeta << endl;
  return fbeta;
}

double W0_full_fun(double *xx, double *par) {

  double W0_full = 0.;
  double x = 0.;
  
  x = xx[0];

  W0_full = (alphapi/x)*(TMath::Log(TMath::Power(sqrtS/me,2))-1)*(2-2*x+x*x);

  //cout << "(par1, par2, par3) = (" << par[0] << ", " << par[1] << ", " << par[2] << ")" << endl;
  //cout << "x = " << x << ", W0_full = " << W0_full << "\n";
  
  return W0_full;
}

double GetISRLumi_exact(double m3pi_lower, double m3pi_upper) {

  double x_upper = 1. - TMath::Power(m3pi_lower * 1e-3 / sqrtS, 2.);
  double x_lower = 1. - TMath::Power(m3pi_upper * 1e-3 / sqrtS, 2.);

  TF1 *func_intW0 = new TF1("fit_intW0", W0_full_fun, x_lower, x_upper, 3);

  // try to set parameters, if there is any
  func_intW0 -> SetParameters(4, 5, 6); // some dummy variables
  
  double integral = func_intW0 -> Integral(x_lower, x_upper);

  double isrlumi = integral * Lumi_int;

  /*
  cout << "theta0 = " << theta0 << ", sqrtS = " << sqrtS << " [GeV/c^2] \n"
       << "m3pi range (" << m3pi_lower << ", " << m3pi_upper << ") [MeV/c^2] \n"
       << "x range (" << x_lower << ", " << x_upper << ") \n"
    //<< "evaluate func_inW0 at m3pi = 1018.82 GeV, x = " << 1. - TMath::Power(1018.82 * 1e-3 / sqrtS, 2.) << ", func_intW0 -> Eval = " << func_intW0 -> Eval(1. - TMath::Power(1018.82 * 1e-3 / sqrtS, 2.)) << "\n" // W0 = 186.672
       << "integral = " << integral << ", Lumi_int = " << Lumi_int << " [nb^-1]\n"
       << "isrlumi = Lumi_int * integral = " << isrlumi << " [nb^-1]\n\n";
  */
  
  return isrlumi;
  
}

double GetISRLumi_apprx(double m3pi, double m3pi_lower, double m3pi_upper, double W0_full) {

  double isrlumi = 0.;
  double Delta_m3pi = (m3pi_upper - m3pi_lower) * 1e-3;
  
  m3pi = m3pi * 1e-3;

  		 
  isrlumi = W0_full * (2. * m3pi / TMath::Power(sqrtS, 2.)) * Lumi_int * Delta_m3pi;
  //isrlumi = W0_full * (2. * m3pi / TMath::Power(sqrtS, 2.)) * 38794.1 * Delta_m3pi;
    
  //cout << "m3pi = " << m3pi << " [GeV/c^2], width = " << Delta_m3pi << " [GeV/c^2], W0_full = " << W0_full << ", isrlumi = " << isrlumi << ", Lumi_int = " << Lumi_int << "\n";

  return isrlumi;
  
}

double Get_W0_full(double * m3pi) {

  double W0_full = 0.;
  double mm = 0., x = 0.;
  
  mm = m3pi[0] * 1e-3;
  x = 1.-TMath::Power(mm/sqrtS, 2.);
  
  W0_full = (alphapi/x)*(TMath::Log(TMath::Power(sqrtS/me,2))-1)*(2-2*x+x*x);

  //cout << mm << endl;
  
  return W0_full;
 
}


double get_ratio(double nb_true, double nb_gen) {

  double ratio = 0.;
  
  if (nb_gen != 0.) {

    ratio = nb_true / nb_gen;
    
  }
  
  return ratio;
}

double GetProbErr_sqr(double prob, double rel_err) {

  double error = prob * rel_err;

  //cout << "prob = " << prob << ", rel_err = " << rel_err << ", prob_err = " << error << endl;
  
  return error * error;

}

double getscale(double Nd, double fra, double N){
  double value = 0.;

  value =  Nd * fra / N;

  return value;
}

double getloglh(double n_d, double mu) {
  double value = 0.;

  //value = n_d * TMath::Log(mu) - mu + n_d - n_d * TMath::Log(n_d);
  //value = n_d * TMath::Log(mu) - mu - (n_d * TMath::Log(n_d) - n_d + 0.5 * TMath::Log(n_d) + 0.5 * TMath::Log(2 * TMath::Pi()));
  value = n_d * TMath::Log(mu) - mu;
  //cout << "n_d = " << n_d << endl;

  return value;
}

void fcn_sfw2d(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag){


  int counter = 0;

  double nb_data = 0.;
  double nb_eeg = 0.;
  double nb_ksl = 0.;
  double nb_omegapi = 0.;
  double nb_etagam = 0.;
  double nb_isr3pi = 0.;
  double nb_mcrest = 0.;

  double f_eeg = par[0];
  double f_isr3pi = par[1];
  double f_omegapi = par[2];
  double f_etagam = par[3];
  double f_ksl = par[4];
  double f_mcrest = par[5];

  double eeg_sfw = 0., eeg_mu = 0.;
  double isr3pi_sfw = 0., isr3pi_mu = 0.;
  double omegapi_sfw = 0., omegapi_mu = 0.;
  double etagam_sfw = 0., etagam_mu = 0.;
  double ksl_sfw = 0., ksl_mu = 0.;
  double mcrest_sfw = 0., mcrest_mu = 0.;

  double mu_tmp = 0.;

  double chi2_sum_tmp = 0.;
  double llh_sum = 0.;

  
  //cout << "nb_data_sum = " << nb_data_sum << endl;
  
  for (Int_t irow = 0; irow < TSFW2D -> GetEntries(); irow++) {// fcn loop
    
    TSFW2D -> GetEntry(irow);

    // data
    nb_data = TSFW2D -> GetLeaf("Br_nb_data") -> GetValue(0);

    // eeg
    nb_eeg = TSFW2D -> GetLeaf("Br_nb_eeg") -> GetValue(0);
    eeg_sfw = getscale(nb_data_sum, f_eeg, nb_eeg_sum);
    eeg_mu = nb_eeg * eeg_sfw;

    // ksl
    nb_ksl = TSFW2D -> GetLeaf("Br_nb_ksl") -> GetValue(0);
    ksl_sfw = getscale(nb_data_sum, f_ksl, nb_ksl_sum);
    ksl_mu = nb_ksl * ksl_sfw;

    // omegapi
    nb_omegapi = TSFW2D -> GetLeaf("Br_nb_omegapi") -> GetValue(0);
    omegapi_sfw = getscale(nb_data_sum, f_omegapi, nb_omegapi_sum);
    omegapi_mu = nb_omegapi * omegapi_sfw;

    // isr3pi
    nb_isr3pi = TSFW2D -> GetLeaf("Br_nb_isr3pi") -> GetValue(0);
    isr3pi_sfw = getscale(nb_data_sum, f_isr3pi, nb_isr3pi_sum);
    isr3pi_mu = nb_isr3pi * isr3pi_sfw;
    
    // etagam
    nb_etagam = TSFW2D -> GetLeaf("Br_nb_etagam") -> GetValue(0);
    etagam_sfw = getscale(nb_data_sum, f_etagam, nb_etagam_sum);
    etagam_mu = nb_etagam * etagam_sfw;

    // mcrest
    nb_mcrest = TSFW2D -> GetLeaf("Br_nb_mcrest") -> GetValue(0);
    mcrest_sfw = getscale(nb_data_sum, f_mcrest, nb_mcrest_sum);
    mcrest_mu = nb_mcrest * mcrest_sfw;

    // mu
    mu_tmp = eeg_mu + ksl_mu + etagam_mu + omegapi_mu + mcrest_mu + isr3pi_mu;

    if (mu_tmp > 0. && nb_data > 0.) {
    
      //
      /*
      if (irow == 55) { 
	cout << "indx = " << irow << "\n"
	     << "nb_data(sum) = " << nb_data << "(" << nb_data_sum << ")\n"
	     << "nb_eeg(sum) = " << nb_eeg << "(" << nb_eeg_sum << "), f = " << f_eeg << ", sfw = " << eeg_sfw << ", mu = " << eeg_mu << "\n"
	     << "nb_ksl(sum) = " << nb_ksl << "(" << nb_ksl_sum << "), f = " << f_ksl << ", sfw = " << ksl_sfw << ", mu = " << ksl_mu << "\n"
	     << "nb_omegapi(sum) = " << nb_omegapi << "(" << nb_omegapi_sum << "), f = " << f_omegapi << ", sfw = " << omegapi_sfw << ", mu = " << omegapi_mu << "\n"
	     << "nb_isr3pi(sum) = " << nb_isr3pi << "(" << nb_isr3pi_sum << "), f_isr3pi = " << f_isr3pi << ", sfw = " << isr3pi_sfw << ", mu = " << isr3pi_mu << "\n"
	     << "nb_etagam(sum) = " << nb_etagam << "(" << nb_etagam_sum << "), f_etagam = " << f_etagam << ", sfw = " << etagam_sfw << ", mu = " << etagam_mu << "\n"
	     << "nb_mcrest(sum) = " << nb_mcrest << "(" << nb_mcrest_sum << "), f_mcrest = " << f_mcrest << ", sfw = " << mcrest_sfw << ", mu = " << mcrest_mu << "\n\n";
	
      }
      */
      
      chi2_sum_tmp += (nb_data - mu_tmp) * (nb_data - mu_tmp) / (nb_data + mu_tmp);
      //chi2_sum_tmp += (nb_data - mu_tmp) * (nb_data - mu_tmp) / nb_data;
      
      counter ++;
      
      llh_sum -= 2. * getloglh(nb_data, mu_tmp);
      
    }

  }// end fcn loop 

  chi2_sfw2d_sum = chi2_sum_tmp;
  residul_size_sfw2d = counter;
   
  f = llh_sum;
  //f = chi2_sfw2d_sum;
    
}

double GetScalError(double N_d, double N, double f, double f_error) {
  double error = 0.;
  double scale = N_d * f / N;
  error = scale * TMath::Sqrt(1. / N_d + 1. / N + TMath::Power(f_error / f, 2));

  //cout << "D = " << N_d << ", N = "<< N << ", f = " << f << "+/-" << f_error << ", scale = " << scale << " +/- " << error << endl;

  return error;

}



double getloglh_exact(double n_d, double mu) {
  double value = 0.;

  value = n_d * TMath::Log(mu) - mu - TMath::Log(TMath::Factorial(TMath::Nint(n_d)));

  //cout << "value = " << value << ", TMath::Log(TMath::Factorial(TMath::Nint(n_d))) = " <<  TMath::Log(TMath::Factorial(TMath::Nint(n_d))) << ", n_d = "<< n_d << endl;
  
  return value;
}


void fcn_crx3pi(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag){

  //const int list_size = 1000;

  double Mass_V = par[0];
  //double Gam_V = par[1] + Delta_FWHM;
  //double BB = par[2] - Delta_BB;

  double Gam_V = par[1];
  double BB = par[2];

  //double M3PI[list_size];
  //double EFFICY[list_size];
  //double ISRLUMI[list_size];
  //double CRX3PI_BW[list_size];

  //double N3PI_TRUE[list_size];
  //double N3PI_GEN[list_size];
  //double N3PI_CORR[list_size];
  //double N3PI_PRE[list_size];
  //double N3PI_SMEAR[list_size];
  //double N3PI_OBS[list_size];
  //double N3PI_OBS_ERR[list_size];
	       
  //TF1 *fun_bw = new TF1("func_bw", Func_bw, 0., 2000., 3);
  fun_bw -> SetParameters(Mass_V, Gam_V, BB);

  int counter = 0;
  
  for (Int_t irow = 0; irow < TRESULT -> GetEntries(); irow++) {// begin for loop

    TRESULT -> GetEntry(irow);

    M3PI[counter] = TRESULT -> GetLeaf("Br_m3pi") -> GetValue(0);
    EFFICY[counter] = TRESULT -> GetLeaf("Br_efficy") -> GetValue(0);
    EFFICY_ERR[counter] = TRESULT -> GetLeaf("Br_efficy_err") -> GetValue(0);
    ISRLUMI[counter] = TRESULT -> GetLeaf("Br_isrlumi_apprx") -> GetValue(0);
    CRX3PI_BW[counter] = fun_bw -> Eval(M3PI[counter]);

    //N3PI_TRUE[counter] = TRESULT -> GetLeaf("Br_nb_isr3pi_true") -> GetValue(0);
    //N3PI_GEN[counter] = TRESULT -> GetLeaf("Br_nb_isr3pi_gen") -> GetValue(0);
    //N3PI_CORR[counter] = TRESULT -> GetLeaf("Br_nb_isr3pi_corred") -> GetValue(0);
    N3PI_OBS[counter] = TRESULT -> GetLeaf("Br_nb_isr3pi_obs") -> GetValue(0);
    N3PI_OBS_ERR[counter] = TRESULT -> GetLeaf("Br_nb_isr3pi_obs_err") -> GetValue(0);
	       
    // predicted number of 3pi events after the acceptance correction
    N3PI_PRE[counter] = ISRLUMI[counter] * CRX3PI_BW[counter] * EFFICY[counter];

    //cout << "ISRLUMI = " << ISRLUMI[counter] << ", CRX3PI_BW = " << CRX3PI_BW[counter] << "\n";
    //cout << "N3PI_OBS = " << N3PI_OBS[counter] << "\n";

    /*
    cout << "bin indx = " << irow + 1 << ", M3PI = " << M3PI[counter] << ", EFFICY [%] = " << EFFICY[counter] * 100. << " +/- " << EFFICY_ERR[counter] * 100. << ", ISRLUMI = " << ISRLUMI[counter] << ", CRX3PI_BW = " << CRX3PI_BW[counter] << ", N3PI_OBS = " << N3PI_OBS[counter] << " +/- " << N3PI_OBS_ERR[counter] << ", N3PI_PRE = " << N3PI_PRE[counter] << "\n";
    */
    
    counter ++;
    
  }

  //cout << "\n"; 

  
  
  // smearing
  
  for (int j = 1; j <= hsmearmatr -> ProjectionY() -> GetNbinsX(); j ++ ) {

    double nb_rec_tmp = 0.;
    double nb_rec_tmp_err2_sum = 0.; 
    double nb_gen = 0., nb_gen_sum = 0.;

    double prob_sum = 0.;
    double prob_tmp = 0.;
      
    for (int i = 1; i <= hsmearmatr -> ProjectionX() -> GetNbinsX(); i ++ ) {

      prob_tmp = hsmearmatr -> GetBinContent(i, j);
      prob_sum += prob_tmp;

      nb_gen = h1d_IM3pi_TISR3PI_SIG_GEN -> GetBinContent(i);
      nb_gen_sum += nb_gen;
      
      //nb_rec_tmp += prob_tmp * N3PI_TRUE[i - 1];
      nb_rec_tmp += prob_tmp * N3PI_PRE[i - 1];
      nb_rec_tmp_err2_sum += nb_gen * prob_tmp * (1 - prob_tmp) * (N3PI_PRE[i - 1] * N3PI_PRE[i - 1]);
      
      //cout << "(rec_indx, true_indx) = (" << j << ", " << i << "), N3pi_true = " << N3PI_TRUE[i - 1] << ", prob_tmp = " << prob_tmp << ", nb_rec_tmp = " << nb_rec_tmp << ", N3pi_corr = " << N3PI_CORR[i - 1] << "\n";

      //cout << "(rec_indx, true_indx) = (" << j << ", " << i << "), prob_tmp = " << prob_tmp << ", N3PI_PRE = " << N3PI_PRE[i - 1] << ", nb_rec_tmp = " << nb_rec_tmp << "\n";
      
      //cout << "(rec_indx, true_indx) = (" << j << ", " << i << "), prob_tmp = " << prob_tmp << ", N3PI_PRE = " << N3PI_PRE[i - 1] << ", nb_rec_tmp = " << nb_rec_tmp << "\n";
      
      
    }

    N3PI_SMEAR[j -1] = nb_rec_tmp;
    //N3PI_SMEAR_ERR[j -1] = TMath::Sqrt(nb_rec_tmp); //0.1 * N3PI_SMEAR[j - 1];
    N3PI_SMEAR_ERR[j -1] = TMath::Sqrt(nb_rec_tmp_err2_sum) / nb_gen_sum; //0.1 * N3PI_SMEAR[j - 1];

    //cout << "bin indx = " << j << ", N3PI_SMEAR = " << N3PI_SMEAR[j - 1] << " +/- " << N3PI_SMEAR_ERR[j - 1] << "\n";

    //cout << "\n"; 

    
  }

  //cout << "\n"; 

  // mass interval
  
      
  //
  double diff_sqr = 0., err_sqr = 0.;
  double chi2_sum_tmp = 0., chi2_tmp = 0., chi2_tmp_err = 0.;
  
  int fit_indx = 0;
  int range_indx = 0, sub_indx = 0;

  for (int i = 0; i < counter; i ++) {
    
    if (M3PI[i] > fit_min && M3PI[i] < fit_max) {
      
      /*
	if (M3PI[i] > M3PI_LOWER[range_indx] && range_indx < mass_steps) {
	
	range_indx ++;
	
	sub_indx = 0;
	
	}
	
	if (M3PI[i] > M3PI_LOWER[range_indx - 1] && M3PI[i] < M3PI_UPPER[range_indx - 1]) {
	
	sub_indx ++;
	
	}
      */
      
      //cout << sub_indx << ", " << M3PI[i] << ", " << M3PI_LOWER[range_indx - 1] << ", " << M3PI_UPPER[range_indx - 1] << endl;
      
      //cout << h1d_IM3pi_TISR3PI_SIG_TRUE -> GetBinContent(i) << endl;
      
      diff_sqr = TMath::Abs((N3PI_SMEAR[i] - N3PI_OBS[i]) * (N3PI_SMEAR[i] - N3PI_OBS[i]));
      
      err_sqr = N3PI_OBS_ERR[i] * N3PI_OBS_ERR[i] + N3PI_SMEAR_ERR[i] * N3PI_SMEAR_ERR[i];

      chi2_tmp = diff_sqr / err_sqr;
      //chi2_tmp_err = 0.1 * chi2_tmp;
      
      chi2_sum_tmp += chi2_tmp;
      
      BIN_INDX[fit_indx] = i;
      

      N3PI_DIFF_ERR[i] = TMath::Sqrt(N3PI_OBS_ERR[i] * N3PI_OBS_ERR[i] + N3PI_SMEAR_ERR[i] * N3PI_SMEAR_ERR[i]);
      N3PI_DIFF[i] = N3PI_OBS[i] - N3PI_SMEAR[i];

      /*
      cout << "\t/////////////////////////////////\n"
	   << "\tbin " << i << "\n"
	   << "\tMass_V = " << Mass_V << ", Gam_V = " << Gam_V << ", BB = " << BB << "\n" 
	   << "\tm3pi = "   << M3PI[i] << " [MeV/c^2]\n"
	   << "\tefficy [%] = " << EFFICY[i] * 100. << " +/- " << EFFICY_ERR[i] << "\n"
	   << "\tisrlumi = " << ISRLUMI[i] << "[nb^-1]\n"
	   << "\tcrx3pi_bw = " << CRX3PI_BW[i] << "\n"
	//<< "\tnb_isr3pi_gen = " << N3PI_TRUE[i] << ", efficy = " << EFFICY[i] << ", gen * efficy (true) = " << N3PI_GEN[i] * EFFICY[i] << "(" << N3PI_TRUE[i] << ")\n"
	//<< "\tnb_isr3pi_corr = " << N3PI_CORR[i] << "\n"
	   << "\tnb_isr3pi_obs = " << N3PI_OBS[i] << " +/- " << N3PI_OBS_ERR[i] << "\n"
	   << "\tnb_isr3pi_fit = " << N3PI_SMEAR[i] << " +/- " << N3PI_SMEAR_ERR[i] << "\n"
	   << "\tn3pi_diff = obs - fit = " << N3PI_DIFF[i] << " +/- " << N3PI_DIFF_ERR[i] << "\n"
	   << "\tchi2_sum_tmp = " << chi2_sum_tmp << "\n";
      */
      //cout << "bin " << i << "\n";
      
      fit_indx ++;
      
    }

    //cout << chi2_sum_tmp << endl;

      
  }

  //cout << "\n";
  
  chi2_sum_crx3pi = chi2_sum_tmp,
  residul_size_crx3pi = fit_indx;


  f = chi2_sum_tmp;
 
  
}

double GetCrx3piMax_new(double BB, double Mass_V) {

  double crx3pi = 12 * pi * BB / Mass_V / Mass_V * 1e6 * crx_scale;

  /*
  cout << "\n\tBB = " << BB << "\n"
       << "\tM_V = " << mass_omega << " [MeV]\n"
       << "\tcrx3pi = " << crx3pi << " [nb]\n";
  */
  
  return crx3pi;
  
}

double GetBB_new(double sigma_max, double Mass_V) {

  double BB = sigma_max * 1e-6 / crx_scale * Mass_V * Mass_V / 12. / pi;

  return BB;

}

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

Double_t DetectorEvent(Double_t mTrue) {
  // smear by double-gaussian
  if(rnd->Rndm()>frac) {
    return rnd->Gaus(mTrue+smallBias,smallSigma);
  } else {
    return rnd->Gaus(mTrue+wideBias,wideSigma);
  }
}
