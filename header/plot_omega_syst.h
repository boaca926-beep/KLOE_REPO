const TString syst_path = "../../result_syst/norm_lumi_nb/";
const TString outputPlot = "../../plot_results/plot_norm_lumi_nb/";

const TString cut_label = "lumi_nb";
const TString cut_title = "#DeltaL_{int} [#deltaL_{int}/L_{int}]";
int err_type = 2;

double XLIST[1000], XLIST_ERR[1000];
double YLIST[1000], YLIST_ERR[1000];
double ZLIST[1000], ZLIST_ERR[1000];
double BAND[1000], BAND_ERR[1000];
double SIGMABAND[1000], SIGMABAND_ERR[1000];
double UNCORR_ERR[1000];
double X_NORM[1], Y_NORM[1], Z_NORM[1];
double X_ERR_NORM[1], Y_ERR_NORM[1], Z_ERR_NORM[1];
double SYST_ERR[2];

TTree* TRESULT = new TTree("TRESULT", "recreate");
  
const int para_indx = 2;
const TString para_label = "BB";
const TString para_title = "B_{ee}B_{3#pi}";
const TString para_unit = "";

TGraphErrors *gf_syst = new TGraphErrors();
TGraphErrors *gf_band = new TGraphErrors();
TGraphErrors *gf_sigmaband_plus = new TGraphErrors();
TGraphErrors *gf_sigmaband_nega = new TGraphErrors();
TGraphErrors *gf_norm = new TGraphErrors();
TGraphErrors *gf_Z = new TGraphErrors();
TGraphErrors *gf_Z_norm = new TGraphErrors();

const double nb_sigma = 2.73;
double band_limit = 0.;


X_ERR_NORM[0] = 0.;

Z_NORM[0] = 0.;
Z_ERR_NORM[0] = nb_sigma;

int norm_indx = -1;
int step1_indx = -1;
int step2_indx = -1;

double step1_diff = 0., step2_diff = 0.;
double step1_Z = 0., step2_Z = 0.;
//double nega_err = 0., plus_err = 0.;


void fill_uncorr(double LIST_TARGET[], double LIST[], double LIST_ERR[], int length) {

  cout << "Uncorr. error" << endl;
  //cout << "step_nb = " << step_nb << endl;
  //cout << "Y_NORM[0] = " << Y_NORM[0] << " +/- " << Y_ERR_NORM[0] << endl;
  
  double uncorr_err = 0.;
  double Z_max = 0.;
  double Y_diff = 0.;

  for (int i = 0; i < length; i++) {

    uncorr_err = TMath::Sqrt(TMath::Abs(LIST_ERR[i] * LIST_ERR[i] - Y_ERR_NORM[0] * Y_ERR_NORM[0]));
    LIST_TARGET[i] = uncorr_err;
    BAND[i] = Y_NORM[0];
    BAND_ERR[i] = Y_ERR_NORM[0];

    // 2sigma band
    SIGMABAND[i] = 0.;
    SIGMABAND_ERR[i] = nb_sigma; //nb_sigma;

    // Z band
    Y_diff = YLIST[i] - Y_NORM[0];
    
    if (YLIST[i] == Y_NORM[0]) {
      ZLIST[i] = 0.;
    }
    else {
      ZLIST[i] = Y_diff / uncorr_err;
    }

    if (TMath::Abs(ZLIST[i]) > Z_max) {
      Z_max = TMath::Abs(ZLIST[i]);
    }
    
    ZLIST_ERR[i] = 0.;

    // determine Z for +/- 1 resolution variation, ang type II error
    if (i == step1_indx) {
      if (TMath::Abs(ZLIST[i]) > nb_sigma) {
	step1_diff = Y_diff;
	step1_Z = ZLIST[i];
	//cout << i << endl;
      }
      else {
	step1_diff = 0.;
      }
    }

    if (i == step2_indx) {
      if (TMath::Abs(ZLIST[i]) > nb_sigma) {
	step2_diff = Y_diff;
	step2_Z = ZLIST[i];
	//cout << i << endl;
      }
      else {
	step2_diff = 0.;
      }
    }
   

    //
    cout << i << ", sigma_delta = " << LIST[i] << "+/-" << LIST_TARGET[i] << ", Y_Err = " << LIST_ERR[i] << ", YNORM_Err = " << Y_ERR_NORM[0] << ", Diff = " << Y_diff << ", Z value = " << ZLIST[i] << ", Z_max = " << Z_max << endl;

  }


  if (Z_max > nb_sigma) {
    band_limit = Z_max * 1.5;
  }
  else {
    band_limit = nb_sigma * 1.5;
  }
  
}


double get_waverage(double LIST[], double LIST_ERR[], int length) {

  double average = 0.; // weighted average
  double w = 0.;
  double w_sum = 0.;
  double wy_sum = 0.;
  
  for (int i = 0; i < length; i++) {

    w = 1 / LIST_ERR[i] / LIST_ERR[i];
    w_sum += w;
    wy_sum += LIST[i] * w;
  }

  average = wy_sum / w_sum;
  
  return average;
  
}

double get_wy_sdv(double LIST[], double LIST_ERR[], int length, double waverage) {

  double w = 0.;
  double wy2_sum = 0.;
  double w_sum = 0.;
  double w2_sum = 0.;
  double N_eff = 0.;
  double corr_factor = 0.;
  double wy_sdv = 0.;
  
  for (int i = 0; i < length; i ++) {
    w = 1 / LIST_ERR[i] / LIST_ERR[i];
    w_sum += w;
    w2_sum += w * w;
    wy2_sum += (LIST[i] - waverage) * (LIST[i] - waverage) * w;
    //cout << LIST[i] << "+/-" << LIST_ERR[i] << endl;
  }

  N_eff = w_sum * w_sum / w2_sum;
  corr_factor = N_eff / (N_eff - 1); // corretion factor
  
  return wy_sdv = TMath::Sqrt(wy2_sum / w_sum * corr_factor);
  
}

void getnorm(double NORM[], double ERR_NORM[], double YLIST[], double YLIST_ERR[], int norm_indx, int length) {

  for (int i = 0; i < length; i++) {

    if (i == norm_indx) {
      NORM[0] = YLIST[i];
      ERR_NORM[0] = YLIST_ERR[i];
    }
    
  }
  

}

double getmaxoftwo(double a, double b) {

  double max = 0.;

  if (a >= b) {

    max = a;
    
  }
  else {

    max = b;
    
  }

  return max;
  
}

double getminoftwo(double a, double b) {

  double min = 0.;

  if (a <= b) {

    min = a;
    
  }
  else {

    min = b;
    
  }

  return min;
  
}

void get_syst_errII(double err1, double err2){

  double max = getmaxoftwo(err1, err2);
  double min = getminoftwo(err1, err2);
  double plus_err = 0., nega_err = 0.;

    
  if (err1 * err2 < 0.) {

    if (err1 == max) {
      plus_err = err1;
      nega_err = err2;
    }
    else {
      plus_err = err2;
      nega_err = err1;
    }
    //cout << max << endl;
    
  }
  else if (err1 * err2 > 0.) {

    //cout << "min = " << min << ", max = " << max << endl;
    
    if (err1 < 0) {
      plus_err = 0.;
      nega_err = min;
    }
    else {
      plus_err = max;
      nega_err = 0.;
    }
    
  }
  else if (err1 == 0. && err2 != 0.) {

    if (err2 > 0.) {
       plus_err = err2;
       nega_err = 0.;
    }
    else {
      plus_err = 0.;
      nega_err = err2;
    }
    
  }
  else if (err2 == 0. && err1 != 0.) {

    if (err1 > 0.) {
       plus_err = err1;
       nega_err = 0.;
    }
    else {
      plus_err = 0.;
      nega_err = err1;
    }
    
  }
  else if (err1 == 0. && err2 == 0.) {

    plus_err = 0.;
    nega_err = 0.;
    
  }
  else {
    cout << "New case !!!" << endl;
  }
  
  
  SYST_ERR[0] = nega_err;
  SYST_ERR[1] = plus_err;
  
  

}

