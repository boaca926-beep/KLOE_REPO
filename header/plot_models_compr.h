TGraphErrors *n3pi_obs_f1, *n3pi_obs_f2;
TGraphErrors *n3pi_fit_f1, *n3pi_fit_f2;
TGraphErrors *n3pi_diff_f1, *n3pi_diff_f2;

TTree *TCRX3PI_f1 = nullptr;
TTree *TCRX3PI_f2 = nullptr;

TH1D *hdata;

double Delta_m3pi = 0.;
double OMEGA_PARA_F1[3], OMEGA_PARA_ERR_F1[3];
double OMEGA_PARA_F2[3], OMEGA_PARA_ERR_F2[3];
const TString PARA_LABEL[3] = {"mass", "width", "BB"};
  
double chi2_sum_f1 = 0., Lumi_int_fit_f1 = 0.;
double chi2_sum_f2 = 0., Lumi_int_fit_f2 = 0.;

int ndf_f1 = 0, ndf_f2 = 0;

const TString file_type1 = "isrlumi_norm";
const TString file_type2 = "norm";

const TString legend_type1 = "Exact";
const TString legend_type2 = "Approx";

const TString outputFolder = "../../plot_models_compr/";
