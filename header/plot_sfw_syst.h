const TString syst_path = "../../presel_syst_egammamin/";
const TString outputPlot = "../../plot_sfw_syst_egammamin/";
const TString pathType = "sfw2d_path_output.txt";

const TString cut_label = "egammamin";
const TString cut_title = "E^{min}_{clust} [MeV]";

double XLIST[1000], XLIST_ERR[1000];
double YLIST[1000], YLIST_ERR[1000];
double BAND[1000], BAND_ERR[1000];
double UNCORR_ERR[1000];
double X_NORM[1], Y_NORM[1];
double X_ERR_NORM[1], Y_ERR_NORM[1];

const int para_indx = 5;
const TString para_label = "mcrest_sfw";
const TString para_title = "MC others scaling factor";
const TString para_unit = "";

TGraphErrors *gf_syst = new TGraphErrors();
TGraphErrors *gf_band = new TGraphErrors();
TGraphErrors *gf_norm = new TGraphErrors();

double upper_band_size = 7.;
double down_band_size = 6.;

  
TGraphErrors *get_graph_norm(double xlist[], double ylist[], double xlist_err[], double ylist_err[], int length) {

  gf_norm = new TGraphErrors(length, xlist, ylist, xlist_err, ylist_err);

  gf_norm -> SetName("gf_norm");
  gf_norm -> SetMarkerStyle(20);
  gf_norm -> SetMarkerSize(1.1);
  gf_norm -> SetLineColor(kBlue);
  gf_norm -> SetMarkerColor(kBlue);
  gf_norm -> SetLineWidth(2);

  return gf_norm;
  
}

TGraphErrors *get_graph_band(double xlist[], double ylist[], double xlist_err[], double ylist_err[], int length) {

  gf_band = new TGraphErrors(length, xlist, ylist, xlist_err, ylist_err);
  
  gf_band -> SetName("gf_band");
  gf_band -> SetMarkerStyle(20);
  gf_band -> SetMarkerSize(1.1);
  gf_band -> SetLineColor(kBlue);
  gf_band -> SetLineWidth(2);
  gf_band -> SetFillStyle(3002);
  gf_band -> SetFillColor(1);
  
  return gf_band;

}

TGraphErrors *get_graph_syst(double xlist[], double ylist[], double xlist_err[], double ylist_err[], int length) {

  TGraphErrors *gf = new TGraphErrors(length, xlist, ylist, xlist_err, ylist_err);

  gf -> SetTitle("Variation");
  gf -> SetName("gf_syst");
  
  gf -> SetMarkerStyle(20);
  gf -> SetMarkerSize(1.2);
  gf -> SetLineColor(kBlack);
  gf -> SetLineWidth(2);
  gf -> SetMarkerColor(kBlack);
  gf -> GetXaxis() -> SetTitle(cut_title);
  gf -> GetXaxis() -> CenterTitle();
  //gf -> GetXaxis() -> SetRangeUser(-0.5., 0.5.);
  gf -> GetXaxis() -> SetTitleOffset(1.);
  gf -> GetXaxis() -> SetTitleSize(0.06);
  gf -> GetXaxis() -> SetLabelSize(0.05);  

  gf -> GetYaxis() -> SetTitleOffset(1.);
  //gf -> GetYaxis() -> SetRangeUser(Y_NORM[0] - down_band_size * Y_ERR_NORM[0], Y_NORM[0] + upper_band_size * Y_ERR_NORM[0]);
  gf -> GetYaxis() -> SetTitleSize(0.08);
  gf -> GetYaxis() -> SetLabelSize(0.05);  
  gf -> GetYaxis() -> SetTitle(para_title + " " + para_unit);
  gf -> GetYaxis() -> CenterTitle();
  
  //gf -> Draw("AP");
  
  return gf;

}

void fill_uncorr(double LIST_TARGET[], double LIST[], double LIST_ERR[], int length) {

  cout << "Uncorr. error" << endl;
  //cout << "step_nb = " << step_nb << endl;
  //cout << "Y_NORM[0] = " << Y_NORM[0] << " +/- " << Y_ERR_NORM[0] << endl;
  
  double uncorr_err = 0.;
  for (int i = 0; i < length; i++) {

    uncorr_err = TMath::Sqrt(TMath::Abs(LIST_ERR[i] * LIST_ERR[i] - Y_ERR_NORM[0] * Y_ERR_NORM[0]));
    LIST_TARGET[i] = uncorr_err;
    BAND[i] = Y_NORM[0];
    BAND_ERR[i] = Y_ERR_NORM[0];
    cout << i << ", " << LIST[i] << "+/-" << LIST_TARGET[i] << ", " << LIST_ERR[i] << ", " << Y_ERR_NORM[0] << endl;

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
