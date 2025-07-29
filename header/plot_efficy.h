const TString input_folder = "../../efficy_evtcls";
const TString systType = "evtcls";

const double mass_min = 760., mass_max = 800.;
cout << "mass_min = " << mass_min << ", mass_max = " << mass_max << endl;


TArrayD get_w_ratio(TGraphErrors *gf, const double mass_min, const double mass_max, const double Delta_m3pi) {

  // Weighted average of gf
  double *x_gf = gf -> GetX();
  double *y_gf = gf -> GetY();
  double *y_gf_err = gf -> GetEY();

  double ratio_tmp = 0.;
  double w_y = 0., w_sum = 0., wy_sum = 0.;
  double err_tmp = 0.;
  
  int count = 0;
  
  for (int i = 0; i < gf -> GetN(); i ++) {
  //for (int i = 0; i < 3; i ++) {

    if (x_gf[i] >= mass_min - Delta_m3pi && x_gf[i] <= mass_max + Delta_m3pi) {
      count ++;
      ratio_tmp = y_gf[i];
      err_tmp = y_gf_err[i];
      w_y = 1. / (err_tmp * err_tmp);
      w_sum += w_y;
      wy_sum += ratio_tmp * w_y;

      //cout << "mass = " << x_gf[i] << ", ratio = " << ratio_tmp << "+/-" << err_tmp  << endl;
   
      //cout << "ratio = " << ratio_tmp << "+/-" << err_tmp << ", w_y =" << w_y << ", wy_sum = " << wy_sum << ", ratio * w_ y = " << ratio_tmp * w_y << ", w_sum = " << w_sum << endl;
    }

  }

  TArrayD W(2);
  
  double ratio_average = wy_sum / w_sum;
  double sigma_ratio = 1. / TMath::Sqrt(w_sum);
  //cout << "ratio_average = " << ratio_average << "+/-" << sigma_ratio << endl;

  W[0] = ratio_average;
  W[1] = sigma_ratio;

  return W;

}

TCanvas *plotting_efficy(const TString cv_title, const TString cv_nm, TGraphErrors *gf_sig, TGraphErrors *gf_ufo, TGraphErrors *gf_ratio, TGraphErrors *gf_ratio_corr, const double ymax, const TString Note){

  double x1 = 0., y1 = 0.;
  double x2 = 0., y2 = 0.;
  
  gf_sig -> GetPoint(0, x1, y1);
  gf_sig -> GetPoint(1, x2, y2);

  double Delta_m3pi = x2 - x1;
  //double ymax = gf_ufo -> GetMaximum() * 1.2;
  
  cout << "Delta_m3pi = " << Delta_m3pi << ", ymax = " << ymax << endl;

  cout << "w_ratio" << endl; 
  TArrayD w_ratio = get_w_ratio(gf_ratio, mass_min, mass_max, Delta_m3pi);

  cout << "w_efficy_sig" << endl; 
  TArrayD w_efficy_sig = get_w_ratio(gf_sig, mass_min, mass_max, Delta_m3pi);

  cout << "w_efficy_ufo" << endl; 
  TArrayD w_efficy_ufo = get_w_ratio(gf_ufo, mass_min, mass_max, Delta_m3pi);

  cout << "mass range [" << mass_min << ", " << mass_max << "] MeV/c^2\n"
       << "ratio average mean = " << w_ratio[0] << "+/-" << w_ratio[1] << endl;
  cout << "efficy_sig average mean = " << w_efficy_sig[0] << "+/-" << w_efficy_sig[1] << endl;
  cout << "efficy_ufo average mean = " << w_efficy_ufo[0] << "+/-" << w_efficy_ufo[1] << endl;
  cout << "w_efficy_ufo[0] / w_efficy_sig[0] = " << w_efficy_ufo[0] / w_efficy_sig[0] << endl;
  cout << "Full mass range: ratio = " << 0.561752 / 0.656811 << endl;
  
  TLine *line = new TLine(mass_min - Delta_m3pi / 2., w_ratio[0], mass_max, w_ratio[0]);
  line -> SetLineColor(kRed);
  line -> SetLineWidth(2);

  TPaveText *pt34 = new TPaveText(0.5, 0.85, 0.7, 0.8, "NDC");
  PteAttr(pt34);
  pt34 -> SetTextSize(0.07);
  pt34 -> AddText(Note);

  // plot
  TCanvas *cv = new TCanvas(cv_title, cv_nm, 1200, 800);

  TPad *p2 = new TPad("p2", "p2", 0., 0., 1., 0.35);
  p2 -> Draw();
  p2 -> SetBottomMargin(0.25);
  p2 -> SetTopMargin(0.04);
  p2 -> SetLeftMargin(0.1);
  //p2 -> SetGrid();

  TPad *p1 = new TPad("p1", "p1", 0., 0.35, 1., 1.);
  p1 -> Draw();
  p1 -> SetBottomMargin(0.02);//0.007
  p1 -> SetLeftMargin(0.1);
  p1 -> SetGrid();

  p1 -> cd();

  gf_ufo -> GetYaxis() -> SetNdivisions(512);
  gf_ufo -> GetYaxis() -> SetTitleFont(43);
  gf_ufo -> GetYaxis() -> SetRangeUser(0., ymax * 1.4);
  //gf_ufo -> GetYaxis() -> SetTitle(TString::Format("Efficiency (#tilde{#varepsilon})/[%0.2f MeV/c^{2}]", Delta_m3pi));
  gf_ufo -> GetYaxis() -> SetTitle("Efficiency #tilde{#varepsilon}");
  gf_ufo -> GetYaxis() -> SetTitleSize(33);
  gf_ufo -> GetYaxis() -> SetTitleOffset(1.5);
  gf_ufo -> GetYaxis() -> SetLabelSize(0.05);
  gf_ufo -> GetYaxis() -> CenterTitle();

  gf_ufo -> GetXaxis() -> SetTitle("M_{3#pi} [MeV/c^{2}]");
  gf_ufo -> GetXaxis() -> SetTitleOffset(1.1);
  gf_ufo -> GetXaxis() -> SetTitleSize(0.04);
  gf_ufo -> GetXaxis() -> SetLabelSize(0.04);
  gf_ufo -> GetXaxis() -> SetLabelOffset(4);
  gf_ufo -> GetXaxis() -> SetRangeUser(mass_min, mass_max);
  gf_ufo -> GetXaxis() -> CenterTitle();
  
  //gf_ufo -> SetLineColor(1);
  gf_ufo -> SetLineWidth(2);
  
  //gf_sig -> SetLineColor(46);
  //gf_sig -> SetMarkerColor(4);
  gf_sig -> SetLineWidth(2);
  
  gf_ufo -> Draw("APZ");
  gf_sig -> Draw("PZ");
  pt34 -> Draw("Same");
  
  p1 -> Update();

  //
  TLegend * legd_cv_p1 = new TLegend(0.2, 0.75, 0.5, 0.85);
  
  SetLegend(legd_cv_p1);
  legd_cv_p1 -> SetTextSize(0.07);
  legd_cv_p1 -> SetNColumns(2);
  
  legd_cv_p1 -> AddEntry(gf_ufo, "#tilde{#varepsilon}_{UFO}", "lep");
  legd_cv_p1 -> AddEntry(gf_sig, "#tilde{#varepsilon}_{sig}", "lep");
  
  legd_cv_p1 -> Draw("Same");

  p2 -> cd();

  //gf_ratio -> SetLineColor(0);
  gf_ratio -> GetYaxis() -> SetNdivisions(505);
  gf_ratio -> GetYaxis() -> SetTitleSize(33);
  gf_ratio -> GetYaxis() -> SetTitleFont(43);
  gf_ratio -> GetYaxis() -> SetTitleOffset(1.4);
  ////gf_ratio -> GetYaxis() -> SetLabelFont(43); // Absolute front size in pixel (precision 3)
  gf_ratio -> GetYaxis() -> SetLabelSize(0.1);
  //gf_ratio -> GetYaxis() -> SetRangeUser(0.5, 1.2);
  //gf_ratio -> GetYaxis() -> SetRangeUser(0.5, 1.7);
  gf_ratio -> GetYaxis() -> SetRangeUser(0.9, 1.7);
  
  //gf_ratio -> GetYaxis() -> SetRangeUser(0., gf_ratio -> GetMaximum() * 1.2);
  gf_ratio -> GetYaxis() -> CenterTitle();
  gf_ratio -> GetYaxis() -> SetTitle("#tilde{#varepsilon}_{UFO}/#tilde{#varepsilon}_{sig}");

  gf_ratio -> GetXaxis() -> SetRangeUser(mass_min, mass_max);
  gf_ratio -> GetXaxis() -> SetTitleOffset(.8);
  gf_ratio -> GetXaxis() -> SetLabelSize(0.1);
  gf_ratio -> GetXaxis() -> SetTitleSize(0.13);
  gf_ratio -> GetXaxis() -> SetTitle("M_{3#pi} [MeV/c^{2}]");
  gf_ratio -> GetXaxis() -> CenterTitle();

  gf_ratio -> SetLineColor(kBlack);
  gf_ratio -> SetLineWidth(2);
  
  gf_ratio -> Draw("APZ");
  gf_ratio_corr -> Draw("P");
  //line -> Draw("Same");
  
  p2 -> SetGrid();
  
  TLegend *legd_cv_p2 = new TLegend(0.2, 0.65, 0.4, 0.9);
  SetLegend(legd_cv_p2);
  legd_cv_p2 -> SetTextSize(0.1);
  legd_cv_p2 -> SetNColumns(1);
  
  legd_cv_p2 -> AddEntry(gf_ratio, "#tilde{#varepsilon}_{sig}/#tilde{#varepsilon}_{UFO}", "E");
  
  //legd_cv_p2 -> Draw("Same");
  
  return cv;

}

TArrayD get_gf_max(TGraphErrors *gf) {

  // Weighted average of gf
  double *x_gf = gf -> GetX();
  double *y_gf = gf -> GetY();
  double *y_gf_err = gf -> GetEY();

  int count = 0;
  double max_value = 0.;
  double value_tmp = 0.;

  double x1 = 0., x2 = 0.;
  double y1 = 0., y2 = 0.;
  
  gf -> GetPoint(0, x1, y1);
  gf -> GetPoint(1, x2, y2);

  double Delta_m3pi = x2 - x1;

  for (int i = 0; i < gf -> GetN(); i ++) {

    if (x_gf[i] >= mass_min - Delta_m3pi && x_gf[i] <= mass_max + Delta_m3pi) {
      
      count ++;
      
      value_tmp = y_gf[i];
      
      if (value_tmp > max_value) {
	max_value = value_tmp;
      }
      
      //cout << value_tmp << ", max_value = " << max_value << endl;

    }
    
  }

  //cout << "max_value = " << max_value << endl;

  TArrayD W(1);
  
  W[0] = max_value;
  
  return W;

}

TCanvas *plotting_nb(const TString cv_title, const TString cv_nm, TGraphErrors *gf_nb_sel, TGraphErrors *gf_nb_evtcls, TGraphErrors *gf_efficy, const TString y_title, const double ymax, const TString Note){

  double x1 = 0., y1 = 0.;
  double x2 = 0., y2 = 0.;
  
  gf_nb_sel -> GetPoint(0, x1, y1);
  gf_nb_sel -> GetPoint(1, x2, y2);

  double Delta_m3pi = x2 - x1;

  cout << "Delta_m3pi = " << Delta_m3pi << ", ymax = " << ymax << endl;

  const double mass_min = 760., mass_max = 800.;
  cout << "mass_min = " << mass_min << ", mass_max = " << mass_max << endl;
  
  TPaveText *pt34 = new TPaveText(0.5, 0.85, 0.7, 0.8, "NDC");
  PteAttr(pt34);
  pt34 -> SetTextSize(0.07);
  pt34 -> AddText(Note);

  TCanvas *cv = new TCanvas(cv_title, cv_nm, 1200, 800);
  
  TPad *p2 = new TPad("p2", "p2", 0., 0., 1., 0.35);
  p2 -> Draw();
  p2 -> SetBottomMargin(0.25);
  p2 -> SetTopMargin(0.04);
  p2 -> SetLeftMargin(0.1);
  //p2 -> SetGrid();

  TPad *p1 = new TPad("p1", "p1", 0., 0.35, 1., 1.);
  p1 -> Draw();
  p1 -> SetBottomMargin(0.02);//0.007
  p1 -> SetLeftMargin(0.1);
  p1 -> cd();

  gf_nb_sel -> GetYaxis() -> SetNdivisions(512);
  gf_nb_sel -> GetYaxis() -> SetTitleFont(43);
  gf_nb_sel -> GetYaxis() -> SetRangeUser(0., ymax *  1.5);
  gf_nb_sel -> GetYaxis() -> SetTitle(TString::Format("Events/[%0.2f MeV/c^{2}]", Delta_m3pi));
  //gf_nb_sel -> GetYaxis() -> SetTitle("");
  gf_nb_sel -> GetYaxis() -> SetTitleSize(33);
  gf_nb_sel -> GetYaxis() -> SetTitleOffset(1.5);
  gf_nb_sel -> GetYaxis() -> SetLabelSize(0.05);
  gf_nb_sel -> GetYaxis() -> CenterTitle();

  gf_nb_sel -> GetXaxis() -> SetTitle("M_{3#pi} [MeV/c^{2}]");
  gf_nb_sel -> GetXaxis() -> SetTitleOffset(1.1);
  gf_nb_sel -> GetXaxis() -> SetTitleSize(0.04);
  gf_nb_sel -> GetXaxis() -> SetLabelSize(0.04);
  gf_nb_sel -> GetXaxis() -> SetLabelOffset(4);
  gf_nb_sel -> GetXaxis() -> SetRangeUser(mass_min, mass_max);
  gf_nb_sel -> GetXaxis() -> CenterTitle();
  
  //gf_nb_sel -> SetLineColor(1);
  gf_nb_sel -> SetLineWidth(2);
  
  //gf_nb_sel -> SetLineColor(46);
  //gf_nb_sel -> SetMarkerColor(4);
  gf_nb_sel -> SetLineWidth(2);
  
  gf_nb_sel -> Draw("APZ");
  gf_nb_evtcls -> Draw("PZ");
  pt34 -> Draw("Same");
  
  //
  TLegend * legd_cv_p1 = new TLegend(0.15, 0.7, 0.5, 0.9);
  
  SetLegend(legd_cv_p1);
  legd_cv_p1 -> SetTextSize(0.06);
  legd_cv_p1 -> SetNColumns(1);
  
  legd_cv_p1 -> AddEntry(gf_nb_evtcls, "After KSL", "lep");
  legd_cv_p1 -> AddEntry(gf_nb_sel, "Before KSL", "lep");
  
  legd_cv_p1 -> Draw("Same");

  p2 -> cd();

  //gf_efficy -> SetLineColor(0);
  gf_efficy -> GetYaxis() -> SetNdivisions(505);
  gf_efficy -> GetYaxis() -> SetTitleSize(33);
  gf_efficy -> GetYaxis() -> SetTitleFont(43);
  gf_efficy -> GetYaxis() -> SetTitleOffset(1.4);
  ////gf_efficy -> GetYaxis() -> SetLabelFont(43); // Absolute front size in pixel (precision 3)
  gf_efficy -> GetYaxis() -> SetLabelSize(0.1);
  //gf_efficy -> GetYaxis() -> SetRangeUser(0.2, 0.8);
  gf_efficy -> GetYaxis() -> SetRangeUser(0., 1.);
  
  //gf_efficy -> GetYaxis() -> SetRangeUser(0., gf_efficy -> GetMaximum() * 1.2);
  gf_efficy -> GetYaxis() -> CenterTitle();
  gf_efficy -> GetYaxis() -> SetTitle(y_title);

  gf_efficy -> GetXaxis() -> SetRangeUser(mass_min, mass_max);
  gf_efficy -> GetXaxis() -> SetTitleOffset(.8);
  gf_efficy -> GetXaxis() -> SetLabelSize(0.1);
  gf_efficy -> GetXaxis() -> SetTitleSize(0.13);
  gf_efficy -> GetXaxis() -> SetTitle("M_{3#pi} [MeV/c^{2}]");
  gf_efficy -> GetXaxis() -> CenterTitle();

  gf_efficy -> SetLineColor(kBlack);
  gf_efficy -> SetLineWidth(2);
  
  gf_efficy -> Draw("APZ");
  p2 -> SetGrid();
  
  TLegend *legd_cv_p2 = new TLegend(0.2, 0.65, 0.4, 0.9);
  SetLegend(legd_cv_p2);
  legd_cv_p2 -> SetTextSize(0.1);
  legd_cv_p2 -> SetNColumns(1);
  
  legd_cv_p2 -> AddEntry(gf_efficy, "#tilde{#varepsilon}_{sig}/#tilde{#varepsilon}_{UFO}", "E");
  
  //legd_cv_p2 -> Draw("Same");

  
  return cv;

}

double ratioErr(double a, double sigma_a, double b, double sigma_b) {// ratio = a / b

  if (b == 0) {

    cout << "Division by zero!" << endl;
    return 0;

  }

  double c = a / b;
  double rela_err = TMath::Sqrt(TMath::Power(sigma_a / a, 2) + TMath::Power(sigma_b / b, 2));
  double sigma_c = c * rela_err;

  //cout << TMath::Power(sigma_a / a, 2) + TMath::Power(sigma_b / b, 2) << ", a = " << a << ", b = " << b << endl;
   
  return sigma_c;
  
}

double get_ratio(double a, double b) {// ratio = a / b

  if (b == 0) {

    cout << "Division by zero!" << endl;
    return 0;

  }

  double c = a / b;

  return c;
  
}

double get_efficy_ratio_err(double *PARA, double *PARA_ERR, double x) {

  double ratio_err = 0.;
  double ratio_err2 = 0.;

  double a_err = PARA_ERR[0];
  double b_err = PARA_ERR[1];
  double c_err = PARA_ERR[2];
  
  double term1 = TMath::Power(x, 4) * TMath::Power(a_err, 2);
  double term2 = TMath::Power(x, 2) * TMath::Power(b_err, 2);
  double term3 = TMath::Power(c_err, 2);

  ratio_err2 = term1 + term2 + term3;
  ratio_err = TMath::Sqrt(ratio_err2);
  
  cout << "a = " << PARA[0] << "+/-" << PARA_ERR[0] << "\n"
       << "b = " << PARA[1] << "+/-" << PARA_ERR[1] << "\n"
       << "c = " << PARA[2] << "+/-" << PARA_ERR[2] << "\n"
       << "x = " << x << "\n"
       << "ratio_err = " << ratio_err << endl;
 
  return ratio_err;
}
