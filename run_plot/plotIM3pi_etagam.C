#include "fitfun.h"
#include "hist.h"

int plotIM3pi_etagam() {

  //gErrorIgnoreLevel = kError;
  TGaxis::SetMaxDigits(4);
  gStyle->SetOptTitle(0);
  gStyle->SetStatBorderSize(0);
  //gStyle->SetOptFit(11);
  //gStyle->SetStatX(0.89);
  //gStyle->SetStatY(0.89);
  gStyle->SetFitFormat("6.4g");
  
  //gStyle->SetOptStat(0);

  cout << "plot IM3pi, input file: IM3pi_etagam, IM3pi_rhopi" << endl;

  TFile *intree = new TFile("prompt_output/IM3pi.root"); // IM3pi.root
  
  TIter next_tree(intree -> GetListOfKeys());
  
  TString objnm_tree, classnm_tree;
  
  int i = 0;

  TKey *key;
  
  while ( (key = (TKey *) next_tree() ) ) {// start tree while lop
    
    i ++;
    
    objnm_tree   =  key -> GetName();
    classnm_tree = key -> GetClassName();
    key -> GetSeekKey();
    
    cout << " tree" << i << ": classnm = " << classnm_tree << ", objnm = " << objnm_tree << endl;
    
  }

  // get histo
  TH1D* h1d = (TH1D*)intree -> Get("IM3pi_etagam"); //IM3pi_etagam
  
  //const double bin_size = h1d -> GetNbinsX(); 
  //double norm_hist = h1d -> Integral(1, bin_size);

  //const double bin_size1 = h1d_resid -> GetNbinsX(); 
  //double norm_hist1 = h1d_resid -> Integral(1, bin_size1);

  //h1d -> Scale(1. / norm_hist);
  //h1d_resid -> Scale(1. / norm_hist1);

  // double Gaussian fit
  TString hist_nm = h1d -> GetName();
  double xmin_tmp = h1d -> GetXaxis() -> GetXmin();
  double xmax_tmp = h1d -> GetXaxis() -> GetXmax();

  cout << "xmin = " << xmin_tmp << ", xmax = " << xmax_tmp << endl;
  
  const double fit_range1 = 542.;
  const double fit_range2 = 555.;

  const int npar = 6, npar_sub = 3;
  double rms = h1d -> GetRMS();
  double mean = h1d -> GetMean();
  double peak = h1d -> GetMaximum();
  double bin_width = getbinwidth(h1d);
  cout << "peak = " << peak << endl;
  
  double fitpara[npar] = {peak, mean, rms, 0.01 * peak, mean, 3. * rms};
  double fitgfpar[npar], fitgfpar1[npar_sub], fitgfpar2[npar_sub];
  TF1 *fitfun = new TF1("fitfun", fun_double, xmin_tmp, xmax_tmp, npar);
  fitfun -> SetParNames("p1","p2","p3","p4","p5","p6");
  fitfun -> SetParameters(fitpara);
  //fitfun -> SetParLimits(4,-.5,.5);
  fitfun -> SetLineWidth(2);
  fitfun -> SetNpx(5000);

  // fitting
  TFitResultPtr r = h1d -> Fit(fitfun, "L0Q", " ", fit_range1, fit_range2);
  
  // gauss1
  TF1 *gauss1_fit = new TF1("gauss1_fit", gauss1d, fit_range1, fit_range2, npar_sub);
  fitgfpar1[0] = fitfun -> GetParameter(0); // peak1
  fitgfpar1[1] = fitfun -> GetParameter(1); // mean1
  fitgfpar1[2] = fitfun -> GetParameter(2); // width1
  
  gauss1_fit -> SetParameters(fitgfpar1);
  gauss1_fit -> SetLineColor(kBlack);
  gauss1_fit -> SetLineWidth(2);
  gauss1_fit -> SetLineStyle(1);
  gauss1_fit -> SetNpx(5000);

  // gauss2
  TF1 *gauss2_fit = new TF1("gauss2_fit", gauss1d, fit_range1, fit_range2, npar_sub);
  fitgfpar2[0] = fitfun -> GetParameter(3); // peak2
  fitgfpar2[1] = fitfun -> GetParameter(4); // mean2
  fitgfpar2[2] = fitfun -> GetParameter(5); // width2
  
  gauss2_fit -> SetParameters(fitgfpar2);
  gauss2_fit -> SetLineColor(kBlue);
  gauss2_fit -> SetLineStyle(3);
  gauss2_fit -> SetLineWidth(2);
  gauss2_fit -> SetNpx(5000); 
  
  // gauss_sum
  TF1 *gauss_sum = new TF1("gauss_sum", fun_double, fit_range1, fit_range2, npar);
  fitgfpar[0] = fitfun -> GetParameter(0); // peak1
  fitgfpar[1] = fitfun -> GetParameter(1); // mean1
  fitgfpar[2] = fitfun -> GetParameter(2); // width1
  fitgfpar[3] = fitfun -> GetParameter(3); // peak2
  fitgfpar[4] = fitfun -> GetParameter(4); // mean2
  fitgfpar[5] = fitfun -> GetParameter(5); // width2
  
  gauss_sum -> SetParameters(fitgfpar);
  gauss_sum -> SetLineColor(kBlack);
  gauss_sum -> SetLineStyle(1);
  gauss_sum -> SetLineWidth(2);
  gauss_sum -> SetNpx(5000);

  // mass difference in %
  const double Meta_pdg = 547.862; // eta mass pdg MeV/c^2
  const double Meta_err_pdg = 0.017; // eta mass err pdg MeV/c^2
  
  const double Momega_pdg = 782.66; // omega mass pdg MeV/c^2
  const double Momega_err_pdg = 0.13; // omega mass err pdg MeV/c^2

  //const double Momega_kloe = 782.732; // omega mass KLOE analysis, dissertation
  //const double Momega_err_kloe = 3.76286e-2; // omega mass err KLOE analysis, dissertation

  const double Momega_kloe = 782.73; // omega mass KLOE analysis, vmd fit
  const double Momega_err_kloe = 1e-2; // omega mass err KLOE analysis, vmd fit
  
  double Meta_fit = fitfun -> GetParameter(1); // fit eta mass
  double Meta_err_fit = fitfun -> GetParError(1); // fit eta mass err
  
  double Meta_diff = TMath::Abs(Meta_pdg - Meta_fit) / Meta_fit * 100.;
  double Momega_diff = TMath::Abs(Momega_pdg - Momega_kloe) / Momega_kloe * 100.;

  double sigma_Meta = TMath::Sqrt(Meta_err_pdg * Meta_err_pdg + Meta_err_fit * Meta_err_fit);
  
  double signif_Meta = TMath::Abs(Meta_pdg - Meta_fit) / sigma_Meta;

  double signif_Momega = TMath::Abs(Momega_pdg - Momega_kloe) / TMath::Sqrt(Momega_err_kloe * Momega_err_kloe + Momega_err_pdg * Momega_err_pdg);
  
    
  //cout << TMath::Abs(Meta_pdg - Meta_fit) << " " << TMath::Sqrt(Meta_fit * Meta_fit + Meta_err_fit * Meta_err_fit) << endl;
  
  cout << "eta mass [MeV/c^2]: pdg = " << Meta_pdg << "+/-" << Meta_err_pdg << ", fit = " << Meta_fit << "+/-" << Meta_err_fit << ", diff. = " << Form("%0.3f", Meta_diff) << "%, signif = " << signif_Meta << " sigma, sigma = " << sigma_Meta << "\n"
       << "omega mass [MeV/c^2]: pdg = " << Momega_pdg << "+/-" << Momega_err_pdg << ", kloe = " << Momega_kloe << "+/-" << Momega_err_kloe << ", diff = " << Form("%0.3f", Momega_diff) << "%, signif = " << signif_Momega << " sigma \n";
  
  // plot
  const TString cv_type = "etagam";

  //
  TPaveText *pt1 = new TPaveText(0.17, 0.82, 0.8, 0.83, "NDC");
  
  pt1 -> SetTextSize(0.06);
  pt1 -> SetFillColor(0);
  pt1 -> SetTextAlign(12);

  TPaveText *pt2 = new TPaveText(0.2, 0.72, 0.5, 0.73, "NDC");
  
  pt2 -> SetTextSize(0.10);
  pt2 -> SetFillColor(0);
  pt2 -> SetTextAlign(12);

  pt1 -> AddText(Form("M_{#eta}=%0.2f#pm%0.2f MeV/c^{2}", fitfun -> GetParameter(1), fitfun -> GetParError(1)));
  
  pt2 -> AddText("(a)");

  TCanvas *cv = new TCanvas("cv_" + cv_type, "cv_" + cv_type, 0, 0, 700, 700);
  cv -> SetBottomMargin(0.15);//0.007
  cv -> SetLeftMargin(0.15);

  const double ymax = h1d -> GetMaximum();
  xmin_tmp = h1d -> GetXaxis() -> GetXmin();
  xmax_tmp = h1d -> GetXaxis() -> GetXmax();

  h1d -> GetYaxis() -> SetTitle("Events");
  h1d -> GetYaxis() -> SetTitleOffset(1.3);
  h1d -> GetYaxis() -> CenterTitle();
  h1d -> GetYaxis() -> SetLabelSize(0.04);
  h1d -> GetYaxis() -> SetTitleSize(0.06);
  h1d -> GetYaxis() -> SetRangeUser(0.,  1.2 * ymax);
  
  h1d -> GetXaxis() -> SetNdivisions(505);
  h1d -> GetXaxis() -> SetTitle("M_{3#pi} [MeV/c^{2}]");
  h1d -> GetXaxis() -> SetRangeUser(xmin_tmp, xmax_tmp);
  h1d -> GetXaxis() -> SetTitleOffset(1.);
  h1d -> GetXaxis() -> CenterTitle();
  h1d -> GetXaxis() -> SetLabelSize(0.05);
  h1d -> GetXaxis() -> SetTitleSize(0.06);

  h1d -> Draw();
  //gauss1_fit -> Draw("Same");
  //gauss2_fit -> Draw("Same");
  gauss_sum -> Draw("Same");
  pt1 -> Draw("Same");
  //pt2 -> Draw("Same");

  /*
  const double fit_range3 = -.1;
  const double fit_range4 = .1;
  const double ymax_resid = h1d_resid -> GetMaximum();
  const double xmin_resid = h1d_resid -> GetXaxis() -> GetXmin();
  const double xmax_resid = h1d_resid -> GetXaxis() -> GetXmax();
  
  double rms = h1d_resid -> GetRMS();
  double mean = h1d_resid -> GetMean();
  double peak = h1d_resid -> GetMaximum();
  bin_width = getbinwidth(h1d_resid);
  fitpara[0] = peak;
  fitpara[1] = mean;
  fitpara[2] = rms;
  fitpara[3] = 0.01 * peak;
  fitpara[4] = mean;
  fitpara[5] = 3. * rms;
  TF1 *fitfun_resid = new TF1("fitfun_resid", fun_double, xmin_resid, xmax_resid, npar);
  fitfun_resid -> SetParNames("p1","p2","p3","p4","p5","p6");
  fitfun_resid -> SetParameters(fitpara);
  //fitfun_resid -> SetParLimits(4,-.5,.5);
  fitfun_resid -> SetLineWidth(2);
  fitfun_resid -> SetNpx(5000);

  // fitting
  TFitResultPtr r1 = h1d_resid -> Fit(fitfun_resid, "L", " ", fit_range3, fit_range4);

  TF1 *gauss_sum_resid = new TF1("gauss_sum_resid", fun_double, fit_range3, fit_range4, npar);
  fitgfpar[0] = fitfun -> GetParameter(0); // peak1
  fitgfpar[1] = fitfun -> GetParameter(1); // mean1
  fitgfpar[2] = fitfun -> GetParameter(2); // width1
  fitgfpar[3] = fitfun -> GetParameter(3); // peak2
  fitgfpar[4] = fitfun -> GetParameter(4); // mean2
  fitgfpar[5] = fitfun -> GetParameter(5); // width2
  
  gauss_sum_resid -> SetParameters(fitgfpar);
  gauss_sum_resid -> SetLineColor(kBlack);
  gauss_sum_resid -> SetLineStyle(1);
  gauss_sum_resid -> SetLineWidth(2);
  gauss_sum_resid -> SetNpx(5000);

  //gauss_sum_resid -> Draw();

  //
  TCanvas *cv_resid = new TCanvas("cv_" + cv_type + "_resid", "cv_" + cv_type + "_resid", 0, 0, 700, 700);
  cv_resid -> SetBottomMargin(0.15);//0.007
  cv_resid -> SetLeftMargin(0.15);

  cv -> SetBottomMargin(0.15);//0.007
  cv -> SetLeftMargin(0.15);

  h1d_resid -> GetYaxis() -> SetTitle("Events fraction");
  h1d_resid -> GetYaxis() -> SetTitleOffset(1.3);
  h1d_resid -> GetYaxis() -> CenterTitle();
  h1d_resid -> GetYaxis() -> SetLabelSize(0.05);
  h1d_resid -> GetYaxis() -> SetTitleSize(0.06);
  h1d_resid -> GetYaxis() -> SetRangeUser(0., 1.2 * ymax_resid);

  h1d_resid -> GetXaxis() -> SetNdivisions(505);
  h1d_resid -> GetXaxis() -> SetTitle("M_{3#pi} [MeV/c^{2}]");
  h1d_resid -> GetXaxis() -> SetRangeUser(xmin_resid, xmax_resid);
  h1d_resid -> GetXaxis() -> SetTitleOffset(1.);
  h1d_resid -> GetXaxis() -> CenterTitle();
  h1d_resid -> GetXaxis() -> SetLabelSize(0.05);
  h1d_resid -> GetXaxis() -> SetTitleSize(0.06);

  h1d_resid -> Draw();
  */
  
  // save
  cv -> SaveAs("IM3pi_" + cv_type + ".pdf");
  
    
  return 0;
  
}
