#include "../header/sm_para.h"
#include "../header/plot.h"
//#include "../header/hist.h"
#include "../header/efficy.h"
#include "../header/plot_efficy.h"
#include "../header/method.h"
#include "../header/graph.h"

int plot_efficy() {

  gErrorIgnoreLevel = kError;
  TGaxis::SetMaxDigits(4);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetFitFormat("6.4g");

  cout << "Plot efficiency ..." << endl;

  cout << "input path: " << input_folder << endl;

  TFile *f_input = new TFile(input_folder + "/efficy.root");
  //getObj(f_input);

  // EFFICY
  TGraphErrors* gf_efficy_sig = (TGraphErrors*)f_input -> Get("gf_efficy_TISR3PI_SIG");
  gf_efficy_sig -> SetLineColor(kBlue);
  gf_efficy_sig -> SetMarkerColor(kBlue);
  gf_efficy_sig -> SetMarkerSize(.8);
  gf_efficy_sig -> SetMarkerStyle(20);
  
  TGraphErrors* gf_efficy_ufo = (TGraphErrors*)f_input -> Get("gf_efficy_TUFO");
  gf_efficy_ufo -> SetLineColor(kBlack);
  gf_efficy_ufo -> SetMarkerColor(kBlack);
  gf_efficy_ufo -> SetMarkerSize(.8);
  gf_efficy_ufo -> SetMarkerStyle(22);
  //gf_efficy_ufo -> Draw("AP");
  
  // NB UFO
  TGraphErrors* gf_nb_sel_ufo = (TGraphErrors*)f_input -> Get("gf_nb_sel_TUFO");
  gf_nb_sel_ufo -> SetLineColor(kBlack);
  gf_nb_sel_ufo -> SetMarkerColor(kBlack);
  gf_nb_sel_ufo -> SetMarkerSize(.8);
  gf_nb_sel_ufo -> SetMarkerStyle(22);
  //gf_nb_sel_ufo -> Draw("AP");

  TGraphErrors* gf_nb_evtcls_ufo = (TGraphErrors*)f_input -> Get("gf_nb_evtcls_TUFO");
  gf_nb_evtcls_ufo -> SetLineColor(kBlue);
  gf_nb_evtcls_ufo -> SetMarkerColor(kBlue);
  gf_nb_evtcls_ufo -> SetMarkerSize(.8);
  gf_nb_evtcls_ufo -> SetMarkerStyle(20);
  //gf_nb_evtcls_ufo -> Draw("P");
  
  // NB SIG
  TGraphErrors* gf_nb_sel_sig = (TGraphErrors*)f_input -> Get("gf_nb_sel_TISR3PI_SIG");
  gf_nb_sel_sig -> SetLineColor(kBlack);
  gf_nb_sel_sig -> SetMarkerColor(kBlack);
  gf_nb_sel_sig -> SetMarkerSize(.8);
  gf_nb_sel_sig -> SetMarkerStyle(22);
  //gf_nb_sel_sig -> Draw("AP");

  TGraphErrors* gf_nb_evtcls_sig = (TGraphErrors*)f_input -> Get("gf_nb_evtcls_TISR3PI_SIG");
  gf_nb_evtcls_sig -> SetLineColor(kBlue);
  gf_nb_evtcls_sig -> SetMarkerColor(kBlue);
  gf_nb_evtcls_sig -> SetMarkerSize(.8);
  gf_nb_evtcls_sig -> SetMarkerStyle(20);
  //gf_nb_evtcls_sig -> Draw("P");
  
  // calcualte ratio
  //gf_efficy_sig->GetN(); //gf_efficy_sig->GetN();
  //bins = 242;
  const int nPoints = 242; //gf_efficy_sig->GetN();
  
  double *x_efficy_sig = gf_efficy_sig -> GetX();
  double *x_efficy_sig_err = gf_efficy_sig -> GetEX();

  //double x_efficy_sig_GeV[nPoints]; //  = gf_efficy_sig -> GetX();
  //double x_efficy_sig_err_GeV[nPoints]; //= gf_efficy_sig -> GetEX();

  double *y_efficy_sig = gf_efficy_sig -> GetY();
  double *y_efficy_sig_err = gf_efficy_sig -> GetEY();
  
  double *y_efficy_ufo = gf_efficy_ufo -> GetY();
  double *y_efficy_ufo_err = gf_efficy_ufo -> GetEY();

  //std::vector<double> RATIO(nPoints);
  //std::vector<double> RATIO_ERR(nPoints);
  double RATIO[nPoints], RATIO_ERR[nPoints];

  cout << "nPoints = " << nPoints << endl;
  
  for (int i = 0; i < nPoints; i ++) {

    if (y_efficy_sig[i] == 0. || y_efficy_ufo[i] == 0.) {
      //cout << y_efficy_ufo[i] << ", " <<  y_efficy_sig[i] << endl;
      RATIO[i] = 1.;
      RATIO_ERR[i] = 0.;
    }
    else {
      RATIO[i] = get_ratio(y_efficy_ufo[i], y_efficy_sig[i]); //y_efficy_ufo[i] / y_efficy_sig[i];
      RATIO_ERR[i] = ratioErr(y_efficy_ufo[i], y_efficy_ufo_err[i], y_efficy_sig[i], y_efficy_sig_err[i]);
  
    }

    //x_efficy_sig_GeV[i] = x_efficy_sig[i] * 1e-3;
    //x_efficy_sig_err_GeV[i] = x_efficy_sig_err[i];
    
    //cout << x_efficy_sig[i] << "+/-" << x_efficy_sig_err[i] << ", in GeV: " << x_efficy_sig_GeV[i] << "+/-" << x_efficy_sig_err_GeV[i] << endl;
    
    //cout << "point " << i << ", mass = " << x_efficy_sig[i] << ", efficy_sig = " << y_efficy_sig[i] << "+/-" << y_efficy_sig_err[i] << ", efficy_ufo = " << y_efficy_ufo[i] << "+/-" << y_efficy_ufo_err[i] << ", efficy ratio (ufo/sig) = " << RATIO[i] << "+/-" << RATIO_ERR[i] << endl;
    
  }
  

  TGraphErrors *gf_ratio = get_graph_syst(x_efficy_sig, RATIO, x_efficy_sig_err, RATIO_ERR, nPoints);
  gf_ratio -> SetName("gf_ratio");
  //gf_ratio -> Draw("AP");

  //TGraphErrors *gf_ratio_GeV = get_graph_syst(x_efficy_sig_GeV, RATIO, x_efficy_sig_err_GeV, RATIO_ERR, nPoints);
  //gf_ratio_GeV -> SetName("gf_ratio_GeV");
  //gf_ratio_GeV -> Draw();
  
  TGraphErrors* gf_ratio_poly2 = (TGraphErrors*)gf_ratio -> Clone("gf_ratio_poly2");
  TGraphErrors* gf_ratio_linear = (TGraphErrors*)gf_ratio -> Clone("gf_ratio_linear");
  TGraphErrors* gf_ratio_const = (TGraphErrors*)gf_ratio -> Clone("gf_ratio_const");
  TGraphErrors* gf_ratio_cloned = (TGraphErrors*)gf_ratio -> Clone("gf_ratio_cloned");
  TGraphErrors* gf_ratio_cloned1 = (TGraphErrors*)gf_ratio -> Clone("gf_ratio_cloned1");
  
  const double fit_range1 = 758., fit_range2 = 800.;
  //const double fit_range1 = 750., fit_range2 = 810.;
  TVirtualFitter::SetDefaultFitter("Minuit2");

  // fit gf_ratio to a const
  TFitResultPtr fitResult_const = gf_ratio_const -> Fit("pol0", "FS", "", fit_range1, fit_range2);
  TF1 *f_ratio_const = gf_ratio_const -> GetFunction("pol0");
  f_ratio_const -> SetLineWidth(2);
  f_ratio_const -> SetLineColor(1);
  f_ratio_const -> SetNpx(5000);
  
  // fit gf_ratio to a linear
  TFitResultPtr fitResult_linear = gf_ratio_linear -> Fit("pol1", "FS", "", fit_range1, fit_range2);
  TF1 *f_ratio_linear = gf_ratio_linear -> GetFunction("pol1");
  f_ratio_linear -> SetLineWidth(2);
  f_ratio_linear -> SetLineColor(1);
  f_ratio_linear -> SetNpx(5000);
  
  // fit gf_ratio to pol2
  TFitResultPtr fitResult = gf_ratio_poly2 -> Fit("pol2", "FS", "", fit_range1, fit_range2);
  TF1 *f_ratio = gf_ratio_poly2 -> GetFunction("pol2");
  f_ratio -> SetLineWidth(2);
  f_ratio -> SetLineColor(1);
  f_ratio -> SetNpx(5000);

  double p0 = f_ratio -> GetParameter(0);
  double p1 = f_ratio -> GetParameter(1);
  double p2 = f_ratio -> GetParameter(2);

  double p0_err = f_ratio -> GetParError(0);
  double p1_err = f_ratio -> GetParError(1);
  double p2_err = f_ratio -> GetParError(2);

  cout << "p0 = " << p0 << "+/-" << p0_err << "\n"
       << "p1 = " << p1 << "+/-" << p1_err << "\n"
       << "p2 = " << p2 << "+/-" << p2_err << "\n";

  
  
  /*
  // Define the sine function
  TF1 *f_sine = new TF1("f_sine", "[0]*sin([1]*x + [2]) + [3]", fit_range1, fit_range2);
  TF1 *f_cos = new TF1("f_cos", "[0]*cos([1]*x + [2]) + [3]", fit_range1, fit_range2);
  
  f_sine->SetParNames("Amplitude", "Frequency", "Phase", "Offset");

  // Set initial parameter guesses (important for convergence)
  f_sine->SetParameter(0, 1.0);   // Amplitude (A)
  f_sine->SetParameter(1, 0.1);   // Frequency (B)
  f_sine->SetParameter(2, 0.0);   // Phase (C)
  f_sine->SetParameter(3, 0.0);   // Offset (D)
 
  gf_ratio_cloned1 -> Fit("f_sine", "F", "", fit_range1, fit_range2);
  TF1 *f_ratio_linear = gf_ratio_cloned1 -> GetFunction("f_sine");
  */

  /*
  // fit gf_ratio to pol3
  gf_ratio_cloned -> Fit("pol3", "S", "", fit_range1, fit_range2);
  TF1 *f_ratio1 = gf_ratio_cloned -> GetFunction("pol3");
  f_ratio1 -> SetLineWidth(2);
  f_ratio1 -> SetLineColor(kRed);
  f_ratio1 -> SetNpx(5000);

  double p30 = f_ratio1 -> GetParameter(0);
  double p31 = f_ratio1 -> GetParameter(1);
  double p32 = f_ratio1 -> GetParameter(2);
  double p33 = f_ratio1 -> GetParameter(3);

  double p30_err = f_ratio1 -> GetParError(0);
  double p31_err = f_ratio1 -> GetParError(1);
  double p32_err = f_ratio1 -> GetParError(2);
  double p33_err = f_ratio1 -> GetParError(3);

  cout << "p30 = " << p30 << "+/-" << p30_err << "\n"
       << "p31 = " << p31 << "+/-" << p31_err << "\n"
       << "p32 = " << p32 << "+/-" << p32_err << "\n"
       << "p33 = " << p33 << "+/-" << p33_err << "\n";
  */
  
  /*
  TCanvas *cv = new TCanvas("cv_title", "cv", 1200, 800);
  f_ratio -> Draw();
  f_ratio1 -> Draw("same");
  */
  
  //
  //double PARA_FIT_POL3[4] = {p30, p31, p32, p33};
  //double PARA_FIT_POL3_ERR[4] = {p30_err, p31_err, p32_err, p33_err};

  // get corrected ratio
  double *x_gf = gf_ratio_const -> GetX();
  double *y_gf = gf_ratio_const -> GetY();
  double *y_gf_err = gf_ratio_const -> GetEY();

  double x1 = 0., y1 = 0.;
  double x2 = 0., y2 = 0.;

  gf_ratio_poly2 -> GetPoint(0, x1, y1);
  gf_ratio_poly2 -> GetPoint(1, x2, y2);

  double Delta_m3pi = x2 - x1;
  double efficy_ratio = 0.;
  double efficy_ratio_corr = 0., efficy_ratio_corr_err = 0.;
  double efficy_ratio_p3 = 0.;
  double EFFICY_RATIO_CORR[nPoints];
  double EFFICY_RATIO_CORR_ERR[nPoints];
  
  cout << "Delta_m3pi = " << Delta_m3pi << ", nPoints = " << nPoints << endl;
  
  for (int i = 0; i < nPoints; i ++) {

    if (x_gf[i] >= mass_min - Delta_m3pi && x_gf[i] <= mass_max) {

      //efficy_ratio_p3 = f_ratio1 -> Eval(x_gf[i]); // 3nd poly fit results
      efficy_ratio_p3 = 0.;
      
      //efficy_ratio_corr = f_ratio -> Eval(x_gf[i]); // 2nd poly fit correction
      efficy_ratio_corr = f_ratio_const -> Eval(x_gf[i]); // const fit correction
      
      //efficy_ratio_corr_err = get_efficy_ratio_err(fitResult, x_gf[i]);
      
      efficy_ratio_corr_err = 0.;// TMath::Abs(efficy_ratio_corr - efficy_ratio_p3); //get_efficy_ratio_err(PARA_FIT, PARA_FIT_ERR, x_gf[i]);

      //cout << "bin = " << i + 1 << "\t mass = " << x_gf[i] << "\t efficy_ratio = " << y_gf[i] << " +/- " << y_gf_err[i] << "\t efficy_ratio_corr = " << efficy_ratio_corr << " +/- " << efficy_ratio_corr_err << endl;

      cout << "bin = " << i + 1 << "\t mass = " << x_gf[i] << "\t efficy_ratio = " << y_gf[i] << " +/- " << y_gf_err[i] << "\t efficy_ratio_p3 = " << efficy_ratio_p3 << "\t efficy_ratio_corr = " << efficy_ratio_corr << " +/- " << efficy_ratio_corr_err  << "\t efficy_diff = " << efficy_ratio_corr_err << endl;

    }
    else {
      efficy_ratio_corr = y_gf[i];
      efficy_ratio_corr_err = y_gf_err[i];
      
    }

    EFFICY_RATIO_CORR[i] = efficy_ratio_corr;
    EFFICY_RATIO_CORR_ERR[i] = efficy_ratio_corr_err;

    //cout << "bin = " << i + 1 << "\t mass = " << x_gf[i] << "\t efficy_ratio = " << y_gf[i] << "+/-" << y_gf_err[i] << ", \t efficy_ratio_corr = " << EFFICY_RATIO_CORR[i] << "+/-" << EFFICY_RATIO_CORR_ERR[i] << endl;

  }

  TGraphErrors *gf_ratio_corr = get_graph_syst(x_efficy_sig, EFFICY_RATIO_CORR, x_efficy_sig_err, EFFICY_RATIO_CORR_ERR, nPoints);
  gf_ratio_corr -> SetName("gf_ratio_corr");
  gf_ratio_corr -> SetLineColor(kRed);
  gf_ratio_corr -> SetMarkerColor(kRed);

  //gf_ratio_corr -> GetXaxis() -> SetRangeUser(mass_min, mass_max);
  //gf_ratio_corr -> Draw("AP");
  //gf_ratio -> Draw("P");
  
  
  
  //const double mass_min = 760., mass_max = 800.;
  
  //gf_efficy_sig -> GetXaxis() -> SetRangeUser(mass_min, mass_max);
  //gf_efficy_sig -> Draw("AP");

  //gf_efficy_ufo -> GetXaxis() -> SetRangeUser(mass_min, mass_max);
  //gf_efficy_ufo -> Draw("P");

  //gf_ratio -> GetXaxis() -> SetRangeUser(mass_min, mass_max);
  //gf_ratio -> Draw("P");
  
  TFile *f_output = new TFile(input_folder + "/efficy_ratio.root", "update");
  gf_ratio_poly2 -> Write();
  gf_ratio_corr -> Write();
  gf_ratio_const -> Write();
  gf_ratio -> Write();
  

  TTree* TRESULT = new TTree("TRESULT", "recreate");
  TRESULT -> SetAutoSave(0);

  
  TRESULT -> Branch("Br_p0", &p0, "Br_p0/D");
  TRESULT -> Branch("Br_p1", &p1, "Br_p1/D");
  TRESULT -> Branch("Br_p2", &p2, "Br_p2/D");
  
  TRESULT -> Fill();
  
  // plot
  TArrayD ymax_nb_sig = get_gf_max(gf_nb_sel_sig);
  TArrayD ymax_nb_ufo = get_gf_max(gf_nb_sel_ufo);
  
  TArrayD ymax_efficy_sig = get_gf_max(gf_efficy_sig);
  TArrayD ymax_efficy_ufo = get_gf_max(gf_efficy_ufo);
  double ymax_efficy = max(ymax_efficy_sig[0], ymax_efficy_ufo[0]);

  //const TString note = "Prescaled (trigger+FILFO)+KSL";
  //const TString note = "trigger+FILFO+Bkg";
  //const TString note = "trigger+FILFO";
  //const TString note = "trigger";
  //const TString note = "2 tracks+3#gamma";
  const TString note = "";
  
  cout << "ymax_nb_sig = " << ymax_nb_sig[0] << ", ymax_nb_ufo = " << ymax_nb_ufo[0] << "\n"
       << "ymax_efficy_sig = " << ymax_efficy_sig[0] << ", ymax_efficy_ufo = " << ymax_efficy_ufo[0] << ", max_efficy = " << ymax_efficy <<  "\n";  
  
  TCanvas *cv_nb_sig = plotting_nb("cv_nb_sig", "Number of signal events", gf_nb_sel_sig, gf_nb_evtcls_sig, gf_efficy_sig, "Efficiency (#tilde{#varepsilon}_{sig})", ymax_nb_sig[0], note);
  
  TCanvas *cv_nb_ufo = plotting_nb("cv_nb_ufo", "Number of data events", gf_nb_sel_ufo, gf_nb_evtcls_ufo, gf_efficy_ufo, "Efficiency (#tilde{#varepsilon}_{ufo})", ymax_nb_ufo[0], note);
  
  //TCanvas *cv_efficy = plotting_efficy("cv_efficy", "Efficiency Comparsion", gf_efficy_sig, gf_efficy_ufo, gf_ratio_const, gf_ratio_corr, ymax_efficy, note);

  // Weighted average of gf_ratio
  //TGraphErrors *gf_ratio_omega_region = (TGraphErrors *)cv_efficy -> FindObject("gf_ratio");

  //int nPoints_ratio = gf_ratio_omega_region -> GetN();
  //double *x_gf_ratio_omega_region = gf_ratio_omega_region -> GetX();
  //double *x_gf_ratio_omega_region_err = gf_ratio_omega_region -> GetEX();

  //cout << nPoints_ratio << endl;
  
  //TCanvas *cv_tmp = new TCanvas("cv_tmp", "gf_ratio_omega_region", 1200, 800);
  //gf_ratio_omega_region -> Draw("AP");
  
  
  // save
  //cv_efficy -> SaveAs(input_folder + "/cv_efficy_" + systType + ".pdf");
  cv_nb_sig -> SaveAs(input_folder + "/cv_nb_sig_" + systType + ".pdf");
  cv_nb_ufo -> SaveAs(input_folder + "/cv_nb_ufo_" + systType + ".pdf");
  
  cout << input_folder << endl;

  TRESULT -> Write();
  
  f_output -> Close();
  
  return 0;
  
}
