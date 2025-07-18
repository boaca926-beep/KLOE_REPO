#include "../header/sm_para.h"
#include "../header/plot.h"
//#include "../header/hist.h"
#include "../header/efficy.h"
#include "../header/plot_efficy.h"
#include "../header/method.h"
#include "../header/graph.h"

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

int plot_efficy() {

  gErrorIgnoreLevel = kError;
  TGaxis::SetMaxDigits(4);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetFitFormat("6.4g");

  cout << "Plot efficiency ..." << endl;

  cout << "input path: " << input_folder << endl;

  TFile *f_input = new TFile(input_folder + "/efficy.root");
  getObj(f_input);

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
  int nPoints = gf_efficy_sig -> GetN();
  double *x_efficy_sig = gf_efficy_sig -> GetX();
  double *x_efficy_sig_err = gf_efficy_sig -> GetEX();

  double *y_efficy_sig = gf_efficy_sig -> GetY();
  double *y_efficy_sig_err = gf_efficy_sig -> GetEY();
  
  double *y_efficy_ufo = gf_efficy_ufo -> GetY();
  double *y_efficy_ufo_err = gf_efficy_ufo -> GetEY();
  
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
    
    //cout << "point " << i << ", mass =" << x_efficy_sig[i] << ", efficy_sig = " << y_efficy_sig[i] << "+/-" << y_efficy_sig_err[i] << ", efficy_ufo = " << y_efficy_ufo[i] << "+/-" << y_efficy_ufo_err[i] << ", efficy ratio (ufo/sig) = " << RATIO[i] << "+/-" << RATIO_ERR[i] << endl;
    
  }

  TGraphErrors *gf_ratio = get_graph_syst(x_efficy_sig, RATIO, x_efficy_sig_err, RATIO_ERR, nPoints);
  gf_ratio -> SetName("gf_ratio");
  //gf_ratio -> Draw("AP");

  /*
  const double mass_min = 760., mass_max = 800.;
  
  gf_efficy_sig -> GetXaxis() -> SetRangeUser(mass_min, mass_max);
  gf_efficy_sig -> Draw("AP");

  gf_efficy_ufo -> GetXaxis() -> SetRangeUser(mass_min, mass_max);
  gf_efficy_ufo -> Draw("P");

  gf_ratio -> GetXaxis() -> SetRangeUser(mass_min, mass_max);
  gf_ratio -> Draw("P");
  */
  
  TFile *f_output = new TFile(input_folder + "/efficy_ratio.root", "update");
  gf_ratio -> Write();
  
  // plot
  //TCanvas *cv_nb_sig = plotting_nb("cv_nb_sig", "Number of signal events", gf_nb_sel_sig, gf_nb_evtcls_sig, gf_efficy_sig, "Efficiency (#tilde{#varepsilon}_{sig})", 5e4);
  
  //TCanvas *cv_nb_ufo = plotting_nb("cv_nb_ufo", "Number of data events", gf_nb_sel_ufo, gf_nb_evtcls_ufo, gf_efficy_ufo, "Efficiency (#tilde{#varepsilon}_{ufo})", 1000.);
  
  TCanvas *cv_efficy = plotting_efficy("cv_efficy", "Efficiency Comparsion", gf_efficy_sig, gf_efficy_ufo, gf_ratio);
  
  // save
  cv_efficy -> SaveAs(input_folder + "/cv_efficy_" + systType + ".pdf");
  //cv_nb_sig -> SaveAs(input_folder + "/cv_nb_sig_" + systType + ".pdf");
  //cv_nb_ufo -> SaveAs(input_folder + "/cv_nb_ufo_" + systType + ".pdf");
  
  cout << input_folder << endl;

  f_output -> Close();
  
  return 0;
  
}
