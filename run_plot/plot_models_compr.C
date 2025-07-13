#include "../header/plot.h"
#include "../header/graph.h"
#include "../header/plot_models_compr.h"
//#include "../header/omega_norm.txt"
//#include "../header/omega_vertex.txt"
#include "../header/file.h"

int plot_models_compr() {

  gErrorIgnoreLevel = kError;
  TGaxis::SetMaxDigits(4);
  //gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetFitFormat("6.4g");
  
  cout << "Comparing " << legend_type1 << " and " << legend_type2 << "\n";

  //// norm result
  //TFile *f1 = new TFile("../../plot_omega/plot_omega_norm.root");
  TFile *f1 = new TFile("../../plot_omega/plot_omega_" + file_type1 + ".root");
  //getObj(f1);

  hdata = (TH1D *)f1->Get("hdata");
  Delta_m3pi = getbinwidth(hdata);

  n3pi_obs_f1 = (TGraphErrors *)f1->Get("gf_n3pi_obs");
  n3pi_fit_f1 = (TGraphErrors *)f1->Get("gf_n3pi_fit");
  f1->GetObject("TRESULT", TCRX3PI_f1);
  
  for (int i = 0; i < TCRX3PI_f1 -> GetEntries(); i++) {
    TCRX3PI_f1 -> GetEntry(i);
    OMEGA_PARA_F1[0] = TCRX3PI_f1 -> GetLeaf("Br_OMEGA_PARA") -> GetValue(0);
    OMEGA_PARA_F1[1] = TCRX3PI_f1 -> GetLeaf("Br_OMEGA_PARA") -> GetValue(1);
    OMEGA_PARA_F1[2] = TCRX3PI_f1 -> GetLeaf("Br_OMEGA_PARA") -> GetValue(2);

    OMEGA_PARA_ERR_F1[0] = TCRX3PI_f1 -> GetLeaf("Br_OMEGA_PARA_ERR") -> GetValue(0);
    OMEGA_PARA_ERR_F1[1] = TCRX3PI_f1 -> GetLeaf("Br_OMEGA_PARA_ERR") -> GetValue(1);
    OMEGA_PARA_ERR_F1[2] = TCRX3PI_f1 -> GetLeaf("Br_OMEGA_PARA_ERR") -> GetValue(2);

    chi2_sum_f1 = TCRX3PI_f1 -> GetLeaf("Br_chi2_sum") -> GetValue(0);
    Lumi_int_fit_f1 = TCRX3PI_f1 -> GetLeaf("Br_Lumi_int_fit") -> GetValue(0);
    ndf_f1 = TCRX3PI_f1 -> GetLeaf("Br_ndf") -> GetValue(0);
    
  }

  //// omega parameters
  cout << f1 -> GetName() << endl;
  cout << "mass = " << OMEGA_PARA_F1[0] << "+/-" << OMEGA_PARA_ERR_F1[0] << ", width = " << OMEGA_PARA_F1[1] << "+/-" << OMEGA_PARA_ERR_F1[1] << ", BB = " << OMEGA_PARA_F1[2] << "+/-" << OMEGA_PARA_ERR_F1[2] << "\n\n";
  //cout << "Lumi_int_fit_f1 = " << Lumi_int_fit_f1 << "\n"
  //<< "chi2_sum_f1 = " << chi2_sum_f1 << "\n"
  //<< "ndf_f1 = " << ndf_f1 << "\n\n";

  //n3pi_obs_f1 -> SetMarkerSize(.8);
  //n3pi_obs_f1 -> GetYaxis() -> SetTitleOffset(1);
  //n3pi_obs_f1 -> GetYaxis() -> SetTitle("#sigma_{3#pi}");
  //n3pi_obs_f1 -> GetYaxis() -> SetTitle(TString::Format("Events/[%0.2f MeV/c^{2}]", Delta_m3pi));
  //n3pi_obs_f1 -> GetXaxis() -> SetTitle("M_{3#pi} [MeV/c^{2}]");
  //n3pi_obs_f1 -> SetName("n3pi_obs_f1");

  //// vertex result
  TFile *f2 = new TFile("../../plot_omega/plot_omega_" + file_type2 + ".root");
  //TFile *f2 = new TFile("../../plot_omega/plot_omega_vertex.root");
  //getObj(f2);

  n3pi_obs_f2  = (TGraphErrors *)f2->Get("gf_n3pi_obs");
  n3pi_fit_f2 = (TGraphErrors *)f2->Get("gf_n3pi_fit");
  f2 -> GetObject("TRESULT", TCRX3PI_f2);

  for (int i = 0; i < TCRX3PI_f2 -> GetEntries(); i++) {
    TCRX3PI_f2 -> GetEntry(i);
    OMEGA_PARA_F2[0] = TCRX3PI_f2 -> GetLeaf("Br_OMEGA_PARA") -> GetValue(0);
    OMEGA_PARA_F2[1] = TCRX3PI_f2 -> GetLeaf("Br_OMEGA_PARA") -> GetValue(1);
    OMEGA_PARA_F2[2] = TCRX3PI_f2 -> GetLeaf("Br_OMEGA_PARA") -> GetValue(2);

    OMEGA_PARA_ERR_F2[0] = TCRX3PI_f2 -> GetLeaf("Br_OMEGA_PARA_ERR") -> GetValue(0);
    OMEGA_PARA_ERR_F2[1] = TCRX3PI_f2 -> GetLeaf("Br_OMEGA_PARA_ERR") -> GetValue(1);
    OMEGA_PARA_ERR_F2[2] = TCRX3PI_f2 -> GetLeaf("Br_OMEGA_PARA_ERR") -> GetValue(2);

    chi2_sum_f2 = TCRX3PI_f2 -> GetLeaf("Br_chi2_sum") -> GetValue(0);
    Lumi_int_fit_f2 = TCRX3PI_f2 -> GetLeaf("Br_Lumi_int_fit") -> GetValue(0);
    ndf_f2 = TCRX3PI_f2 -> GetLeaf("Br_ndf") -> GetValue(0);
    
  }

  //// omega parameters
  cout << f2 -> GetName() << endl;
  cout << PARA_LABEL[0] << ": " << OMEGA_PARA_F2[0] << "+/-" << PARA_LABEL[1] << "; " << OMEGA_PARA_ERR_F2[0] << ": " << OMEGA_PARA_F2[1] << "+/-" << OMEGA_PARA_ERR_F2[1] << "; " << PARA_LABEL[2] << "; " << OMEGA_PARA_F2[2] << "+/-" << OMEGA_PARA_ERR_F2[2] << "\n\n";

  /*
  cout << "Lumi_int_fit_f2 = " << Lumi_int_fit_f2 << "\n"
       << "chi2_sum_f2 = " << chi2_sum_f2 << "\n"
       << "ndf_f2 = " << ndf_f2 << "\n\n";
  */

  //// Calculate errors
  double y_diff = 0.;
  double y_uncorr_err = 0.;
  double R = 0.;
  
  for (int i = 0; i < 3; i ++) {// loop over omega parameters

    y_diff = TMath::Abs(OMEGA_PARA_F2[i] - OMEGA_PARA_F1[i]);

    y_uncorr_err = TMath::Sqrt(TMath::Abs(OMEGA_PARA_ERR_F2[i] * OMEGA_PARA_ERR_F2[i] - OMEGA_PARA_ERR_F1[i] * OMEGA_PARA_ERR_F1[i]));

    R = y_diff / y_uncorr_err;
    
    cout << PARA_LABEL[i] << ", |y_diff| = " << y_diff << ", y_uncorr_err = " << y_uncorr_err << ", significance = " << R << endl;

  }
  //// Pull; difference between BW & VMD result
  const int nbins = n3pi_fit_f1 -> GetN();
  //cout << "nbins = " << nbins << endl;

  double ZVALUE[nbins], ZVALUE_ERR[nbins];
  double MVALUE[nbins], MVALUE_ERR[nbins];

  double m3pi = 0.;
  double nb_f1 = 0., nb_err_f1 = 0.;
  double nb_f2 = 0., nb_err_f2 = 0.;

  double Z = 0.;
  double err_diff = 0.;

  for (Int_t  k = 0; k < nbins; k++) {

    n3pi_fit_f1 -> GetPoint(k, m3pi, nb_f1);
    nb_err_f1 = n3pi_fit_f1 -> GetErrorY(k);

    n3pi_fit_f2 -> GetPoint(k, m3pi, nb_f2);
    nb_err_f2 = n3pi_fit_f2 -> GetErrorY(k);

    // calculate significance Z (sigma)
    err_diff = TMath::Sqrt(nb_err_f1 * nb_err_f1 + nb_err_f2 * nb_err_f2);
    Z = TMath::Abs(nb_f1 - nb_f2) / err_diff;

    ZVALUE[k] = Z;
    ZVALUE_ERR[k] = 0.;

    MVALUE[k] = m3pi;
    MVALUE_ERR[k] = 0.;
    
    //cout << k << ": nb_f1 = " << nb_f1 << "+/-" << nb_err_f1 << ": nb_f2 = " << nb_f2 << "+/-" << nb_err_f2 << "; significance Z = " << Z << " [sigma], err_diff = " << err_diff << "\n";
    
  }

  TGraphErrors *gf_Z = new TGraphErrors(nbins, MVALUE, ZVALUE, MVALUE_ERR, ZVALUE_ERR);
  gf_Z -> Draw();

  //// Plot
  TCanvas *cv_n3pi_fit = new TCanvas("cv_n3pi_fit", "Nobs 3pi", 0, 0, 1200, 800);

  cv_n3pi_fit -> SetLeftMargin(0.1);
  cv_n3pi_fit -> SetBottomMargin(0.1);//0.007
  
  char display1[50], display2[50], display3[50], display4[50], display5[50];

  TPaveText *pt1 = new TPaveText(0.12, 0.70, 0.43, 0.72, "NDC");
  TPaveText *pt2 = new TPaveText(0.12, 0.64, 0.43, 0.66, "NDC");
  TPaveText *pt3 = new TPaveText(0.12, 0.58, 0.43, 0.60, "NDC");
  TPaveText *pt4 = new TPaveText(0.12, 0.52, 0.43, 0.54, "NDC");

  TPaveText *pt11 = new TPaveText(0.65, 0.70, 0.80, 0.72, "NDC");
  TPaveText *pt22 = new TPaveText(0.65, 0.64, 0.80, 0.66, "NDC");
  TPaveText *pt33 = new TPaveText(0.65, 0.58, 0.80, 0.60, "NDC");
  TPaveText *pt44 = new TPaveText(0.65, 0.52, 0.80, 0.54, "NDC");
  
  PteAttr(pt1); pt1 -> SetTextSize(0.04); pt1 -> SetTextColor(kBlack);
  PteAttr(pt2); pt2 -> SetTextSize(0.04); pt2 -> SetTextColor(kBlack);
  PteAttr(pt3); pt3 -> SetTextSize(0.04); pt3 -> SetTextColor(kBlack);
  PteAttr(pt4); pt4 -> SetTextSize(0.04); pt4 -> SetTextColor(kBlack);

  PteAttr(pt11); pt11 -> SetTextSize(0.04); pt11 -> SetTextColor(kBlue);
  PteAttr(pt22); pt22 -> SetTextSize(0.04); pt22 -> SetTextColor(kBlue);
  PteAttr(pt33); pt33 -> SetTextSize(0.04); pt33 -> SetTextColor(kBlue);
  PteAttr(pt44); pt44 -> SetTextSize(0.04); pt44 -> SetTextColor(kBlue);
  
  //sprintf(display1,"#chi^{2}/ndf=%0.2f", chi2_sum / ndf);
  pt1 -> AddText(Form("M_{#omega}=%0.2f#pm%0.2f MeV/c^{2}", OMEGA_PARA_F2[0], OMEGA_PARA_ERR_F2[0]));
  pt2 -> AddText(Form("#Gamma_{#omega}=%0.2f#pm%0.2f MeV", OMEGA_PARA_F2[1], OMEGA_PARA_ERR_F2[1]));
  pt3 -> AddText(Form("B_{ee}B_{3#pi}=(%0.2f#pm%0.2f)#times10^{-5}", OMEGA_PARA_F2[2] * 1e5, OMEGA_PARA_ERR_F2[2] * 1e5));
  //pt44 -> AddText(Form("#chi^{2}/ndf=%0.2f", chi2_sum_f2 / ndf2));  
  pt4 -> AddText(Form("#chi^{2}/ndf=%0.1f/%0.2d", chi2_sum_f2, ndf_f2)); 
 
  pt11 -> AddText(Form("M_{#omega}=%0.2f#pm%0.2f MeV/c^{2}", OMEGA_PARA_F1[0], OMEGA_PARA_ERR_F1[0]));
  pt22 -> AddText(Form("#Gamma_{#omega}=%0.2f#pm%0.2f MeV", OMEGA_PARA_F1[1], OMEGA_PARA_ERR_F1[1]));
  pt33 -> AddText(Form("B_{ee}B_{3#pi}=(%0.2f#pm%0.2f)#times10^{-5}", OMEGA_PARA_F1[2] * 1e5, OMEGA_PARA_ERR_F1[2] * 1e5));
  //pt4 -> AddText(Form("#chi^{2}/ndf=%0.2f", chi2_sum_f1 / ndf_f1)); 
  pt44 -> AddText(Form("#chi^{2}/ndf=%0.1f/%0.2d", chi2_sum_f1, ndf_f1)); 

  SetGFAttr(n3pi_obs_f2, "M_{3#pi} [MeV/c^{2}]", " ");
  n3pi_obs_f2 -> SetMarkerStyle(21);
  //n3pi_obs_f2 -> SetMarkerColor(kBlue);
  n3pi_obs_f2 -> SetLineColor(kBlack);
  n3pi_obs_f2 -> SetLineWidth(2);

  n3pi_obs_f2 -> GetXaxis() -> SetTitleOffset(1);
  n3pi_obs_f2 -> GetXaxis() -> SetLabelOffset(1.2);
  n3pi_obs_f2 -> GetXaxis() -> SetLabelSize(0.04);
  
  n3pi_obs_f2 -> GetYaxis() -> SetTitle(TString::Format("Events/[%0.2f MeV/c^{2}]", Delta_m3pi));
  n3pi_obs_f2 -> GetYaxis() -> SetTitleSize(0.06);

  //
  SetGFAttr(n3pi_obs_f1, "M_{3#pi} [MeV/c^{2}]", " ");
  n3pi_obs_f1 -> SetMarkerStyle(22);
  n3pi_obs_f1 -> SetMarkerColor(kBlue);
  n3pi_obs_f1 -> SetLineColor(kBlue);
  n3pi_obs_f1 -> SetLineWidth(2);

  SetGFAttr(n3pi_fit_f1, "M_{3#pi} [MeV/c^{2}]", " ");
  n3pi_fit_f1 -> SetMarkerStyle(22);
  n3pi_fit_f1 -> SetMarkerColor(kBlue);
  n3pi_fit_f1 -> SetLineColor(kBlue);
  n3pi_fit_f1 -> SetLineWidth(2);
 
  //n3pi_obs_f1 -> GetXaxis() -> SetTitleOffset(1.2);
  //n3pi_obs_f1 -> GetXaxis() -> SetLabelOffset(1.2);
  //n3pi_obs_f1 -> GetXaxis() -> SetLabelSize(0.03);
  //n3pi_obs_f1 -> GetYaxis() -> SetTitle(TString::Format("Events/[%0.2f MeV/c^{2}]", Delta_m3pi));
  //n3pi_obs_f1 -> GetYaxis() -> SetTitleSize(0.06);

  TPad *p2 = new TPad("p2", "p2", 0., 0., 1., 0.25);
  p2 -> Draw();
  p2 -> SetBottomMargin(0.4);
  p2 -> SetLeftMargin(0.1);
  p2 -> SetGrid();

  TPad *p1 = new TPad("p1", "p1", 0., 0.25, 1., 1.);
  p1 -> Draw();
  p1 -> SetBottomMargin(0.02);//0.007
  p1 -> SetLeftMargin(0.1);
  p1 -> cd();
  
  n3pi_fit_f1 -> SetLineStyle(7);
  
  n3pi_obs_f2 -> Draw("PA");
  n3pi_fit_f2 -> Draw("cSame");
  n3pi_obs_f1 -> Draw("pSame");
  n3pi_fit_f1 -> Draw("cSame");
  pt1 -> Draw("Same"), pt11 -> Draw("Same");
  pt2 -> Draw("Same"), pt22 -> Draw("Same");
  pt3 -> Draw("Same"), pt33 -> Draw("Same");
  pt4 -> Draw("Same"), pt44 -> Draw("Same");

  TLegend * legd_cv = new TLegend(0.15, 0.8, 0.6, 0.9);
  
  SetLegend(legd_cv);
  legd_cv -> SetNColumns(3);

  legd_cv -> AddEntry(n3pi_obs_f2, "Data " + legend_type2, "lep");
  legd_cv -> AddEntry(n3pi_fit_f2, legend_type2, "l");
  //legd_cv -> AddEntry(n3pi_obs_f1, "Data " + legend_type1, "lep");
  legd_cv -> AddEntry(n3pi_fit_f1, legend_type1, "l");
  
  legd_cv -> Draw("Same");

  legtextsize(legd_cv, 0.04);
  
  p2 -> cd();

  gf_Z -> SetLineColor(0);
  gf_Z -> GetYaxis() -> SetNdivisions(505);
  gf_Z -> GetYaxis() -> SetTitleSize(33);
  gf_Z -> GetYaxis() -> SetTitleFont(43);
  gf_Z -> GetYaxis() -> SetTitleOffset(1.2);
  gf_Z -> GetYaxis() -> SetLabelFont(43); // Absolute front size in pixel (precision 3)
  gf_Z -> GetYaxis() -> SetLabelSize(25);
  gf_Z -> GetYaxis() -> CenterTitle();
  gf_Z -> GetYaxis() -> SetTitle("Z [#sigma]");
  gf_Z -> GetYaxis() -> SetRangeUser(-1., 5.);
  //gf_Z -> GetYaxis() -> SetTitle("(N_{d} - N_{mc})/#sqrt{N_{d}}");

  //gf_Z -> GetXaxis() -> SetLabelFont(43); // Absolute font size in pixel (precision 3)
  gf_Z -> GetXaxis() -> SetTitleOffset(0.85);
  //gf_Z -> GetXaxis() -> SetTitleSize(45);
  gf_Z -> GetXaxis() -> SetLabelSize(0.15);
  //gf_Z -> GetXaxis() -> SetLabelOffset(0.2);
  gf_Z -> GetXaxis() -> SetTitleSize(0.2);
  gf_Z -> GetXaxis() -> SetTitle("M_{3#pi} [MeV/c^{2}]");
  gf_Z -> GetXaxis() -> CenterTitle();
  gf_Z -> SetMarkerStyle(21);
  gf_Z -> SetMarkerSize(0.7);
  //gf_Z -> Draw("E0");

  //TMultiGraph *mg = new TMultiGraph();
  //mg -> SetTitle("Exclusion graphs");
  
  //mg -> Add(gf_n3pi_diff);
  gf_Z -> SetFillColor(1);
  gf_Z -> SetLineStyle(9);
  gf_Z -> SetFillStyle(3001);

  //gf_n3pi_vmd_diff -> SetFillColor(4);
  //gf_n3pi_vmd_diff -> SetLineStyle(9);
  //gf_n3pi_vmd_diff -> SetFillStyle(3001);
  
  gf_Z -> Draw();
  //gf_n3pi_vmd_diff -> Draw("3");

  //gPad -> SetLogy();

  // save
  cv_n3pi_fit -> SaveAs(outputFolder + legend_type1 + "_"  + legend_type2 + ".pdf");
  
  return 0;
  
}
