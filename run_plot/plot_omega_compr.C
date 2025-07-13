#include "../header/plot.h"
#include "../header/graph.h"
#include "../header/plot_omega_compr.h"
//#include "../header/omega_norm.txt"
//#include "../header/omega_vertex.txt"
#include "../header/file.h"

int plot_omega_compr() {

  gErrorIgnoreLevel = kError;
  TGaxis::SetMaxDigits(4);
  //gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetFitFormat("6.4g");
  
  cout << "Comparing norm and vertex" << endl;

  //// norm result
  //TFile *f1 = new TFile("../../plot_omega/plot_omega_norm.root");
  TFile *f1 = new TFile("../../plot_omega/plot_omega_" + file_type1 + ".root");
  getObj(f1);

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
  cout << "OMEGA_PARA_F1[0] = " << OMEGA_PARA_F1[0] << "+/-" << OMEGA_PARA_ERR_F1[0] << "\n"
       << "OMEGA_PARA_F1[1] = " << OMEGA_PARA_F1[1] << "+/-" << OMEGA_PARA_ERR_F1[1] << "\n"
       << "OMEGA_PARA_F1[2] = " << OMEGA_PARA_F1[2] << "+/-" << OMEGA_PARA_ERR_F1[2] << "\n"
       << "Lumi_int_fit_f1 = " << Lumi_int_fit_f1 << "\n"
       << "chi2_sum_f1 = " << chi2_sum_f1 << "\n"
       << "ndf_f1 = " << ndf_f1 << "\n\n";

  //n3pi_obs_f1 -> SetMarkerSize(.8);
  //n3pi_obs_f1 -> GetYaxis() -> SetTitleOffset(1);
  //n3pi_obs_f1 -> GetYaxis() -> SetTitle("#sigma_{3#pi}");
  //n3pi_obs_f1 -> GetYaxis() -> SetTitle(TString::Format("Events/[%0.2f MeV/c^{2}]", Delta_m3pi));
  //n3pi_obs_f1 -> GetXaxis() -> SetTitle("M_{3#pi} [MeV/c^{2}]");
  //n3pi_obs_f1 -> SetName("n3pi_obs_f1");

  //// vertex result
  TFile *f2 = new TFile("../../plot_omega/plot_omega_" + file_type2 + ".root");
  //TFile *f2 = new TFile("../../plot_omega/plot_omega_vertex.root");
  getObj(f2);

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
  cout << "OMEGA_PARA_F2[0] = " << OMEGA_PARA_F2[0] << "+/-" << OMEGA_PARA_ERR_F2[0] << "\n"
       << "OMEGA_PARA_F2[1] = " << OMEGA_PARA_F2[1] << "+/-" << OMEGA_PARA_ERR_F2[1] << "\n"
       << "OMEGA_PARA_F2[2] = " << OMEGA_PARA_F2[2] << "+/-" << OMEGA_PARA_ERR_F2[2] << "\n"
       << "Lumi_int_fit_f2 = " << Lumi_int_fit_f2 << "\n"
       << "chi2_sum_f2 = " << chi2_sum_f2 << "\n"
       << "ndf_f2 = " << ndf_f2 << "\n\n";

  //// Plot
  TCanvas *cv_n3pi_fit = new TCanvas("cv_n3pi_fit", "Nobs 3pi", 0, 0, 1500, 800);

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
  pt1 -> AddText(Form("#Gamma_{#omega}=%0.2f#pm%0.2f MeV", OMEGA_PARA_F1[1], OMEGA_PARA_ERR_F1[1])); 
  pt2 -> AddText(Form("B_{ee}B_{3#pi}=(%0.2f#pm%0.2f)#times10^{-5}", OMEGA_PARA_F1[2] * 1e5, OMEGA_PARA_ERR_F1[2] * 1e5));
  pt3 -> AddText(Form("M_{#omega}=%0.2f#pm%0.2f MeV/c^{2}", OMEGA_PARA_F1[0], OMEGA_PARA_ERR_F1[0])); 
  //pt4 -> AddText(Form("#chi^{2}/ndf=%0.2f", chi2_sum_f1 / ndf_f1)); 
  pt4 -> AddText(Form("#chi^{2}/ndf=%0.1f/%0.2d", chi2_sum_f1, ndf_f1)); 

  pt11 -> AddText(Form("#Gamma_{#omega}=%0.2f#pm%0.2f MeV", OMEGA_PARA_F2[1], OMEGA_PARA_ERR_F2[1])); 
  pt22 -> AddText(Form("B_{ee}B_{3#pi}=(%0.2f#pm%0.2f)#times10^{-5}", OMEGA_PARA_F2[2] * 1e5, OMEGA_PARA_ERR_F2[2] * 1e5));
  pt33 -> AddText(Form("M_{#omega}=%0.2f#pm%0.2f MeV/c^{2}", OMEGA_PARA_F2[0], OMEGA_PARA_ERR_F2[0])); 
  //pt44 -> AddText(Form("#chi^{2}/ndf=%0.2f", chi2_sum_f2 / ndf2));  
  pt44 -> AddText(Form("#chi^{2}/ndf=%0.1f/%0.2d", chi2_sum_f2, ndf_f2)); 
 
  SetGFAttr(n3pi_obs_f2, "M_{3#pi} [MeV/c^{2}]", " ");
  n3pi_obs_f2 -> SetMarkerStyle(21);
  //n3pi_obs_f2 -> SetMarkerColor(kBlue);
  n3pi_obs_f2 -> SetLineColor(kBlack);
  n3pi_obs_f2 -> SetLineWidth(2);

  n3pi_obs_f2 -> GetXaxis() -> SetTitleOffset(1);
  //n3pi_obs_f2 -> GetXaxis() -> SetLabelOffset(1.2);
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
 
  n3pi_obs_f2 -> Draw("PA");
  n3pi_fit_f2 -> Draw("cSame");
  n3pi_obs_f1 -> Draw("pSame");
  n3pi_fit_f1 -> Draw("cSame");
  pt1 -> Draw("Same"), pt11 -> Draw("Same");
  pt2 -> Draw("Same"), pt22 -> Draw("Same");
  pt3 -> Draw("Same"), pt33 -> Draw("Same");
  pt4 -> Draw("Same"), pt44 -> Draw("Same");

  TLegend * legd_cv = new TLegend(0.15, 0.8, 0.9, 0.9);
  
  SetLegend(legd_cv);
  legd_cv -> SetNColumns(4);
  
  legd_cv -> AddEntry(n3pi_obs_f1, "Data " + legend_type1, "lep");
  legd_cv -> AddEntry(n3pi_fit_f1, "Fit " + legend_type1, "l");
  legd_cv -> AddEntry(n3pi_obs_f2, "Data (vertex cut)", "lep");
  legd_cv -> AddEntry(n3pi_fit_f2, "Fit (vertex cut)", "l");
  
  legd_cv -> Draw("Same");

  legtextsize(legd_cv, 0.04);
  
  // save
  cv_n3pi_fit -> SaveAs(outputFolder + "omegaPara_" + file_type1 + "_"  + file_type2 + ".pdf");
  
  return 0;
  
}
