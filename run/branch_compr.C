#include "../header/plot.h"
#include "../header/compr.h"

int branch_compr(){

  gErrorIgnoreLevel = kError;
  TGaxis::SetMaxDigits(4);
  //gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  // fill list
  const int nb_point = 7;
  double BRANCH_BAND[nb_point];
  double BRANCH_ERR_BAND[nb_point];
  double XLIST[nb_point], XLIST_ERR[nb_point];
  
  // PDG
  const double branch_pdg = 6.36;
  const double branch_pdg_err = 0.14;

  for (int i = 0; i < nb_point; i ++) {

    BRANCH_BAND[i] = branch_pdg;
    BRANCH_ERR_BAND[i] = branch_pdg_err;

    XLIST[i] = i;
    XLIST_ERR[i] = 0.;
    
  }

  // KLOE ISR (this work)
  const double branch1 = 5.47389;
  const double branch1_err = 0.06;
  const double branch1_syst_err_plus = 0.11;
  const double branch1_syst_err_minus = 0.08;

  double branch1_exl = 0.; // lower error for x
  double branch1_exh = 0.; // higher error for x
  double branch1_eyl = TMath::Sqrt(branch1_err * branch1_err + branch1_syst_err_minus * branch1_syst_err_minus); // lower error for y
  double branch1_eyh = TMath::Sqrt(branch1_err * branch1_err + branch1_syst_err_plus * branch1_syst_err_plus); // higher error for y

  const double BRANCH_EXP1[1] = {branch1};
  const double BRANCH_EXP1_ERR[1] = {branch1_err};
  const double XLIST_EXP1[1] = {0.1};
  const double XLIST_EXP1_ERR[1] = {0};

  const double BRANCH1_EXL[1] = {branch1_exl};
  const double BRANCH1_EXH[1] = {branch1_exh};
  const double BRANCH1_EYL[1] = {branch1_eyl};
  const double BRANCH1_EYH[1] = {branch1_eyh};
  

  // BABR 2004
  const double branch2 = 6.70;
  const double branch2_err = TMath::Sqrt(0.06 * 0.06 + 0.27 * 0.27);

  const double BRANCH_EXP2[1] = {branch2};
  const double BRANCH_EXP2_ERR[1] = {branch2_err};
  const double XLIST_EXP2[1] = {0.2};
  const double XLIST_EXP2_ERR[1] = {0};

  // BABR 2021 (primary)
  // P_omega=Gamma_{omega}(e+ e- -> omega)Branch(omega -> pi+ pi- pi0)=(0.5698 ± 0.0031 ± 0.0082) keV
  const double Gamma_tot = 8.49; // MeV
  const double Gamma_tot_err = 0.08; // MeV
  const double conv_factor = 1e3 * 1e-6; // 1 keV = 1e3 * 1e-6 MeV
  double P_omega = 0.5698 * conv_factor; // MeV
  double P_omega_err = TMath::Sqrt(0.0031 * 0.0031  + 0.0082 * 0.0082) * conv_factor;
    
  const double branch3 = P_omega / Gamma_tot * 1e5;
  const double branch3_err = branch3 * TMath::Sqrt((P_omega_err / P_omega) * (P_omega_err / P_omega) + (Gamma_tot_err / Gamma_tot) * (Gamma_tot_err / Gamma_tot));

  cout << branch3 << " +/- " << branch3_err << endl;    
  
  const double BRANCH_EXP3[1] = {branch3};
  const double BRANCH_EXP3_ERR[1] = {branch3_err};
  const double XLIST_EXP3[1] = {0.3}; 
  const double XLIST_EXP3_ERR[1] = {0};

  // BESIII (primary)
  const double branch4 = 6.94;
  const double branch4_err = TMath::Sqrt(0.08 * 0.08  + 0.16 * 0.16);

  const double BRANCH_EXP4[1] = {branch4};
  const double BRANCH_EXP4_ERR[1] = {branch4_err};
  const double XLIST_EXP4[1] = {0.4};
  const double XLIST_EXP4_ERR[1] = {0};

  // CMD-2 
  const double branch5 = 6.24;
  const double branch5_err = TMath::Sqrt(0.11 * 0.11 + 0.08 * 0.08);

  const double BRANCH_EXP5[1] = {branch5};
  const double BRANCH_EXP5_ERR[1] = {branch5_err};
  const double XLIST_EXP5[1] = {0.5};
  const double XLIST_EXP5_ERR[1] = {0};

  // RVUE 
  const double branch6 = 6.74;
  const double branch6_err = TMath::Sqrt(0.04 * 0.04 + 0.24 * 0.24);

  const double BRANCH_EXP6[1] = {branch6};
  const double BRANCH_EXP6_ERR[1] = {branch6_err};
  const double XLIST_EXP6[1] = {0.6};
  const double XLIST_EXP6_ERR[1] = {0};

  // ND  e+ e- -> pi+ pi- pi0 
  const double branch7 = 6.37;
  const double branch7_err = 0.35;
  
  const double BRANCH_EXP7[1] = {branch7};
  const double BRANCH_EXP7_ERR[1] = {branch7_err};
  const double XLIST_EXP7[1] = {0.7};
  const double XLIST_EXP7_ERR[1] = {0};

  // label
  const double LABEL_INDX[nb_point] = {XLIST_EXP1[0], XLIST_EXP2[0], XLIST_EXP3[0], XLIST_EXP4[0], XLIST_EXP5[0], XLIST_EXP6[0], XLIST_EXP7[0]};
  const double LABEL_GF[nb_point] = {0., 0., 0., 0., 0., 0., 0.};
  const TString EXP_STR[nb_point] = {"KLOE-2", "BABR04", "BABR20", "BESIII", "CMD-2", "RVUE", "ND"};

  // graphs
  TGraphAsymmErrors *gf_exp1 = new TGraphAsymmErrors(1, XLIST_EXP1, BRANCH_EXP1, BRANCH1_EXL, BRANCH1_EXH, BRANCH1_EYL, BRANCH1_EYH);
  //TGraphErrors *gf_exp1 = new TGraphErrors(1, XLIST_EXP1, BRANCH_EXP1, XLIST_EXP1_ERR, BRANCH_EXP1_ERR);
  gf_exp1 -> SetName("gf_exp1");
  gf_exp1 -> SetMarkerStyle(21);
  gf_exp1 -> SetMarkerSize(1.3);
  gf_exp1 -> SetMarkerColor(kRed);

  TGraphErrors *gf_exp2 = new TGraphErrors(1, XLIST_EXP2, BRANCH_EXP2, XLIST_EXP2_ERR, BRANCH_EXP2_ERR);
  gf_exp2 -> SetName("gf_exp2");
  gf_exp2 -> SetMarkerStyle(8);
  gf_exp2 -> SetMarkerSize(1.3);
  gf_exp2 -> SetMarkerColor(kBlue);

  TGraphErrors *gf_exp3 = new TGraphErrors(1, XLIST_EXP3, BRANCH_EXP3, XLIST_EXP3_ERR, BRANCH_EXP3_ERR);
  gf_exp3 -> SetName("gf_exp3");
  gf_exp3 -> SetMarkerStyle(8);
  gf_exp3 -> SetMarkerSize(1.3);
  gf_exp3 -> SetMarkerColor(kBlue);

  TGraphErrors *gf_exp4 = new TGraphErrors(1, XLIST_EXP4, BRANCH_EXP4, XLIST_EXP4_ERR, BRANCH_EXP4_ERR);
  gf_exp4 -> SetName("gf_exp4");
  gf_exp4 -> SetMarkerStyle(8);
  gf_exp4 -> SetMarkerSize(1.3);
  gf_exp4 -> SetMarkerColor(kBlue);

  TGraphErrors *gf_exp5 = new TGraphErrors(1, XLIST_EXP5, BRANCH_EXP5, XLIST_EXP5_ERR, BRANCH_EXP5_ERR);
  gf_exp5 -> SetName("gf_exp5");
  gf_exp5 -> SetMarkerStyle(22);
  gf_exp5 -> SetMarkerSize(1.3);
  gf_exp5 -> SetMarkerColor(kBlack);

  TGraphErrors *gf_exp6 = new TGraphErrors(1, XLIST_EXP6, BRANCH_EXP6, XLIST_EXP6_ERR, BRANCH_EXP6_ERR);
  gf_exp6 -> SetName("gf_exp6");
  gf_exp6 -> SetMarkerStyle(22);
  gf_exp6 -> SetMarkerSize(1.3);
  gf_exp6 -> SetMarkerColor(kBlack);

  TGraphErrors *gf_exp7 = new TGraphErrors(1, XLIST_EXP7, BRANCH_EXP7, XLIST_EXP7_ERR, BRANCH_EXP7_ERR);
  gf_exp7 -> SetName("gf_exp7");
  gf_exp7 -> SetMarkerStyle(33);
  gf_exp7 -> SetMarkerSize(1.3);
  gf_exp7 -> SetMarkerColor(46);

  TGraphErrors *gf_band = new TGraphErrors(nb_point, XLIST, BRANCH_BAND, XLIST_ERR, BRANCH_ERR_BAND);
  gf_band -> SetName("gf_band");
  gf_band -> SetMarkerStyle(20);
  gf_band -> SetMarkerSize(1.1);
  gf_band -> SetLineColor(kBlack);
  gf_band -> SetLineWidth(4);
  gf_band -> SetFillColor(12);
  gf_band -> SetFillStyle(3002);

  TGraph *gf_label = new TGraph(nb_point, LABEL_INDX, LABEL_GF);
  gf_label -> SetName("gf_label");

  TMultiGraph *mg = new TMultiGraph();
  mg -> SetTitle("Exclusion graphs");
  

  const double band_limit = BRANCH_ERR_BAND[0];

  int bin_indx = 0;
  for (int i = 0; i < nb_point; i ++) {
    bin_indx = gf_label -> GetXaxis() -> FindBin(LABEL_INDX[i]);
    cout << bin_indx << ", EXP_STR = " << EXP_STR[i] << endl;
    gf_label -> GetXaxis() -> SetBinLabel(bin_indx, EXP_STR[i]);
  }

  // plot
  TCanvas *cv = new TCanvas("cv_branch_compr", "cv_branch_compr", 0, 0, 1000, 700);
  //cv -> SetGrid();
  //cv -> SetLeftMargin(0.2);
  cv -> SetBottomMargin(0.17);

  TPaveText *pt1 = new TPaveText(0.11, 0.8, 0.65, 0.82, "NDC");
  PteAttr(pt1);
  pt1 -> SetTextSize(0.045);
  pt1 -> SetTextColor(kRed);
  pt1 -> AddText(Form("B_{ee}B_{3#pi}=%0.5g#pm%0.2g (this work)", branch1, branch1_err)); 

  TPaveText *pt2 = new TPaveText(0.11, 0.72, 0.65, 0.74, "NDC");
  PteAttr(pt2);
  pt2 -> SetTextSize(0.045);
  pt2 -> SetTextColor(kBlack);
  pt2 -> AddText(Form("B_{ee}B_{3#pi}=%0.5g#pm%0.2g (PDG)", branch_pdg, branch_pdg_err)); 

  gf_label -> GetXaxis() -> SetLabelSize(0.06);
  gf_label -> GetXaxis() -> SetLabelOffset(0.01);
  gf_label -> GetYaxis() -> SetNdivisions(512);
  gf_label -> GetYaxis() -> SetRangeUser(BRANCH_BAND[0] - 8. * band_limit, BRANCH_BAND[0] + 8. * band_limit);
  gf_label -> GetYaxis() -> SetTitleOffset(.9);
  gf_label -> GetYaxis() -> SetTitle("B(e^{+}e^{-}#rightarrow#omega)B(#omega#rightarrow#pi^{+}#pi^{-}#pi^{0}) [10^{-5}]");
  gf_label -> GetYaxis() -> CenterTitle();
  gf_label -> GetYaxis() -> SetLabelSize(0.04);
  gf_label -> GetYaxis() -> SetTitleSize(0.05);

  mg -> Add(gf_band);

  gf_label -> Draw("AP");
  mg -> Draw("C3");
  gf_exp1 -> Draw("P");
  gf_exp2 -> Draw("P");
  gf_exp3 -> Draw("P");
  gf_exp4 -> Draw("P");
  gf_exp5 -> Draw("P");
  gf_exp6 -> Draw("P");
  gf_exp7 -> Draw("P");

  //pt1 -> Draw("Same");
  //pt2 -> Draw("Same");

  TLegend * legd_gf = new TLegend(0.65, 0.75, .85, 0.85);

  SetLegend(legd_gf);
  legd_gf -> SetNColumns(1);
  
  legd_gf -> AddEntry(gf_band, "PDG", "lf");

  legd_gf -> Draw("Same");

  legtextsize(legd_gf, 0.07);

  gPad -> Update();
  
  // save
  cv -> SaveAs(outputFile + "/branch_compr.pdf");
  

  return 0;

}
