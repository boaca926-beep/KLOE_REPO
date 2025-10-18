#include "../header/plot.h"
#include "../header/compr.h"

int mass_compr(){

  gErrorIgnoreLevel = kError;
  TGaxis::SetMaxDigits(4);
  //gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  // fill list
  const int nb_point = 12;
  double MASS_BAND[nb_point];
  double MASS_ERR_BAND[nb_point];
  double XLIST[nb_point], XLIST_ERR[nb_point];
  
  // PDG
  const double mass_pdg = 782.66;
  const double mass_pdg_err = 0.13;

  for (int i = 0; i < nb_point; i ++) {

    MASS_BAND[i] = mass_pdg;
    MASS_ERR_BAND[i] = mass_pdg_err;

    XLIST[i] = i;
    XLIST_ERR[i] = 0.;
    
  }

  
  // H1 e  p -> e pi+ pi- p
  const double mass1 = 777.9;
  const double mass1_syst_err_averg = (4.3 + 2.2) / 2.;
  const double mass1_err = TMath::Sqrt(2.2 * 2.2 + mass1_syst_err_averg * mass1_syst_err_averg);

  //double MASS_EXP1[1], MASS_EXP1_ERR[1];
  //double XLIST_EXP1[1], XLIST_EXP1_ERR[1];
  
  const double MASS_EXP1[1] = {mass1};
  const double MASS_EXP1_ERR[1] = {mass1_err};
  const double XLIST_EXP1[1] = {0.7};
  const double XLIST_EXP1_ERR[1] = {0};

  // KLOE ISR (this work)
  const double mass2 = 782.73;
  const double mass2_err = 0.0394615;
  const double mass2_syst_err_plus = 0.05;
  const double mass2_syst_err_minus = 0.07;

  double mass2_exl = 0.; // lower error for x
  double mass2_exh = 0.; // higher error for x
  double mass2_eyl = TMath::Sqrt(mass2_err * mass2_err + mass2_syst_err_minus * mass2_syst_err_minus); // lower error for y
  double mass2_eyh = TMath::Sqrt(mass2_err * mass2_err + mass2_syst_err_plus * mass2_syst_err_plus); // higher error for y

  cout << "mass2_eyl = " << mass2_eyl << ", mass2_eyh = " << mass2_eyh << endl; 
    
  const double MASS_EXP2[1] = {mass2};
  const double MASS_EXP2_ERR[1] = {mass2_err};
  const double XLIST_EXP2[1] = {0.1};
  const double XLIST_EXP2_ERR[1] = {0};

  const double MASS2_EXL[1] = {mass2_exl};
  const double MASS2_EXH[1] = {mass2_exh};
  const double MASS2_EYL[1] = {mass2_eyl};
  const double MASS2_EYH[1] = {mass2_eyh};

  // BABAR ISR
  const double mass3 = 782.59 - 0.2; //pdg: Phys. Lett.B 592, 1 (2004)
  const double mass3_err = TMath::Sqrt(0.11 * 0.11 + 0.1 * 0.1);
  //cout << mass3_err << endl;

  const double MASS_EXP3[1] = {mass3};
  const double MASS_EXP3_ERR[1] = {mass3_err};
  const double XLIST_EXP3[1] = {0.2};
  const double XLIST_EXP3_ERR[1] = {0};

  // BABAR ISR (primary)
  const double mass4 = 782.65 - 0.04; //pdg: Phys. 2020, 083C01 (2020) 
  const double mass4_err = TMath::Sqrt(0.12 * 0.12 + 0.06 * 0.06);

  const double MASS_EXP4[1] = {mass4};
  const double MASS_EXP4_ERR[1] = {mass4_err};
  const double XLIST_EXP4[1] = {0.3};
  const double XLIST_EXP4_ERR[1] = {0};

  // BESIII ISR (primary)
  const double mass5 = 783.20; //arXiv:1912.11208 
  const double mass5_err = TMath::Sqrt(0.07 * 0.07 + 0.24 * 0.24);

  const double MASS_EXP5[1] = {mass5};
  const double MASS_EXP5_ERR[1] = {mass5_err};
  const double XLIST_EXP5[1] = {0.4};
  const double XLIST_EXP5_ERR[1] = {0};

  // CMD-2 2004
  const double mass6 = 783.20; // arXiv:hep-ex/0409030 e+ e- -> pi0 gamma
  const double mass6_err = TMath::Sqrt(0.13 * 0.13 + 0.16 * 0.16);

  const double MASS_EXP6[1] = {mass6};
  const double MASS_EXP6_ERR[1] = {mass6_err};
  const double XLIST_EXP6[1] = {0.5};
  const double XLIST_EXP6_ERR[1] = {0};

  // CMD-2 2003
  const double mass7 = 782.68; // arXiv:hep-ex/0308008 e+ e- -> pi+ pi- pi0 
  const double mass7_err = TMath::Sqrt(0.09 * 0.09 + 0.04 * 0.04);

  const double MASS_EXP7[1] = {mass7};
  const double MASS_EXP7_ERR[1] = {mass7_err};
  const double XLIST_EXP7[1] = {0.6};
  const double XLIST_EXP7_ERR[1] = {0};

  // RVUE
  const double mass8 = 782.79; // arXiv:hep-ex/0305049 e+ e- -> pi+ pi- pi0  
  const double mass8_err = TMath::Sqrt(0.09 * 0.09 + 0.08 * 0.08);

  const double MASS_EXP8[1] = {mass8};
  const double MASS_EXP8_ERR[1] = {mass8_err};
  const double XLIST_EXP8[1] = {0.7};
  const double XLIST_EXP8_ERR[1] = {0};

  // CBAR 11k 
  const double mass9 = 781.96; // Phys. Lett.B 327, 425 (1994) p pbar -> omega eta pi0, omega -> pi0 gamma  
  const double mass9_err = TMath::Sqrt(0.17 * 0.17 + 0.80 * 0.80);

  const double MASS_EXP9[1] = {mass9};
  const double MASS_EXP9_ERR[1] = {mass9_err};
  const double XLIST_EXP9[1] = {0.8};
  const double XLIST_EXP9_ERR[1] = {0};

  // CBAR
  const double mass10 = 782.08; // Phys. Lett.B 327, 425 (1994) p pbar -> omega eta pi0, omega -> pi0 gamma  
  const double mass10_err = TMath::Sqrt(0.36 * 0.36 + 0.82 * 0.82);

  const double MASS_EXP10[1] = {mass10};
  const double MASS_EXP10_ERR[1] = {mass10_err};
  const double XLIST_EXP10[1] = {0.9};
  const double XLIST_EXP10_ERR[1] = {0};

  // CBAR  
  const double mass11 = 781.96; // Phys. Lett.B 311, 362 (1993) p pbar -> omega pi0 pi0  
  const double mass11_err = TMath::Sqrt(0.13 * 0.13 + 0.17 * 0.17);

  const double MASS_EXP11[1] = {mass11};
  const double MASS_EXP11_ERR[1] = {mass11_err};
  const double XLIST_EXP11[1] = {1.0};
  const double XLIST_EXP11_ERR[1] = {0};

  // ASTE 
  const double mass12 = 782.4; // Z.Phys.C 59 (1993) 387-398 p pbar -> 2pi+ 2pi- pi0  
  const double mass12_err = 0.2;
  
  const double MASS_EXP12[1] = {mass12};
  const double MASS_EXP12_ERR[1] = {mass12_err};
  const double XLIST_EXP12[1] = {1.1};
  const double XLIST_EXP12_ERR[1] = {0};

  // KEYNE 
  const double mass13 = 782.4; //   
  const double mass13_err = 0.5;
  
  const double MASS_EXP13[1] = {mass13};
  const double MASS_EXP13_ERR[1] = {mass13_err};
  const double XLIST_EXP13[1] = {1.2};
  const double XLIST_EXP13_ERR[1] = {0};
  
  // label
  const double LABEL_INDX[nb_point] = {XLIST_EXP2[0], XLIST_EXP3[0], XLIST_EXP4[0], XLIST_EXP5[0], XLIST_EXP6[0], XLIST_EXP7[0], XLIST_EXP8[0], XLIST_EXP9[0], XLIST_EXP10[0], XLIST_EXP11[0], XLIST_EXP12[0], XLIST_EXP13[0]};
  const double LABEL_GF[nb_point] = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};
  const TString EXP1_STR[nb_point] = {"", "", "", "", "", "", "", "", "", "", "", ""};
  const TString EXP_STR[nb_point] = {"This Work", "BABR04", "BABR20", "BESIII", "CMD-2a", "CMD-2b", "RVUE", "CBARa", "CBARb", "CBARc", "ASTE", "CNTR"};

  const double LABEL_INDX1[nb_point] = {XLIST_EXP1[0]};
  const double LABEL_GF1[1] = {0.};
  
  // graphs
  TGraph *gf_label = new TGraph(nb_point, LABEL_INDX, LABEL_GF);
  gf_label -> SetName("gf_label");

  TGraph *gf_label1 =  (TGraphErrors*) gf_label -> Clone(); 
  gf_label1 -> SetName("gf_label1");

  TGraph *gf_label2 = new TGraph(1, LABEL_INDX1, LABEL_GF1);
  gf_label2 -> SetName("gf_label2");

  
  TGraphErrors *gf_exp1 = new TGraphErrors(1, XLIST_EXP1, MASS_EXP1, XLIST_EXP1_ERR, MASS_EXP1_ERR);
  gf_exp1 -> SetName("gf_exp1");
  gf_exp1 -> SetMarkerStyle(21);
  gf_exp1 -> SetMarkerSize(2.3);
  gf_exp1 -> SetMarkerColor(46);

  TGraphAsymmErrors *gf_exp2 = new TGraphAsymmErrors(1, XLIST_EXP2, MASS_EXP2, MASS2_EXL, MASS2_EXH, MASS2_EYL, MASS2_EYH);
  //TGraphErrors *gf_exp2 = new TGraphErrors(1, XLIST_EXP2, MASS_EXP2, XLIST_EXP2_ERR, MASS_EXP2_ERR);
  gf_exp2 -> SetName("gf_exp2"); // KLOE-2
  gf_exp2 -> SetMarkerStyle(21);
  gf_exp2 -> SetMarkerSize(1.3);
  gf_exp2 -> SetMarkerColor(kRed);

  TGraphErrors *gf_exp3 = new TGraphErrors(1, XLIST_EXP3, MASS_EXP3, XLIST_EXP3_ERR, MASS_EXP3_ERR);
  gf_exp3 -> SetName("gf_exp3");
  gf_exp3 -> SetMarkerStyle(8); //22
  gf_exp3 -> SetMarkerSize(1.3);
  gf_exp3 -> SetMarkerColor(kBlue);

  TGraphErrors *gf_exp4 = new TGraphErrors(1, XLIST_EXP4, MASS_EXP4, XLIST_EXP4_ERR, MASS_EXP4_ERR);
  gf_exp4 -> SetName("gf_exp4");
  gf_exp4 -> SetMarkerStyle(8); //23
  gf_exp4 -> SetMarkerSize(1.3);
  gf_exp4 -> SetMarkerColor(kBlue);

  TGraphErrors *gf_exp5 = new TGraphErrors(1, XLIST_EXP5, MASS_EXP5, XLIST_EXP5_ERR, MASS_EXP5_ERR);
  gf_exp5 -> SetName("gf_exp5");
  gf_exp5 -> SetMarkerStyle(8);
  gf_exp5 -> SetMarkerSize(1.3);
  gf_exp5 -> SetMarkerColor(kBlue);

  TGraphErrors *gf_exp6 = new TGraphErrors(1, XLIST_EXP6, MASS_EXP6, XLIST_EXP6_ERR, MASS_EXP6_ERR);
  gf_exp6 -> SetName("gf_exp6");
  gf_exp6 -> SetMarkerStyle(22);//21
  gf_exp6 -> SetMarkerSize(1.3);
  gf_exp6 -> SetMarkerColor(kBlack);

  TGraphErrors *gf_exp7 = new TGraphErrors(1, XLIST_EXP7, MASS_EXP7, XLIST_EXP7_ERR, MASS_EXP7_ERR);
  gf_exp7 -> SetName("gf_exp7");
  gf_exp7 -> SetMarkerStyle(22);
  gf_exp7 -> SetMarkerSize(1.3);
  gf_exp7 -> SetMarkerColor(kBlack);

  TGraphErrors *gf_exp8 = new TGraphErrors(1, XLIST_EXP8, MASS_EXP8, XLIST_EXP8_ERR, MASS_EXP8_ERR);
  gf_exp8 -> SetName("gf_exp8");
  gf_exp8 -> SetMarkerStyle(22);//23
  gf_exp8 -> SetMarkerSize(1.3);
  gf_exp8 -> SetMarkerColor(kBlack);

  TGraphErrors *gf_exp9 = new TGraphErrors(1, XLIST_EXP9, MASS_EXP9, XLIST_EXP9_ERR, MASS_EXP9_ERR);
  gf_exp9 -> SetName("gf_exp9");
  gf_exp9 -> SetMarkerStyle(23);
  gf_exp9 -> SetMarkerSize(1.3);
  gf_exp9 -> SetMarkerColor(28);

  TGraphErrors *gf_exp10 = new TGraphErrors(1, XLIST_EXP10, MASS_EXP10, XLIST_EXP10_ERR, MASS_EXP10_ERR);
  gf_exp10 -> SetName("gf_exp10");
  gf_exp10 -> SetMarkerStyle(23);
  gf_exp10 -> SetMarkerSize(1.3);
  gf_exp10 -> SetMarkerColor(28);

  TGraphErrors *gf_exp11 = new TGraphErrors(1, XLIST_EXP11, MASS_EXP11, XLIST_EXP11_ERR, MASS_EXP11_ERR);
  gf_exp11 -> SetName("gf_exp11");
  gf_exp11 -> SetMarkerStyle(23);
  gf_exp11 -> SetMarkerSize(1.3);
  gf_exp11 -> SetMarkerColor(28);

  TGraphErrors *gf_exp12 = new TGraphErrors(1, XLIST_EXP12, MASS_EXP12, XLIST_EXP12_ERR, MASS_EXP12_ERR);
  gf_exp12 -> SetName("gf_exp12");
  gf_exp12 -> SetMarkerStyle(23);
  gf_exp12 -> SetMarkerSize(1.3);
  gf_exp12 -> SetMarkerColor(28);

  TGraphErrors *gf_exp13 = new TGraphErrors(1, XLIST_EXP13, MASS_EXP13, XLIST_EXP13_ERR, MASS_EXP13_ERR);
  gf_exp13 -> SetName("gf_exp13");
  gf_exp13 -> SetMarkerStyle(33);
  gf_exp13 -> SetMarkerSize(1.3);
  gf_exp13 -> SetMarkerColor(46);

  TGraphErrors *gf_band = new TGraphErrors(nb_point, XLIST, MASS_BAND, XLIST_ERR, MASS_ERR_BAND);
  gf_band -> SetName("gf_band");
  //gf_band -> SetMarkerStyle(33);
  gf_band -> SetMarkerSize(1.1);
  gf_band -> SetLineColor(kBlack);
  gf_band -> SetLineWidth(2);
  gf_band -> SetFillColor(12);
  gf_band -> SetFillStyle(3002);

  TMultiGraph *mg = new TMultiGraph();
  mg -> SetTitle("Exclusion graphs");
  
  cout << "PDG: " << mass_pdg << " +/- " << mass_pdg_err << "\n"
       << "H1: " << mass1 << " +/- " << mass1_err << "\n"
       << "H1-PDG: " << TMath::Abs(mass_pdg - mass1) / mass1_err << " sigma\n\n"
       << "1. this work: " << mass2 << " +/- " << mass2_err << "\n"
       << "2. BABAR (2004): " << mass3 << " +/- " << mass3_err << "\n"
       << "3. BABAR (2021): " << mass4 << " +/- " << mass4_err << "\n"
       << "4. BESIII (2019): " << mass5 << " +/- " << mass5_err << "\n"
       << "5. CMD-2 (2004): " << mass6 << " +/- " << mass6_err << "\n"
       << "6. CMD-2 (2003): " << mass7 << " +/- " << mass7_err << "\n"
       << "7. RVUE: " << mass8 << " +/- " << mass8_err << "\n"
       << "8. CBAR (11k): " << mass9 << " +/- " << mass9_err << "\n"
       << "9. CBAR (3463): " << mass10 << " +/- " << mass10_err << "\n"
       << "10. CBAR (15k): " << mass11 << " +/- " << mass11_err << "\n";

  const double band_limit = MASS_ERR_BAND[0];

  int bin_indx = 0;
  for (int i = 0; i < nb_point; i ++) {
    bin_indx = gf_label -> GetXaxis() -> FindBin(LABEL_INDX[i]);
    cout << bin_indx << ", EXP_STR = " << EXP_STR[i] << endl;
    gf_label -> GetXaxis() -> SetBinLabel(bin_indx, EXP_STR[i]);
    gf_label1 -> GetXaxis() -> SetBinLabel(bin_indx, EXP1_STR[i]);
 
    //gf_label -> GetXaxis() -> SetBinLabel(18, "wfew");
  }

  gf_label2 -> GetXaxis() -> SetBinLabel(9, "H1");
  
  // plot
  mg -> Add(gf_band);
  
  TCanvas *cv_h1 = new TCanvas("cv_mass_H1", "cv_mass_H1", 0, 0, 700, 700);
  cv_h1 -> SetBottomMargin(0.15);
  cv_h1 -> SetLeftMargin(0.15);

  gf_label1 -> GetYaxis() -> SetNdivisions(502);
  gf_label1 -> GetYaxis() -> SetRangeUser(MASS_EXP1[0] - 40. * band_limit, MASS_EXP1[0] + 50. * band_limit);
  gf_label1 -> GetYaxis() -> SetTitleOffset(0.4);
  //gf_label1 -> GetYaxis() -> SetTitle("M_{#omega} [MeV/c^{2}]");
  gf_label1 -> GetYaxis() -> CenterTitle();
  gf_label1 -> GetYaxis() -> SetLabelSize(0.065);
  gf_label1 -> GetYaxis() -> SetTitleSize(0.03);
  
  gf_label1 -> Draw("AP");  
  mg -> Draw("C3");
  gf_exp1 -> Draw("P");

  TLegend * legd_gf1 = new TLegend(0.65, 0.8, .85, 0.9);

  SetLegend(legd_gf1);
  legd_gf1 -> SetNColumns(1);
  
  legd_gf1 -> AddEntry(gf_exp1, "H1", "lep");

  legd_gf1 -> Draw("Same");

  legtextsize(legd_gf1, 0.1);
  
  gPad -> Update();
  
  //
  TCanvas *cv = new TCanvas("cv_mass_compr", "cv_mass_compr", 0, 0, 1000, 700);
  //cv -> SetGrid();
  //cv -> SetLeftMargin(0.2);
  cv -> SetBottomMargin(0.19);

  /*
  TPad *p2 = new TPad("p2", "p2", 0., 0., 1., 0.35);
  p2 -> SetBottomMargin(0.2);
  //p2 -> SetGrid();
  p2 -> Draw();
  
  TPad *p1 = new TPad("p1", "p1", 0., 0.35, 1., 1.);
  p1 -> Draw();
  p1 -> SetBottomMargin(0.2);//0.007
  p1 -> cd();
  */
  
  char display1[50];

  TPaveText *pt1 = new TPaveText(0.11, 0.8, 0.65, 0.82, "NDC");
  PteAttr(pt1);
  pt1 -> SetTextSize(0.045);
  pt1 -> SetTextColor(kRed);
  pt1 -> AddText(Form("M_{#omega}=%0.5g#pm%0.2g MeV/c^{2} (this work)", mass2, mass2_err)); 

  TPaveText *pt2 = new TPaveText(0.11, 0.72, 0.65, 0.74, "NDC");
  PteAttr(pt2);
  pt2 -> SetTextSize(0.045);
  pt2 -> SetTextColor(kBlack);
  pt2 -> AddText(Form("M_{#omega}=%0.5g#pm%0.2g MeV/c^{2} (PDG)", mass_pdg, mass_pdg_err)); 

  gf_label -> GetXaxis() -> SetLabelSize(0.06);
  gf_label -> GetXaxis() -> SetLabelOffset(0.01);
  gf_label -> GetYaxis() -> SetNdivisions(512);
  //gf_label -> GetYaxis() -> SetRangeUser(MASS_BAND[0] - 40. * band_limit, MASS_BAND[0] + 30. * band_limit);
  gf_label -> GetYaxis() -> SetRangeUser(MASS_BAND[0] - 15. * band_limit, MASS_BAND[0] + 10. * band_limit);
  gf_label -> GetYaxis() -> SetTitleOffset(1.0);
  gf_label -> GetYaxis() -> SetTitle("M_{#omega} [MeV/c^{2}]");
  gf_label -> GetYaxis() -> CenterTitle();
  gf_label -> GetYaxis() -> SetLabelSize(0.035);
  gf_label -> GetYaxis() -> SetTitleSize(0.05);

  gf_label -> Draw("AP");
  mg -> Draw("C3");
  gf_exp2 -> Draw("P");
  gf_exp3 -> Draw("P");
  gf_exp4 -> Draw("P");
  gf_exp5 -> Draw("P");
  gf_exp6 -> Draw("P");
  gf_exp7 -> Draw("P");
  gf_exp8 -> Draw("P");
  gf_exp9 -> Draw("P");
  gf_exp10 -> Draw("P");
  gf_exp11 -> Draw("P");
  gf_exp12 -> Draw("P");
  gf_exp13 -> Draw("P");
  
  //pt1 -> Draw("Same");
  //pt2 -> Draw("Same");

  TLegend * legd_gf = new TLegend(0.65, 0.75, .85, 0.85);

  SetLegend(legd_gf);
  legd_gf -> SetNColumns(1);
  
  //
  //legd_gf -> AddEntry(gf_exp1,  "CMD-2", "lep");
  //legd_gf -> AddEntry(gf_exp2,  "RVUE", "lep");
  //legd_gf -> AddEntry(gf_exp3,  "KLOE-2", "lep");
  //legd_gf -> AddEntry(gf_exp4,  "SPEC", "lep");
  //legd_gf -> AddEntry(gf_exp5,  "CMD", "lep");

  legd_gf -> AddEntry(gf_band, "PDG", "lf");

  legd_gf -> Draw("Same");

  legtextsize(legd_gf, 0.07);

  
  
  //gPad -> Update();

  /*
  p2 -> cd();

  gf_label1 -> GetXaxis() -> SetLabelSize(0.15);
  gf_label1 -> GetXaxis() -> SetLabelOffset(0.02);
  //gf_label1 -> GetXaxis() -> SetRangeUser(0.15, 0.85);
  
  gf_label1 -> GetYaxis() -> SetNdivisions(502);
  gf_label1 -> GetYaxis() -> SetRangeUser(MASS_EXP1[0] - 80. * band_limit, MASS_EXP1[0] + 80. * band_limit);
  gf_label1 -> GetYaxis() -> SetTitleOffset(0.4);
  gf_label1 -> GetYaxis() -> SetTitle("M_{#omega} [MeV/c^{2}]");
  gf_label1 -> GetYaxis() -> CenterTitle();
  gf_label1 -> GetYaxis() -> SetLabelSize(0.08);
  gf_label1 -> GetYaxis() -> SetTitleSize(0.13);
  
  gf_label1 -> Draw("AP");  
  gf_exp1 -> Draw("P");

  gPad -> Update();
  */
  
  // save
  cv -> SaveAs(outputFile + "/mass_compr.pdf");
  //cv_h1 -> SaveAs(outputFile + "/mass_h1.pdf");
  
  return 0;
}
