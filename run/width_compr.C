#include "../header/plot.h"
#include "../header/compr.h"

int width_compr(){

  gErrorIgnoreLevel = kError;
  TGaxis::SetMaxDigits(4);
  //gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  // fill list
  const int nb_point = 4;
  double WIDTH_BAND[nb_point];
  double WIDTH_ERR_BAND[nb_point];
  double XLIST[nb_point], XLIST_ERR[nb_point];
  
  // PDG
  const double width_pdg = 8.68;
  const double width_pdg_err = 0.13;

  for (int i = 0; i < nb_point; i ++) {

    WIDTH_BAND[i] = width_pdg;
    WIDTH_ERR_BAND[i] = width_pdg_err;

    XLIST[i] = i;
    XLIST_ERR[i] = 0.;
    
  }

  // KLOE-2 ISR (this work)
  const double width1 = 8.73;
  const double width1_err = 0.11;
  const double width1_syst_err_plus = 0.12;
  const double width1_syst_err_minus = 0.16;

  double width1_exl = 0.; // lower error for x
  double width1_exh = 0.; // higher error for x
  double width1_eyl = TMath::Sqrt(width1_err * width1_err + width1_syst_err_minus * width1_syst_err_minus); // lower error for y
  double width1_eyh = TMath::Sqrt(width1_err * width1_err + width1_syst_err_plus * width1_syst_err_plus); // higher error for y

  const double WIDTH_EXP1[1] = {width1};
  const double WIDTH_EXP1_ERR[1] = {width1_err};
  const double XLIST_EXP1[1] = {0.1};
  const double XLIST_EXP1_ERR[1] = {0};
  
  const double WIDTH1_EXL[1] = {width1_exl};
  const double WIDTH1_EXH[1] = {width1_exh};
  const double WIDTH1_EYL[1] = {width1_eyl};
  const double WIDTH1_EYH[1] = {width1_eyh};

  // BaBar 2pi  arXiv:1205.2228
  const double width2 = 8.68;
  const double width2_err = TMath::Sqrt(0.36 * 0.36 + 0.27 * 0.27);

  const double WIDTH_EXP2[1] = {width2};
  const double WIDTH_EXP2_ERR[1] = {width2_err};
  const double XLIST_EXP2[1] = {0.2};
  const double XLIST_EXP2_ERR[1] = {0};

  cout << "width BaBar = " << width2 << "+/-" << width2_err << endl;
  
  // RVUE
  const double width3 = 8.68;
  const double width3_err = TMath::Sqrt(0.04 * 0.04 + 0.15 * 0.15);

  const double WIDTH_EXP3[1] = {width3};
  const double WIDTH_EXP3_ERR[1] = {width3_err};
  const double XLIST_EXP3[1] = {0.3};
  const double XLIST_EXP3_ERR[1] = {0};

  // CMD2
  const double width4 = 8.68;
  const double width4_err = TMath::Sqrt(0.23 * 0.23 + 0.10 * 0.10);

  const double WIDTH_EXP4[1] = {width4};
  const double WIDTH_EXP4_ERR[1] = {width3_err};
  const double XLIST_EXP4[1] = {0.4}; 
  const double XLIST_EXP4_ERR[1] = {0};

  // label
  const double LABEL_INDX[nb_point] = {XLIST_EXP1[0], XLIST_EXP2[0], XLIST_EXP3[0], XLIST_EXP4[0]};
  const double LABEL_GF[nb_point] = {0., 0., 0., 0.};
  const TString EXP_STR[nb_point] = {"KLOE-2", "BABR", "RVUE", "CMD-2"};

  // graphs
  TGraphAsymmErrors *gf_exp1 = new TGraphAsymmErrors(1, XLIST_EXP1, WIDTH_EXP1, WIDTH1_EXL, WIDTH1_EXH, WIDTH1_EYL, WIDTH1_EYH);
  //TGraphErrors *gf_exp1 = new TGraphErrors(1, XLIST_EXP1, WIDTH_EXP1, XLIST_EXP1_ERR, WIDTH_EXP1_ERR);
  gf_exp1 -> SetName("gf_exp1");
  gf_exp1 -> SetMarkerStyle(21);
  gf_exp1 -> SetMarkerSize(1.3);
  gf_exp1 -> SetMarkerColor(kRed);

  TGraphErrors *gf_exp2 = new TGraphErrors(1, XLIST_EXP2, WIDTH_EXP2, XLIST_EXP2_ERR, WIDTH_EXP2_ERR);
  gf_exp2 -> SetName("gf_exp2");
  gf_exp2 -> SetMarkerStyle(8);
  gf_exp2 -> SetMarkerSize(1.3);
  gf_exp2 -> SetMarkerColor(kBlue);

  TGraphErrors *gf_exp3 = new TGraphErrors(1, XLIST_EXP3, WIDTH_EXP3, XLIST_EXP3_ERR, WIDTH_EXP3_ERR);
  gf_exp3 -> SetName("gf_exp3");
  gf_exp3 -> SetMarkerStyle(22);
  gf_exp3 -> SetMarkerSize(1.3);
  gf_exp3 -> SetMarkerColor(kBlack);

  TGraphErrors *gf_exp4 = new TGraphErrors(1, XLIST_EXP4, WIDTH_EXP4, XLIST_EXP4_ERR, WIDTH_EXP4_ERR);
  gf_exp4 -> SetName("gf_exp4");
  gf_exp4 -> SetMarkerStyle(22);
  gf_exp4 -> SetMarkerSize(1.3);
  gf_exp4 -> SetMarkerColor(kBlack);

  TGraphErrors *gf_band = new TGraphErrors(nb_point, XLIST, WIDTH_BAND, XLIST_ERR, WIDTH_ERR_BAND);
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
  

  const double band_limit = WIDTH_ERR_BAND[0];

  int bin_indx = 0;
  for (int i = 0; i < nb_point; i ++) {
    bin_indx = gf_label -> GetXaxis() -> FindBin(LABEL_INDX[i]);
    cout << bin_indx << ", EXP_STR = " << EXP_STR[i] << endl;
    gf_label -> GetXaxis() -> SetBinLabel(bin_indx, EXP_STR[i]);
  }

  // plot
  TCanvas *cv = new TCanvas("cv_width_compr", "cv_width_compr", 0, 0, 1000, 700);
  //cv -> SetGrid();
  //cv -> SetLeftMargin(0.2);
  cv -> SetBottomMargin(0.17);

  TPaveText *pt1 = new TPaveText(0.11, 0.8, 0.65, 0.82, "NDC");
  PteAttr(pt1);
  pt1 -> SetTextSize(0.045);
  pt1 -> SetTextColor(kRed);
  pt1 -> AddText(Form("#Gamma_{#omega}=%0.5g#pm%0.2g MeV/c^{2} (this work)", width1, width1_err)); 

  TPaveText *pt2 = new TPaveText(0.11, 0.72, 0.65, 0.74, "NDC");
  PteAttr(pt2);
  pt2 -> SetTextSize(0.045);
  pt2 -> SetTextColor(kBlack);
  pt2 -> AddText(Form("#Gamma_{#omega}=%0.5g#pm%0.2g MeV/c^{2} (PDG)", width_pdg, width_pdg_err)); 

  gf_label -> GetXaxis() -> SetLabelSize(0.06);
  gf_label -> GetXaxis() -> SetLabelOffset(0.01);
  gf_label -> GetYaxis() -> SetNdivisions(512);
  gf_label -> GetYaxis() -> SetRangeUser(WIDTH_BAND[0] - 5. * band_limit, WIDTH_BAND[0] + 5. * band_limit);
  gf_label -> GetYaxis() -> SetTitleOffset(0.9);
  gf_label -> GetYaxis() -> SetTitle("#Gamma_{#omega} [MeV]");
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
  cv -> SaveAs(outputFile + "/width_compr.pdf");
  

  return 0;

}
