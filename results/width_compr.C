#include "../plot.h"
#include "compr.h"
#include "../crx3pi/crx3pi.h"

int width_compr(){

  gErrorIgnoreLevel = kError;
  TGaxis::SetMaxDigits(3);
  //gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);  

  // get kloe2 results

  //output8373_finebinned
  const double sample_frac = Lumi_int / 1724470 * 100.;
  cout << "Getting input file from " << infile_tmp << "\n"
       << "Lumi. int [%] = " << sample_frac << "\n";
  
  TFile *intree = new TFile(infile_tmp);

  TIter next_tree(intree -> GetListOfKeys());
  
  TString objnm_tree, classnm_tree;
  
  int i = 0;
  
  while ( (key = (TKey *) next_tree() ) ) {// start tree while lop
    
    i ++;
    
    objnm_tree   =  key -> GetName();
    classnm_tree = key -> GetClassName();
    key -> GetSeekKey();
    
    //cout << " tree/histo" << i << ": classnm = " << classnm_tree << ", objnm = " << objnm_tree << endl;
    
  }

  double Gam_omega_fit = 0., Gam_omega_err_fit = 0.;
  
  for (Int_t irow = 0; irow < TCRX3PI -> GetEntries(); irow++) {// start fill
    
    TCRX3PI -> GetEntry(irow); //cout << irow << endl;

    Gam_omega_fit = TCRX3PI -> GetLeaf("Br_Gam_omega_fit") -> GetValue(0);
    Gam_omega_err_fit = TCRX3PI -> GetLeaf("Br_Gam_omega_err_fit") -> GetValue(0);
  
  }

  cout << Gam_omega_fit << endl;
  
  cout << "Plotting omega width comparison ...\n";
  
  const int para_size = 1;
  const double PARA_PDG[para_size] = {8.68};
  const double PARA_ERR_PDG[para_size] = {0.13};
  
  const int nb_point = 10;
  
  double WIDTH_BAND[nb_point];
  double WIDTH_ERR_BAND[nb_point];
  double XLIST[nb_point], XLIST_ERR[nb_point];
  
  double WIDTH_EXP1[1], WIDTH_EXP1_ERR[1];
  double XLIST_EXP1[1], XLIST_EXP1_ERR[1];

  double WIDTH_EXP2[1], WIDTH_EXP2_ERR[1];
  double XLIST_EXP2[1], XLIST_EXP2_ERR[1];

  double WIDTH_EXP3[1], WIDTH_EXP3_ERR[1];
  double XLIST_EXP3[1], XLIST_EXP3_ERR[1];

  double WIDTH_EXP4[1], WIDTH_EXP4_ERR[1];
  double XLIST_EXP4[1], XLIST_EXP4_ERR[1];

  double WIDTH_EXP5[1], WIDTH_EXP5_ERR[1];
  double XLIST_EXP5[1], XLIST_EXP5_ERR[1];

  
  // CMD2, 8.68 +/- 0.23 +/- 0.10
  WIDTH_EXP1[0] = 8.68, WIDTH_EXP1_ERR[0] = TMath::Sqrt(0.23 * 0.23 + 0.10 * 0.10);
  XLIST_EXP1[0] = 1.1, XLIST_EXP1_ERR[0] = 0;

  // RVUE, 8.68 +/- 0.04 +/- 0.15
  WIDTH_EXP2[0] = 8.68, WIDTH_EXP2_ERR[0] = TMath::Sqrt(0.04 * 0.04 + 0.15 * 0.15);
  XLIST_EXP2[0] = 1.3, XLIST_EXP2_ERR[0] = 0;

  // KLOE-2, 8.49671 +/- 0.154557
  double kloe_value = Gam_omega_fit, kloe_err = Gam_omega_err_fit;
  WIDTH_EXP3[0] = kloe_value, WIDTH_EXP3_ERR[0] = TMath::Sqrt(kloe_err * kloe_err);
  XLIST_EXP3[0] = 1.5, XLIST_EXP3_ERR[0] = 0;

  // SPEC, 8.2 +/- 0.3
  WIDTH_EXP4[0] = 8.2, WIDTH_EXP4_ERR[0] = TMath::Sqrt(0.3 * 0.3);
  XLIST_EXP4[0] = 1.7, XLIST_EXP4_ERR[0] = 0;

  // BABR, 8.13 +/- 0.45
  WIDTH_EXP5[0] = 8.13, WIDTH_EXP5_ERR[0] = TMath::Sqrt(0.45 * 0.45);
  XLIST_EXP5[0] = 1.9, XLIST_EXP5_ERR[0] = 0;
  
  const double LABEL_INDX[5] = {XLIST_EXP1[0], XLIST_EXP2[0], XLIST_EXP3[0], XLIST_EXP4[0], XLIST_EXP5[0]};
  const double LABEL_GF[5] = {0., 0., 0., 0., 0.};
  const TString EXP_STR[5] = {"CMD2", "RVUE", "KLOE-2", "SPEC^{*}", "BABR^{*}"};

  for (int i = 0; i < nb_point; i ++) {

    WIDTH_BAND[i] = PARA_PDG[0];
    WIDTH_ERR_BAND[i] = PARA_ERR_PDG[0];

    XLIST[i] = i;
    XLIST_ERR[i] = 0.;
    
    //cout << "WIDTH_BAND = " << WIDTH_BAND[i] << " +/- " << WIDTH_ERR_BAND[i] << "\n";
	  
  }

  // test combine two measurements: CMD2 and RVUE
  double measure1 = WIDTH_EXP1[0], measure1_staterr = 0.23, measure1_systerr = 0.10;
  double measure1_errsum = TMath::Sqrt(measure1_staterr * measure1_staterr + measure1_systerr * measure1_systerr);
  double weight1 = TMath::Power(1. / measure1_errsum, 2);
  
  double measure2 = WIDTH_EXP2[0], measure2_staterr = 0.04, measure2_systerr = 0.15;
  double measure2_errsum = TMath::Sqrt(measure2_staterr * measure2_staterr + measure2_systerr * measure2_systerr);
  double weight2 = TMath::Power(1. / measure2_errsum, 2);

  double weight1_ratio = weight1 / (weight1 + weight2);
  double weight2_ratio = weight2 / (weight1 + weight2);
  double weight_sum = weight1 + weight2;
  
  double measure_averge = (measure1 * weight1 + measure2 * weight2) / (weight1 + weight2);
  double measure_averge_err = TMath::Sqrt(weight1_ratio * weight1_ratio * measure1_errsum * measure1_errsum + weight2_ratio * weight2_ratio * measure2_errsum * measure2_errsum);
  //(measure1_errsum * weight1 + measure2_errsum * weight2) / (weight1 + weight2);
  

  cout << "measure1 = " << measure1 << " +/-" << measure1_errsum << ", weight1 = " << weight1 << "\n";
  cout << "measure2 = " << measure2 << " +/-" << measure2_errsum << ", weight2 = " << weight2 << "\n";
  cout << "weight1_ratio = " << weight1_ratio << "\n";
  cout << "weight2_ratio = " << weight2_ratio << "\n";
  cout << "average = " << measure_averge << " +/- " << measure_averge_err << "\n";
  cout << "1/sqrt(weight_sum) = " << 1. / TMath::Sqrt(weight_sum) << "\n";

  //

  TGraph *gf_label = new TGraph(5, LABEL_INDX, LABEL_GF);
  gf_label -> SetName("gf_label");
  
  TGraphErrors *gf_exp1 = new TGraphErrors(1, XLIST_EXP1, WIDTH_EXP1, XLIST_EXP1_ERR, WIDTH_EXP1_ERR);
  gf_exp1 -> SetName("gf_exp1");
  gf_exp1 -> SetMarkerStyle(20);
  gf_exp1 -> SetMarkerSize(1.1);

  TGraphErrors *gf_exp2 = new TGraphErrors(1, XLIST_EXP2, WIDTH_EXP2, XLIST_EXP2_ERR, WIDTH_EXP2_ERR);
  gf_exp2 -> SetName("gf_exp2");
  gf_exp2 -> SetMarkerStyle(22);
  gf_exp2 -> SetMarkerSize(1.1);

  TGraphErrors *gf_exp3 = new TGraphErrors(1, XLIST_EXP3, WIDTH_EXP3, XLIST_EXP3_ERR, WIDTH_EXP3_ERR);
  gf_exp3 -> SetName("gf_exp3");
  gf_exp3 -> SetMarkerStyle(21);
  gf_exp3 -> SetMarkerSize(1.1);
  gf_exp3 -> SetMarkerColor(kRed);

  TGraphErrors *gf_exp4 = new TGraphErrors(1, XLIST_EXP4, WIDTH_EXP4, XLIST_EXP4_ERR, WIDTH_EXP4_ERR);
  gf_exp4 -> SetName("gf_exp3");
  gf_exp4 -> SetMarkerStyle(23);
  gf_exp4 -> SetMarkerSize(1.1);

  TGraphErrors *gf_exp5 = new TGraphErrors(1, XLIST_EXP5, WIDTH_EXP5, XLIST_EXP5_ERR, WIDTH_EXP5_ERR);
  gf_exp5 -> SetName("gf_exp3");
  gf_exp5 -> SetMarkerStyle(24);
  gf_exp5 -> SetMarkerSize(1.1);

  TGraphErrors *gf_band = new TGraphErrors(nb_point, XLIST, WIDTH_BAND, XLIST_ERR, WIDTH_ERR_BAND);
  gf_band -> SetName("gf_band");
  gf_band -> SetMarkerStyle(20);
  gf_band -> SetMarkerSize(1.1);
  gf_band -> SetLineColor(kGreen);
  gf_band -> SetLineWidth(2);
  gf_band -> SetFillColor(kBlue);
  gf_band -> SetFillStyle(3002);

  TMultiGraph *mg = new TMultiGraph();
  mg -> SetTitle("Exclusion graphs");
  
  //
  const double band_limit = 12. * WIDTH_ERR_BAND[0];
  //gf_exp1 -> GetXaxis() -> SetLabelOffset(999);
  //gf_exp1 -> GetXaxis() -> SetLabelSize(0);
  //gf_exp1 -> GetXaxis() -> SetTickLength(0);
  
  gf_label -> RemovePoint(1);
  gf_label -> RemovePoint(2);
  gf_label -> RemovePoint(3);
  gf_label -> RemovePoint(4);
  gf_label -> RemovePoint(5);

  gf_label -> GetXaxis() -> SetRangeUser(0., 10.);
  gf_label -> GetYaxis() -> SetRangeUser(WIDTH_BAND[0] - band_limit, WIDTH_BAND[0] + band_limit);
  gf_label -> GetYaxis() -> SetTitleOffset(1.2);
  gf_label -> GetYaxis() -> SetTitle("#Gamma_{#omega} [MeV/c^{2}]");
  gf_label -> GetYaxis() -> CenterTitle();
  gf_label -> GetYaxis() -> SetLabelSize(0.025);
  gf_label -> GetYaxis() -> SetTitleSize(0.04);

  
  int bin_indx = 0;
  for (int i = 0; i < 5; i ++) {
    bin_indx = gf_label -> GetXaxis() -> FindBin(LABEL_INDX[i]);
    //cout << bin_indx << endl;
    gf_label -> GetXaxis() -> SetBinLabel(bin_indx, EXP_STR[i]);
    //gf_label -> GetXaxis() -> SetBinLabel(18, "wfew");
  }


  TCanvas *cv = new TCanvas("cv", "", 0, 0, 700, 700);

  char display1[50];

  TPaveText *pt1 = new TPaveText(0.35, 0.84, 0.85, 0.86, "NDC");
  PteAttr(pt1);
  pt1 -> SetTextColor(kRed);
  pt1 -> AddText(Form("#Gamma_{#omega}=%0.3g#pm%0.2g MeV/c^{2}", kloe_value, kloe_err)); 

  gf_label -> Draw("AP");
  gf_exp1 -> Draw("P");
  gf_exp2 -> Draw("P");
  gf_exp3 -> Draw("P");
  gf_exp4 -> Draw("P");
  gf_exp5 -> Draw("P");
 
  mg -> Add(gf_band);
  mg -> Draw("C3");

  pt1 -> Draw("Same");
  
  TLegend * legd_gf = new TLegend(0.15, 0.8, .5, 0.9);

  SetLegend(legd_gf);
  legd_gf -> SetNColumns(1);
  
  //
  //legd_gf -> AddEntry(gf_exp1,  "CMD2", "lep");
  //legd_gf -> AddEntry(gf_exp2,  "RVUE", "lep");
  //legd_gf -> AddEntry(gf_exp3,  "KLOE-2", "lep");
  //legd_gf -> AddEntry(gf_exp4,  "SPEC", "lep");
  //legd_gf -> AddEntry(gf_exp5,  "CMD", "lep");

  legd_gf -> AddEntry(gf_band,  "PDG", "lf");

  legd_gf -> Draw("Same");

  legtextsize(legd_gf, 0.03);

  gPad -> Update();

  
  // save
  cv -> SaveAs("plots/width_compr.pdf");
  
  return 0;
  
}
