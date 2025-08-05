#include "../plot.h"
#include "compr.h"
#include "../crx3pi/crx3pi.h"

int branch_compr(){

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

  double BB_fit = 0., BB_err_fit = 0.;

  for (Int_t irow = 0; irow < TCRX3PI -> GetEntries(); irow++) {// start fill
    
    TCRX3PI -> GetEntry(irow); //cout << irow << endl;

    BB_fit = TCRX3PI -> GetLeaf("Br_BB_fit") -> GetValue(0);
    BB_err_fit = TCRX3PI -> GetLeaf("Br_BB_err_fit") -> GetValue(0);
    
  }
  
  cout << "Plotting omega branc comparison ...\n";
  //pdg, 6.38 +/- 0.10
  
  const int para_size = 1;
  const double PARA_PDG[para_size] = {6.41};
  const double PARA_ERR_PDG[para_size] = {0.13};
  
  const int nb_point = 12;
  
  double BRANCH_BAND[nb_point];
  double BRANCH_ERR_BAND[nb_point];
  double XLIST[nb_point], XLIST_ERR[nb_point];

  double BRANCH_EXP1[1], BRANCH_EXP1_ERR[1];
  double XLIST_EXP1[1], XLIST_EXP1_ERR[1];

  double BRANCH_EXP2[1], BRANCH_EXP2_ERR[1];
  double XLIST_EXP2[1], XLIST_EXP2_ERR[1];

  double BRANCH_EXP3[1], BRANCH_EXP3_ERR[1];
  double XLIST_EXP3[1], XLIST_EXP3_ERR[1];

  double BRANCH_EXP4[1], BRANCH_EXP4_ERR[1];
  double XLIST_EXP4[1], XLIST_EXP4_ERR[1];

  double BRANCH_EXP5[1], BRANCH_EXP5_ERR[1];
  double XLIST_EXP5[1], XLIST_EXP5_ERR[1];

  double BRANCH_EXP6[1], BRANCH_EXP6_ERR[1];
  double XLIST_EXP6[1], XLIST_EXP6_ERR[1];

  
  // CMD2, 6.24 +/- 0.11 +/- 0.08
  BRANCH_EXP1[0] = 6.24, BRANCH_EXP1_ERR[0] = TMath::Sqrt(0.11 * 0.11 + 0.08 * 0.08);
  XLIST_EXP1[0] = 1.1, XLIST_EXP1_ERR[0] = 0;

  // BABR, 6.70 +/- 0.06 +/- 0.27
  BRANCH_EXP2[0] = 6.70, BRANCH_EXP2_ERR[0] = TMath::Sqrt(0.06 * 0.06 + 0.27 * 0.27);
  XLIST_EXP2[0] = 1.3, XLIST_EXP2_ERR[0] = 0;

  // KLOE-2, 6.54398 +/- 0.0636307 (fine binned)
  double kloe_value = BB_fit * 1e5, kloe_err = BB_err_fit * 1e5;
  
  BRANCH_EXP3[0] = kloe_value, BRANCH_EXP3_ERR[0] = TMath::Sqrt(kloe_err * kloe_err);
  XLIST_EXP3[0] = 1.5, XLIST_EXP3_ERR[0] = 0;

  // RVUE, 6.74 +/- 0.04 +/- 0.24
  BRANCH_EXP4[0] = 6.74, BRANCH_EXP4_ERR[0] = TMath::Sqrt(0.04 * 0.04 + 0.24 * 0.24);
  XLIST_EXP4[0] = 1.7, XLIST_EXP4_ERR[0] = 0;

  // ND, 6.37 +/- 0.35
  BRANCH_EXP5[0] = 6.37, BRANCH_EXP5_ERR[0] = TMath::Sqrt(0.35 * 0.35);
  XLIST_EXP5[0] = 1.9, XLIST_EXP5_ERR[0] = 0;

  // BaBar, 6.45 +/- 0.24
  BRANCH_EXP6[0] = 6.45, BRANCH_EXP6_ERR[0] = TMath::Sqrt(0.24 * 0.24);
  XLIST_EXP6[0] = 1.9, XLIST_EXP6_ERR[0] = 0;
 
  const int exp_nb = 5;
  const double LABEL_INDX[exp_nb] = {XLIST_EXP1[0], XLIST_EXP2[0], XLIST_EXP3[0], XLIST_EXP4[0], XLIST_EXP5[0]};
  const double LABEL_GF[exp_nb] = {0., 0., 0., 0., 0.};
  const TString EXP_STR[exp_nb] = {"CMD2", "BABR", "KLOE-2", "RVUE", "ND"};

  for (int i = 0; i < nb_point; i ++) {

    BRANCH_BAND[i] = PARA_PDG[0];
    BRANCH_ERR_BAND[i] = PARA_ERR_PDG[0];

    XLIST[i] = i;
    XLIST_ERR[i] = 0.;
    
    cout << "BRANCH_BAND = " << BRANCH_BAND[i] << " +/- " << BRANCH_ERR_BAND[i] << "\n";
	  
  }

  // test combine two measurements: CMD2 and RVUE
  const double EXPVALUE[exp_nb - 1] = {BRANCH_EXP1[0], BRANCH_EXP2[0], BRANCH_EXP4[0], BRANCH_EXP5[0]};
  const double ERROR[exp_nb - 1] = {BRANCH_EXP1_ERR[0], BRANCH_EXP2_ERR[0], BRANCH_EXP4_ERR[0], BRANCH_EXP5_ERR[0]};

  double weight_tmp = 0., weight_sum = 0.;
  double exp_weight_sum = 0.;
  double exp_weight_average = 0.;

  cout << "BB [10^{5}], error includes scale factor of 1.2\n";
  
  for (int i = 0; i < exp_nb - 1; i ++) {

    weight_tmp = GetWeight(ERROR[i] * 1.2);
    weight_sum += weight_tmp;

    exp_weight_sum += weight_tmp * EXPVALUE[i];
    exp_weight_average = exp_weight_sum / weight_sum;

    cout << EXPVALUE[i] << " +/- " << ERROR[i] << ", weight_tmp = " << weight_tmp << ", weight_sum = " << weight_sum << "\n"
	 << "exp_weight_average = " << exp_weight_average << "\n"
	 << "1/sqrt(weight_sum) = " << 1. / TMath::Sqrt(weight_sum) << "\n";

  }
  
  //
  TGraph *gf_label = new TGraph(exp_nb, LABEL_INDX, LABEL_GF);
  gf_label -> SetName("gf_label");
  
  TGraphErrors *gf_exp1 = new TGraphErrors(1, XLIST_EXP1, BRANCH_EXP1, XLIST_EXP1_ERR, BRANCH_EXP1_ERR);
  gf_exp1 -> SetName("gf_exp1");
  gf_exp1 -> SetMarkerStyle(20);
  gf_exp1 -> SetMarkerSize(1.1);

  TGraphErrors *gf_exp2 = new TGraphErrors(1, XLIST_EXP2, BRANCH_EXP2, XLIST_EXP2_ERR, BRANCH_EXP2_ERR);
  gf_exp2 -> SetName("gf_exp2");
  gf_exp2 -> SetMarkerStyle(21);
  gf_exp2 -> SetMarkerSize(1.1);
  gf_exp2 -> SetMarkerColor(kBlue);

  TGraphErrors *gf_exp3 = new TGraphErrors(1, XLIST_EXP3, BRANCH_EXP3, XLIST_EXP3_ERR, BRANCH_EXP3_ERR);
  gf_exp3 -> SetName("gf_exp3");
  gf_exp3 -> SetMarkerStyle(21);
  gf_exp3 -> SetMarkerSize(1.1);
  gf_exp3 -> SetMarkerColor(kRed);

  TGraphErrors *gf_exp4 = new TGraphErrors(1, XLIST_EXP4, BRANCH_EXP4, XLIST_EXP4_ERR, BRANCH_EXP4_ERR);
  gf_exp4 -> SetName("gf_exp4");
  gf_exp4 -> SetMarkerStyle(23);
  gf_exp4 -> SetMarkerSize(1.1);

  TGraphErrors *gf_exp5 = new TGraphErrors(1, XLIST_EXP5, BRANCH_EXP5, XLIST_EXP5_ERR, BRANCH_EXP5_ERR);
  gf_exp5 -> SetName("gf_exp5");
  gf_exp5 -> SetMarkerStyle(24);
  gf_exp5 -> SetMarkerSize(1.1);

  TGraphErrors *gf_exp6 = new TGraphErrors(1, XLIST_EXP6, BRANCH_EXP6, XLIST_EXP6_ERR, BRANCH_EXP6_ERR);
  gf_exp6 -> SetName("gf_exp6");
  gf_exp6 -> SetMarkerStyle(25);
  gf_exp6 -> SetMarkerSize(1.1);

  TGraphErrors *gf_band = new TGraphErrors(nb_point, XLIST, BRANCH_BAND, XLIST_ERR, BRANCH_ERR_BAND);
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
  const double band_limit = 12. * BRANCH_ERR_BAND[0];
  gf_label -> RemovePoint(1);
  gf_label -> RemovePoint(2);
  gf_label -> RemovePoint(3);
  gf_label -> RemovePoint(4);
  gf_label -> RemovePoint(5);

  gf_label -> GetXaxis() -> SetRangeUser(0., 10.);
  gf_label -> GetYaxis() -> SetRangeUser(BRANCH_BAND[0] - band_limit, BRANCH_BAND[0] + band_limit);
  gf_label -> GetYaxis() -> SetTitleOffset(1.0);
  gf_label -> GetYaxis() -> SetTitle("B(e^{+}e^{-}#rightarrow#omega)B(#omega#rightarrow3#pi)#times10^{5}");
  gf_label -> GetYaxis() -> CenterTitle();
  gf_label -> GetYaxis() -> SetLabelSize(0.025);
  gf_label -> GetYaxis() -> SetTitleSize(0.04);

  int bin_indx = 0;
  for (int i = 0; i < exp_nb; i ++) {
    bin_indx = gf_label -> GetXaxis() -> FindBin(LABEL_INDX[i]);
    cout << bin_indx << endl;
    gf_label -> GetXaxis() -> SetBinLabel(bin_indx, EXP_STR[i]);
    //gf_label -> GetXaxis() -> SetBinLabel(18, "wfew");
  }
 


  TCanvas *cv = new TCanvas("cv", "", 0, 0, 700, 700);

  char display1[50];

  TPaveText *pt1 = new TPaveText(0.4, 0.84, 0.85, 0.86, "NDC");
  PteAttr(pt1);
  pt1 -> SetTextColor(kRed);
  pt1 -> AddText(Form("B_{ee}B_{3#pi}=%0.3g#pm%0.1g", kloe_value, kloe_err)); 

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
  //legd_gf -> AddEntry(gf_exp2,  "BABR", "lep");
  //legd_gf -> AddEntry(gf_exp3,  "KLOE-2", "lep");
  //legd_gf -> AddEntry(gf_exp4,  "RVUE", "lep");
  //legd_gf -> AddEntry(gf_exp5,  "CMD", "lep");

  legd_gf -> AddEntry(gf_band,  "PDG", "lf");

  legd_gf -> Draw("Same");

  legtextsize(legd_gf, 0.03);

  gPad -> Update();

  // save
  cv -> SaveAs("plots/branch_compr.pdf");
  
  return 0;
  
}

double GetWeight(double error) {

  double weight = TMath::Power(1. / error, 2);
    
  return weight;
  
}
