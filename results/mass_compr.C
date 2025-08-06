#include "plot.h"
#include "compr.h"
//#include "crx3pi.h"

int mass_compr(){

  gErrorIgnoreLevel = kError;
  TGaxis::SetMaxDigits(4);
  //gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
 
  // get kloe2 results

  
  //output8373_finebinned
  /*
  const double sample_frac = Lumi_int / 1724470 * 100.;
  cout << "Getting input file from " << infile_tmp << "\n"
       << "Lumi. int [%] = " << sample_frac << "\n";
  */
  
  TFile *intree = new TFile(infile_tmp);

  TIter next_tree(intree -> GetListOfKeys());
  
  TString objnm_tree, classnm_tree;
  
  int i = 0;
  TKey *key;
  
  while ( (key = (TKey *) next_tree() ) ) {// start tree while lop
    
    i ++;
    
    objnm_tree   =  key -> GetName();
    classnm_tree = key -> GetClassName();
    key -> GetSeekKey();
    
    //cout << " tree/histo" << i << ": classnm = " << classnm_tree << ", objnm = " << objnm_tree << endl;
    
  }

  double Mass_omega_fit = 0., Mass_omega_err_fit = 0.;
    
  for (Int_t irow = 0; irow < TCRX3PI -> GetEntries(); irow++) {// start fill
    
    TCRX3PI -> GetEntry(irow); //cout << irow << endl;

    Mass_omega_fit = TCRX3PI -> GetLeaf("Br_Mass_omega_fit") -> GetValue(0);
    Mass_omega_err_fit = TCRX3PI -> GetLeaf("Br_Mass_omega_err_fit") -> GetValue(0);

  }
 
  cout << "Plotting omega mass comparison ...\n";
  //pdg, 782.65 +/- 0.12 MeV [world average]
  
  const int para_size = 1;
  const double PARA_PDG[para_size] = {782.66};
  const double PARA_ERR_PDG[para_size] = {0.13};
  
  const int nb_point = 12;
  
  double MASS_BAND[nb_point];
  double MASS_ERR_BAND[nb_point];
  double XLIST[nb_point], XLIST_ERR[nb_point];

  double MASS_EXP1[1], MASS_EXP1_ERR[1];
  double XLIST_EXP1[1], XLIST_EXP1_ERR[1];

  double MASS_EXP2[1], MASS_EXP2_ERR[1];
  double XLIST_EXP2[1], XLIST_EXP2_ERR[1];

  double MASS_EXP3[1], MASS_EXP3_ERR[1];
  double XLIST_EXP3[1], XLIST_EXP3_ERR[1];

  double MASS_EXP4[1], MASS_EXP4_ERR[1];
  double XLIST_EXP4[1], XLIST_EXP4_ERR[1];

  double MASS_EXP5[1], MASS_EXP5_ERR[1];
  double XLIST_EXP5[1], XLIST_EXP5_ERR[1];

  double MASS_EXP6[1], MASS_EXP6_ERR[1];
  double XLIST_EXP6[1], XLIST_EXP6_ERR[1];

  
  // AKHMETSHIN 05, CMD2, 783.20 +/- 0.13 +/- 0.16
  MASS_EXP1[0] = 783.20, MASS_EXP1_ERR[0] = TMath::Sqrt(0.13 * 0.13 + 0.16 * 0.16);
  XLIST_EXP1[0] = 0.1, XLIST_EXP1_ERR[0] = 0;

  // AKHMETSHIN 04, CMD2, 782.68 +/- 0.09 +/- 0.04
  MASS_EXP2[0] = 782.68, MASS_EXP2_ERR[0] = TMath::Sqrt(0.09 * 0.09 + 0.04 * 0.04);
  XLIST_EXP2[0] = 0.3, XLIST_EXP2_ERR[0] = 0;

  // KLOE-2, 782.828 +/- 0.0511992  
  double kloe_value = Mass_omega_fit, kloe_err = Mass_omega_err_fit;
  MASS_EXP3[0] = kloe_value, MASS_EXP3_ERR[0] = TMath::Sqrt(kloe_err * kloe_err);
  XLIST_EXP3[0] = .4, XLIST_EXP3_ERR[0] = 0;

  // RVUE, 782.79 +/- 0.08 +/- 0.09
  MASS_EXP4[0] = 782.79, MASS_EXP4_ERR[0] = TMath::Sqrt(0.08 * 0.08 + 0.09 * 0.09);
  XLIST_EXP4[0] = .5, XLIST_EXP4_ERR[0] = 0;

  // CBAR, 781.96 +/- 0.13 +/- 0.17
  MASS_EXP5[0] = 781.96, MASS_EXP5_ERR[0] = TMath::Sqrt(0.13 * 0.13 + 0.17 * 0.17);
  XLIST_EXP5[0] = .6, XLIST_EXP5_ERR[0] = 0;

  // BABR, 782.45 +/- 0.22
  MASS_EXP6[0] = 782.45, MASS_EXP6_ERR[0] = TMath::Sqrt(0.22 * 0.22 + 0.22 * 0.22);
  XLIST_EXP6[0] = .2, XLIST_EXP6_ERR[0] = 0;

  const int exp_nb = 6;
  const double LABEL_INDX[exp_nb] = {XLIST_EXP1[0], XLIST_EXP2[0], XLIST_EXP3[0], XLIST_EXP4[0], XLIST_EXP5[0], XLIST_EXP6[0]};
  const double LABEL_GF[exp_nb] = {0., 0., 0., 0., 0., 0.};
  const TString EXP_STR[exp_nb] = {"CMD2", "CMD2", "KLOE-2", "RVUE", "CBAR", "BABR^{*}"};

  for (int i = 0; i < nb_point; i ++) {

    MASS_BAND[i] = PARA_PDG[0];
    MASS_ERR_BAND[i] = PARA_ERR_PDG[0];

    XLIST[i] = i;
    XLIST_ERR[i] = 0.;
    
    //cout << "MASS_BAND = " << MASS_BAND[i] << " +/- " << MASS_ERR_BAND[i] << "\n";
	  
  }

  //
  TGraph *gf_label = new TGraph(exp_nb, LABEL_INDX, LABEL_GF);
  gf_label -> SetName("gf_label");
  
  TGraphErrors *gf_exp1 = new TGraphErrors(1, XLIST_EXP1, MASS_EXP1, XLIST_EXP1_ERR, MASS_EXP1_ERR);
  gf_exp1 -> SetName("gf_exp1");
  gf_exp1 -> SetMarkerStyle(20);
  gf_exp1 -> SetMarkerSize(1.1);

  TGraphErrors *gf_exp2 = new TGraphErrors(1, XLIST_EXP2, MASS_EXP2, XLIST_EXP2_ERR, MASS_EXP2_ERR);
  gf_exp2 -> SetName("gf_exp2");
  gf_exp2 -> SetMarkerStyle(22);
  gf_exp2 -> SetMarkerSize(1.1);

  TGraphErrors *gf_exp3 = new TGraphErrors(1, XLIST_EXP3, MASS_EXP3, XLIST_EXP3_ERR, MASS_EXP3_ERR);
  gf_exp3 -> SetName("gf_exp3");
  gf_exp3 -> SetMarkerStyle(21);
  gf_exp3 -> SetMarkerSize(1.1);
  gf_exp3 -> SetMarkerColor(kRed);

  TGraphErrors *gf_exp4 = new TGraphErrors(1, XLIST_EXP4, MASS_EXP4, XLIST_EXP4_ERR, MASS_EXP4_ERR);
  gf_exp4 -> SetName("gf_exp4");
  gf_exp4 -> SetMarkerStyle(23);
  gf_exp4 -> SetMarkerSize(1.1);

  TGraphErrors *gf_exp5 = new TGraphErrors(1, XLIST_EXP5, MASS_EXP5, XLIST_EXP5_ERR, MASS_EXP5_ERR);
  gf_exp5 -> SetName("gf_exp5");
  gf_exp5 -> SetMarkerStyle(24);
  gf_exp5 -> SetMarkerSize(1.1);

  TGraphErrors *gf_exp6 = new TGraphErrors(1, XLIST_EXP6, MASS_EXP6, XLIST_EXP6_ERR, MASS_EXP6_ERR);
  gf_exp6 -> SetName("gf_exp6");
  gf_exp6 -> SetMarkerStyle(21);
  gf_exp6 -> SetMarkerSize(1.1);
  gf_exp6 -> SetMarkerColor(kBlue);
 
  TGraphErrors *gf_band = new TGraphErrors(nb_point, XLIST, MASS_BAND, XLIST_ERR, MASS_ERR_BAND);
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
  const double band_limit = 12. * MASS_ERR_BAND[0];
  //gf_exp1 -> GetXaxis() -> SetLabelOffset(999);
  //gf_exp1 -> GetXaxis() -> SetLabelSize(0);
  //gf_exp1 -> GetXaxis() -> SetTickLength(0);

  //gf_label -> RemovePoint(1);
  //gf_label -> RemovePoint(2);
  //gf_label -> RemovePoint(3);
  //gf_label -> RemovePoint(4);
  //gf_label -> RemovePoint(5);
  //gf_label -> RemovePoint(6);

  gf_label -> GetXaxis() -> SetLabelSize(0.06);
  
  gf_label -> GetYaxis() -> SetRangeUser(MASS_BAND[0] - band_limit, MASS_BAND[0] + band_limit);
  gf_label -> GetYaxis() -> SetTitleOffset(1.2);
  gf_label -> GetYaxis() -> SetTitle("M_{#omega} [MeV/c^{2}]");
  gf_label -> GetYaxis() -> CenterTitle();
  gf_label -> GetYaxis() -> SetLabelSize(0.04);
  gf_label -> GetYaxis() -> SetTitleSize(0.055);

  int bin_indx = 0;
  for (int i = 0; i < exp_nb; i ++) {
    bin_indx = gf_label -> GetXaxis() -> FindBin(LABEL_INDX[i]);
    cout << bin_indx << ", EXP_STR = " << EXP_STR[i] << endl;
    gf_label -> GetXaxis() -> SetBinLabel(bin_indx, EXP_STR[i]);
    //gf_label -> GetXaxis() -> SetBinLabel(18, "wfew");
  }
  
    
  

  TCanvas *cv = new TCanvas("cv", "", 0, 0, 1000, 700);
  //cv -> SetGrid();
  cv -> SetLeftMargin(0.2);
  cv -> SetBottomMargin(0.2);
  
  char display1[50];

  TPaveText *pt1 = new TPaveText(0.3, 0.84, 0.85, 0.86, "NDC");
  PteAttr(pt1);
  pt1 -> SetTextSize(0.05);
  pt1 -> SetTextColor(kRed);
  pt1 -> AddText(Form("M_{#omega}=%0.5g#pm%0.1g MeV/c^{2}", kloe_value, kloe_err)); 

  gf_label -> Draw("AP");
  gf_exp1 -> Draw("P");
  gf_exp6 -> Draw("P");
  gf_exp2 -> Draw("P");
  gf_exp3 -> Draw("P");
  gf_exp4 -> Draw("P");
  gf_exp5 -> Draw("P");
  
  mg -> Add(gf_band);
  mg -> Draw("C3");

  pt1 -> Draw("Same");

  TLegend * legd_gf = new TLegend(0.45, 0.7, .75, 0.8);

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

  legtextsize(legd_gf, 0.05);

  gPad -> Update();

  // save
  cv -> SaveAs("plots/mass_compr.pdf");
  
  return 0;
  
}
