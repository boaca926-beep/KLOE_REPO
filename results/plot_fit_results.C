#include "compr.h"
#include "../plot.h"

int plot_fit_results() {

  gErrorIgnoreLevel = kError;
  TGaxis::SetMaxDigits(4);
  //gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetFitFormat("6.4g");
  
  cout << "Plotting fit results BW v.s. VMD" << endl;

  TFile *f_input = new TFile("result_norm.root");

  TIter next_tree(f_input -> GetListOfKeys());
  
  TString objnm_tree, classnm_tree;
  
  int i = 0;
  
  while ( (key = (TKey *) next_tree() ) ) {// start tree while lop
    
    i ++;
    
    objnm_tree   =  key -> GetName();
    classnm_tree = key -> GetClassName();
    key -> GetSeekKey();
    
    cout << " tree/histo" << i << ": classnm = " << classnm_tree << ", objnm = " << objnm_tree << endl;
    
  }

  TGraphErrors * gf_data_KLOE2 = (TGraphErrors *) n3pi_obs -> Clone();
  TGraphErrors * gf_fit_bw = (TGraphErrors *) n3pi_fit -> Clone();
  TGraphErrors * gf_n3pi_diff = (TGraphErrors *) n3pi_diff -> Clone();
  gf_n3pi_diff -> Draw();
  
  double Mass_bw = 0., Mass_err_bw = 0.;
  double Gamma_bw = 0., Gamma_err_bw = 0.;
  double BB_bw = 0., BB_err_bw = 0.;
  double chi2_sum_bw = 0.;
  double Delta_m3pi = 0.;
  int ndf_bw = 0;
  
  TTree* TRESULT = (TTree*)f_input -> Get("TRESULT");

  for (Int_t irow = 0; irow < TRESULT -> GetEntries(); irow ++) {// loop trees

    TRESULT -> GetEntry(irow);

    Mass_bw = TRESULT -> GetLeaf("Br_Mass") -> GetValue(0);
    Mass_err_bw = TRESULT -> GetLeaf("Br_Mass") -> GetValue(1);
    
    Gamma_bw = TRESULT -> GetLeaf("Br_Gamma") -> GetValue(0);
    Gamma_err_bw = TRESULT -> GetLeaf("Br_Gamma") -> GetValue(1);
    
    BB_bw = TRESULT -> GetLeaf("Br_BB") -> GetValue(0);
    BB_err_bw = TRESULT -> GetLeaf("Br_BB") -> GetValue(1);

    chi2_sum_bw = TRESULT -> GetLeaf("Br_chi2_sum") -> GetValue(0);
    ndf_bw = TRESULT -> GetLeaf("Br_ndf") -> GetValue(0);

    Delta_m3pi = TRESULT -> GetLeaf("Br_Delta_m3pi") -> GetValue(0);
    

    cout << "Mass_bw = " << Mass_bw << " +/- " << Mass_err_bw << "\n"
	 << "Gamma_bw = " << Gamma_bw << " +/- " << Gamma_err_bw << "\n"
	 << "BB_bw = " << BB_bw << " +/- " << BB_err_bw << "\n"
	 << "chi2_sum_bw = " << chi2_sum_bw << "\n"
	 << "ndf_bw = " << ndf_bw << "\n";

  }

  f_input -> Close();

  TFile *f_input = new TFile("result_vmd.root");

  TGraphErrors * gf_fit_vmd = (TGraphErrors *) n3pi_fit -> Clone();
  TGraphErrors * gf_n3pi_vmd_diff = (TGraphErrors *) n3pi_diff -> Clone();
  //gf_n3pi_vmd_diff -> Draw();
  //gf_n3pi_diff -> Draw();
  
  double Mass_vmd = 0., Mass_err_vmd = 0.;
  double Gamma_vmd = 0., Gamma_err_vmd = 0.;
  double BB_vmd = 0., BB_err_vmd = 0.;
  double chi2_sum_vmd = 0.;
  int ndf_vmd = 0;
  int fit_indx = 0;

  const int list_size = 1000;
  double M3PI_FIT[list_size], M3PI_FIT_ERR[list_size];
  
  TTree* TRESULT = (TTree*)f_input -> Get("TRESULT");

  for (Int_t irow = 0; irow < TRESULT -> GetEntries(); irow ++) {// loop trees

    TRESULT -> GetEntry(irow);

    Mass_vmd = TRESULT -> GetLeaf("Br_Mass") -> GetValue(0);
    Mass_err_vmd = TRESULT -> GetLeaf("Br_Mass") -> GetValue(1);
    
    Gamma_vmd = TRESULT -> GetLeaf("Br_Gamma") -> GetValue(0);
    Gamma_err_vmd = TRESULT -> GetLeaf("Br_Gamma") -> GetValue(1);
    
    BB_vmd = TRESULT -> GetLeaf("Br_BB") -> GetValue(0);
    BB_err_vmd = TRESULT -> GetLeaf("Br_BB") -> GetValue(1);

    chi2_sum_vmd = TRESULT -> GetLeaf("Br_chi2_sum") -> GetValue(0);
    ndf_vmd = TRESULT -> GetLeaf("Br_ndf") -> GetValue(0);

    fit_indx ++;
    
    cout << "Mass_vmd = " << Mass_vmd << " +/- " << Mass_err_vmd << "\n"
	 << "Gamma_vmd = " << Gamma_vmd << " +/- " << Gamma_err_vmd << "\n"
	 << "BB_vmd = " << BB_vmd << " +/- " << BB_err_vmd << "\n"
	 << "chi2_sum_vmd = " << chi2_sum_vmd << "\n"
	 << "ndf_vmd = " << ndf_vmd << "\n";

  }

  
  // Plot
  TCanvas *cv_n3pi_fit = new TCanvas("cv_n3pi_fit", "Nobs 3pi", 0, 0, 1000, 800);

  char display1[50], display2[50], display3[50], display4[50], display5[50];

  TPaveText *pt1 = new TPaveText(0.12, 0.70, 0.43, 0.72, "NDC");
  TPaveText *pt2 = new TPaveText(0.12, 0.64, 0.43, 0.66, "NDC");
  TPaveText *pt3 = new TPaveText(0.12, 0.58, 0.43, 0.60, "NDC");
  TPaveText *pt4 = new TPaveText(0.12, 0.52, 0.43, 0.54, "NDC");

  TPaveText *pt11 = new TPaveText(0.64, 0.70, 0.80, 0.72, "NDC");
  TPaveText *pt22 = new TPaveText(0.64, 0.64, 0.80, 0.66, "NDC");
  TPaveText *pt33 = new TPaveText(0.64, 0.58, 0.80, 0.60, "NDC");
  TPaveText *pt44 = new TPaveText(0.64, 0.52, 0.80, 0.54, "NDC");
  
  PteAttr(pt1); pt1 -> SetTextSize(0.04); pt1 -> SetTextColor(kBlack);
  PteAttr(pt2); pt2 -> SetTextSize(0.04); pt2 -> SetTextColor(kBlack);
  PteAttr(pt3); pt3 -> SetTextSize(0.04); pt3 -> SetTextColor(kBlack);
  PteAttr(pt4); pt4 -> SetTextSize(0.04); pt4 -> SetTextColor(kBlack);

  PteAttr(pt11); pt11 -> SetTextSize(0.04); pt11 -> SetTextColor(kBlue);
  PteAttr(pt22); pt22 -> SetTextSize(0.04); pt22 -> SetTextColor(kBlue);
  PteAttr(pt33); pt33 -> SetTextSize(0.04); pt33 -> SetTextColor(kBlue);
  PteAttr(pt44); pt44 -> SetTextSize(0.04); pt44 -> SetTextColor(kBlue);
  
  //sprintf(display1,"#chi^{2}/ndf=%0.2f", chi2_sum / ndf);
  pt1 -> AddText(Form("#Gamma_{#omega}=%0.2f#pm%0.2f MeV", Gamma_bw, Gamma_err_bw)); 
  pt2 -> AddText(Form("B_{ee}B_{3#pi}=(%0.2f#pm%0.2f)#times10^{-5}", BB_bw * 1e5, BB_err_bw * 1e5));
  pt3 -> AddText(Form("M_{#omega}=%0.2f#pm%0.2f MeV/c^{2}", Mass_bw, Mass_err_bw)); 
  //pt4 -> AddText(Form("#chi^{2}/ndf=%0.2f", chi2_sum_bw / ndf_bw)); 
  pt4 -> AddText(Form("#chi^{2}/ndf=%0.1f/%0.2d", 13.4, ndf_bw)); 

  pt11 -> AddText(Form("#Gamma_{#omega}=%0.2f#pm%0.2f", Gamma_vmd, Gamma_err_vmd)); 
  pt22 -> AddText(Form("B_{ee}B_{3#pi}=%0.2f#pm%0.2f", BB_vmd * 1e5, BB_err_vmd * 1e5));
  pt33 -> AddText(Form("M_{#omega}=%0.2f#pm%0.2f", Mass_vmd, Mass_err_vmd)); 
  //pt44 -> AddText(Form("#chi^{2}/ndf=%0.2f", chi2_sum_vmd / ndf_vmd));  
  pt44 -> AddText(Form("#chi^{2}/ndf=%0.1f/%0.2d", chi2_sum_vmd, ndf_vmd)); 
 
  SetGFAttr(gf_data_KLOE2, "M_{3#pi} [MeV/c^{2}]", " ");
  gf_data_KLOE2 -> SetMarkerStyle(21);
  //gf_data_KLOE2 -> SetMarkerColor(kBlue);
  gf_data_KLOE2 -> SetLineColor(kBlack);
  gf_data_KLOE2 -> SetLineWidth(2);

  gf_data_KLOE2 -> GetXaxis() -> SetTitleOffset(1.2);
  gf_data_KLOE2 -> GetXaxis() -> SetLabelOffset(1.2);
  gf_data_KLOE2 -> GetXaxis() -> SetLabelSize(0.03);
  gf_data_KLOE2 -> GetYaxis() -> SetTitle(TString::Format("Events/[%0.2f MeV/c^{2}]", Delta_m3pi));
  gf_data_KLOE2 -> GetYaxis() -> SetTitleSize(0.06);
  
  gf_fit_vmd -> SetLineStyle(kDashed);
  gf_fit_vmd -> SetLineColor(kBlue);
  
  gf_fit_bw -> SetLineColor(kBlack);

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
  
  gf_data_KLOE2 -> Draw("PA");
  gf_fit_bw -> Draw("cSame");
  gf_fit_vmd -> Draw("cSame");
  pt1 -> Draw("Same"), pt11 -> Draw("Same");
  pt2 -> Draw("Same"), pt22 -> Draw("Same");
  pt3 -> Draw("Same"), pt33 -> Draw("Same");
  pt4 -> Draw("Same"), pt44 -> Draw("Same");

  TLegend * legd_cv = new TLegend(0.15, 0.8, 0.6, 0.9);
  
  SetLegend(legd_cv);
  legd_cv -> SetNColumns(3);
  
  legd_cv -> AddEntry(gf_data_KLOE2, "Data", "lep");
  legd_cv -> AddEntry(gf_fit_bw, "BW", "l");
  legd_cv -> AddEntry(gf_fit_vmd, "VMD", "l");
  
  legd_cv -> Draw("Same");

  legtextsize(legd_cv, 0.06);
  
  p2 -> cd();
  
  gf_n3pi_diff -> SetLineColor(0);
  gf_n3pi_diff -> GetYaxis() -> SetNdivisions(505);
  gf_n3pi_diff -> GetYaxis() -> SetTitleSize(33);
  gf_n3pi_diff -> GetYaxis() -> SetTitleFont(43);
  gf_n3pi_diff -> GetYaxis() -> SetTitleOffset(1.2);
  gf_n3pi_diff -> GetYaxis() -> SetLabelFont(43); // Absolute front size in pixel (precision 3)
  gf_n3pi_diff -> GetYaxis() -> SetLabelSize(25);
  gf_n3pi_diff -> GetYaxis() -> CenterTitle();
  gf_n3pi_diff -> GetYaxis() -> SetTitle("Residual");
  //gf_n3pi_diff -> GetYaxis() -> SetRangeUser(diff_min * 3., diff_max * 3.);
  //gf_n3pi_diff -> GetYaxis() -> SetTitle("(N_{d} - N_{mc})/#sqrt{N_{d}}");
  //gf_n3pi_diff -> GetXaxis() -> SetLabelFont(43); // Absolute font size in pixel (precision 3)
  gf_n3pi_diff -> GetXaxis() -> SetTitleOffset(0.85);
  //gf_n3pi_diff -> GetXaxis() -> SetTitleSize(45);
  gf_n3pi_diff -> GetXaxis() -> SetLabelSize(0.15);
  //gf_n3pi_diff -> GetXaxis() -> SetLabelOffset(0.2);
  gf_n3pi_diff -> GetXaxis() -> SetTitleSize(0.2);
  gf_n3pi_diff -> GetXaxis() -> SetTitle("M_{3#pi} [MeV/c^{2}]");
  gf_n3pi_diff -> GetXaxis() -> CenterTitle();
  gf_n3pi_diff -> SetMarkerStyle(21);
  gf_n3pi_diff -> SetMarkerSize(0.7);
  //gf_n3pi_diff -> Draw("E0");

  //TMultiGraph *mg = new TMultiGraph();
  //mg -> SetTitle("Exclusion graphs");
  
  //mg -> Add(gf_n3pi_diff);
  gf_n3pi_diff -> SetFillColor(1);
  gf_n3pi_diff -> SetLineStyle(9);
  gf_n3pi_diff -> SetFillStyle(3001);

  gf_n3pi_vmd_diff -> SetFillColor(4);
  gf_n3pi_vmd_diff -> SetLineStyle(9);
  gf_n3pi_vmd_diff -> SetFillStyle(3001);
  
  gf_n3pi_diff -> Draw("a3");
  gf_n3pi_vmd_diff -> Draw("3");
  
  // save
  cv_n3pi_fit -> SaveAs("./plots/crx3pi_bw_vmd.pdf");
  
  
  return 0;
  
}
