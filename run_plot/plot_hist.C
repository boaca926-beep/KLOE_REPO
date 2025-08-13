#include "../header/bkg_compr.h"
//#include "../hist.h"
//#include "../fitfun.h"
#include "../header/plot.h"

void SetCVAttr(TH1D *h1d) {

  const double ymax_normed = h1d -> GetMaximum();
  cout << "ymax = " << ymax_normed << endl;
  
  h1d -> GetXaxis() -> SetTitle(var_symb + " " + unit);
  h1d -> GetXaxis() -> CenterTitle();
  h1d -> GetXaxis() -> SetTitleSize(0.07);
  h1d -> GetXaxis() -> SetTitleOffset(.9);
  h1d -> GetXaxis() -> SetLabelOffset(.01);
  h1d -> GetXaxis() -> SetLabelSize(0.05);//0.03
  //h1d -> GetXaxis() -> SetRangeUser(110., 450.);
  //h1d -> GetXaxis() -> SetRangeUser(-500., 50.);
  //h1d -> GetXaxis() -> SetRangeUser(20., 180.);
  
  //h1d -> GetYaxis() -> SetLabelFont(43); // Absolute front size in pixel (precision 3)
  h1d -> GetYaxis() -> SetLabelSize(0.05);
  //h1d -> GetYaxis() -> SetRangeUser(1e-8, ymax_normed * 60.); //*1e2, *1.3, *30 *60
  h1d -> GetYaxis() -> SetRangeUser(0., ymax_normed * 1.5); //*1e2, *1.3, *30 *60
  
  //h1d -> GetYaxis() -> SetTitle(TString::Format("Entries/%0.2f", binwidth) + " " + unit);
  //h1d -> GetYaxis() -> SetTitle("Events");
  h1d -> GetYaxis() -> CenterTitle();
  h1d -> GetYaxis() -> SetTitleSize(0.07);
  h1d -> GetYaxis() -> SetTitleOffset(1.1);
  
  //hist_data -> SetStats(0);
  
}

int plot_hist(){

  //gROOT->SetBatch(kTRUE);  
  gErrorIgnoreLevel = kError;
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  //TGaxis::SetMaxDigits(2);
  //gROOT->SetStyle("Plain");
  gROOT->ForceStyle();
  
  TFile* intree = new TFile(output_folder + "/hist_" + var_nm + ".root");
  
  TIter next_tree(intree -> GetListOfKeys());

  TString objnm_tree, classnm_tree;

  int i = 0;
  TKey *key;
  
  while ( (key = (TKey *) next_tree() ) ) {
    
    i ++;
    
    objnm_tree   =  key -> GetName();
    classnm_tree = key -> GetClassName();
    key -> GetSeekKey();
    
    //cout << "tree" << i << ": classnm = " << classnm_tree << ", objnm = " << objnm_tree << endl;
    
  }

  TH1D* hist_data = (TH1D *)intree -> Get("hist_data");
  TH1D* hist_rhopi = (TH1D *)intree -> Get("hist_rhopi");
  TH1D* hist_ksl = (TH1D *)intree -> Get("hist_ksl");
  TH1D* hist_kpm = (TH1D *)intree -> Get("hist_kpm");
  
  TH1D* hist_bkgsum_sc = (TH1D *)intree -> Get("hist_bkgsum_sc");
  TH1D* hist_isr3pi_sc = (TH1D *)intree -> Get("hist_isr3pi_sc");
  TH1D* hist_etagam_sc = (TH1D *)intree -> Get("hist_etagam_sc");
  
  // mcsum 
  TH1D* hist_mcsum_sc = (TH1D*) hist_bkgsum_sc -> Clone();
  hist_mcsum_sc -> Add(hist_isr3pi_sc, 1.);
  hist_mcsum_sc -> SetName("hist_mcsum_sc");
  format_h(hist_mcsum_sc, 1, 2);

  // plot distri.
  const double binwidth = getbinwidth(hist_data);
  double ymax = hist_data -> GetMaximum();
  const double residul_min = -25., residul_max = 25.;
  
  TH1D * hresidul = new TH1D("hresidul", "", binsize, var_min, var_max);
  TH1D * hresidul_distr = new TH1D("hresidul_distr", "", 200, residul_min, residul_max);

  double nb_data = 0., evnt_err = 0., Nb_data = 0., nb_mcsum = 0.;
  double residul = 0.;
  
  for (int j = 1; j <= binsize; j ++ ) {

    nb_data =  hist_data -> GetBinContent(j);
    Nb_data += nb_data;

    nb_mcsum = hist_mcsum_sc -> GetBinContent(j);  
    evnt_err = TMath::Sqrt(nb_data + nb_mcsum);
    
    // residual
    residul = (nb_data - nb_mcsum) / evnt_err;

    if (nb_data > 0. && nb_mcsum > 0.) {
      
     hresidul -> SetBinContent(j, residul);
     hresidul_distr -> Fill(residul);

    }

    //cout << j << ": nb_data = " << nb_data << " +/- " << evnt_err << ", mcsum = " << nb_mcsum << ", residul= " << residul << "\n";

    
  }

  //
  hist_data -> SetMarkerStyle(21);
  hist_data -> SetMarkerSize(0.7);

  //hist_data -> GetYaxis() -> SetNdivisions(505);
  //hist_data -> GetYaxis() -> SetLabelSize(20);
  //hist_data -> GetYaxis() -> SetTitleFont(43);
  //hist_data -> GetYaxis() -> SetLabelFont(43); // Absolute font size in pixel (precision 3)

  //hist_etagam_sc -> Draw();
  
  
  TCanvas *cv = new TCanvas("cv_" + var_nm, " ", 700, 700);
  cv -> SetBottomMargin(0.15);//0.007
  cv -> SetLeftMargin(0.15);

  //TPad *p2 = new TPad("p2", "p2", 0., 0., 1., 0.25);
  //p2 -> Draw();
  //p2 -> SetBottomMargin(0.3);
  //p2 -> SetGrid();
  
  //TPad *p1 = new TPad("p1", "p1", 0., 0.24, 1., 1.);
  //p1 -> Draw();
  //p1 -> SetBottomMargin(0.03);//0.007
  //p1 -> cd();

  // Normalization
  
  //TH1D* hbkg_normed = (TH1D*) hist_eeg_sc -> Clone();
  //TH1D* hbkg_normed = (TH1D*) hist_ksl_sc -> Clone();
  TH1D* hbkg_normed = (TH1D*) hist_etagam_sc -> Clone();
  //TH1D* hbkg_normed = (TH1D*) hist_omegapi_sc -> Clone();
  //TH1D* hbkg_normed = (TH1D*) hist_rhopi -> Clone();

  TH1D* hsig_normed = (TH1D*) hist_isr3pi_sc -> Clone();
   
  const double bin_size = hist_data -> GetNbinsX();
  //cout << bin_size << endl;
  double normbkg_tmp = hbkg_normed -> Integral(1, bin_size);
  double normsig_tmp = hsig_normed -> Integral(1, bin_size);
  hbkg_normed -> Scale(1. / normbkg_tmp);
  hsig_normed -> Scale(1. / normsig_tmp);
  cout << hist_kpm -> Integral(1, bin_size) << endl;
  cout << hist_ksl -> Integral(1, bin_size) << endl;
   
  //formatfill_h(hbkg_normed, 2, 4050);
  formatfill_h(hsig_normed, 13, 3001);
  format_h(hbkg_normed, 2, 2);

  //
  const double pvalue=TMath::Prob(38., 7);

  cout << "pvalue = " << pvalue << endl;
  
  // arrow
  const double cut_tmp = 320;
  //TArrow arrow(cut_tmp, 8e4, cut_tmp, 2e-3, 0.02,"|>");
  TArrow arrow(cut_tmp, 0.03, cut_tmp, 2e-3, 0.02,"|>");

  arrow.SetLineWidth(4);
  arrow.SetLineColor(1);
  arrow.SetFillStyle(1001);
  arrow.SetFillColor(kBlack);

  SetCVAttr(hbkg_normed);
  
  format_h(hist_mcsum_sc, 1, 2);

  hbkg_normed -> Draw("Hist");
  hsig_normed -> Draw("SameHist");
  arrow.DrawClone();
  //hist_mcsum_sc -> Draw("SameHist");
  //hist_isr3pi_sc -> Draw("SameHist");
  //hist_omegapi_sc -> Draw("SameHist");
  //hist_etagam_sc -> Draw("SameHist");
  //hist_ksl_sc -> Draw("SameHist");
  //hist_rhopi -> Draw("SameHist");
  //hist_mcrest_sc -> Draw("SameHist");
  //hist_eeg_sc -> Draw("SameHist");
  
  //gPad->SetLogy();

  //TLegend *legd_cv = new TLegend(0.25, 0.79, 0.8, 0.87);
  TLegend *legd_cv = new TLegend(0.2, 0.7, 0.6, 0.9);
  
  legd_cv -> SetTextFont(132);
  legd_cv -> SetFillStyle(0);
  legd_cv -> SetBorderSize(0);
  //legd_cv -> SetNColumns(2);

  //const TString bkg_str = "e^{+}e^{-}#gamma";
  //const TString bkg_str = "#omega#pi^{0}";
  //const TString bkg_str = "K_{L}K_{S}";
  //const TString bkg_str = "#eta#gamma";
  //const TString bkg_str = "#phi#rightarrow#rho#pi#rightarrow#pi^{+}#pi^{-}#pi^{0}#gamma";
  const TString bkg_str = "#phi#rightarrow#eta#gamma";
  //const TString bkg_str = "#rho#pi";
  
  
  //"#rho#pi#rightarrow#pi^{+}#pi^{-}#pi^{0}#gamma", "#eta#gamma", "#omega#pi^{0}", "e^{+}e^{-}#gamma", "K_{L}K_{S}"
  legd_cv -> AddEntry(hsig_normed, "e^{+}e^{-}#rightarrow#pi^{+}#pi^{-}#pi^{0}#gamma", "f");
  //legd_cv -> AddEntry(hsig_normed, "#pi^{+}#pi^{-}#pi^{0}#gamma", "f");
  legd_cv -> AddEntry(hbkg_normed, bkg_str, "l");
  
  //legd_cv -> AddEntry(hist_data, "Data", "lep");
  //legd_cv -> AddEntry(hist_mcsum_sc, "MC sum", "l");
  //legd_cv -> AddEntry(hist_isr3pi_sc, "#pi^{+}#pi^{-}#pi^{0}#gamma", "l");
  //legd_cv -> AddEntry(hist_omegapi_sc, "#varphi#rightarrow#omega#pi^{0}", "l");
  //legd_cv -> AddEntry(hist_ksl_sc, "#varphi#rightarrowK_{L}K_{S}", "l");
  //legd_cv -> AddEntry(hist_ksl_sc, "K_{L}K_{S}", "l");
  //legd_cv -> AddEntry(hist_etagam_sc, "#eta#gamma", "l");
  //legd_cv -> AddEntry(hist_rhopi, "#rho#pi", "l");
  
  //legd_cv -> AddEntry(hist_etagam_sc, "#eta#gamma#rightarrow#pi^{+}#pi^{-}#pi^{0}#gamma", "l");
  //legd_cv -> AddEntry(hist_eeg_sc, "e^{+}e^{-}#gamma", "l");
  //legd_cv -> AddEntry(hist_rhopi, "#rho#pi#rightarrow#pi^{+}#pi^{-}#pi^{0}#gamma", "l");
  //legd_cv -> AddEntry(hist_mcrest_sc, "MC rest", "l"); 
  
  legd_cv -> Draw("Same");
  
  legtextsize(legd_cv, 0.07);
  // save
  
  //cv_residul -> SaveAs("./output/cv_residul_" + var_nm + ".pdf");
  cv -> SaveAs(output_folder + "/cv_" + var_nm + ".pdf");
  //cv -> Write();

 
  return 0;
  
}

