#include "../header/plot.h"

const TString output_folder = "../../bkg_compr_IM3pi_7C";
const TString var_nm = "IM3pi_7C";
const TString unit = "[MeV]";
const TString var_symb = "M_{3#pi}";

const int binsize = 150;
const double var_min = 520;
const double var_max = 580;

const double sfw1d_isr3pi = 4.60022e-02;

//  f(x) = p0*exp(-0.5*((x-p1)/p2)^2)

// Linear background function
Double_t background(Double_t *x, Double_t *par) {
  return par[0] + par[1]*x[0];
}

// Gaussian peak function with 3 parameters
double gaussianPeak(double *x, double *par) {
  double r1 = (x[0] - par[1]) / par[2];
  double fitval = par[0] * TMath::Exp(-0.5 * r1 * r1);
  //cout << r1 << endl;
  return fitval;
}
// Sum of background and peak function
Double_t fitFunction(Double_t *x, Double_t *par) {
  double *p1 = &par[0];
  double *p2 = &par[2];
  return gaussianPeak(x, p1);
  //return background(x, p1) + gaussianPeak(x, p2);
}

int etapeak() {

  //gROOT->SetBatch(kTRUE);  
  gErrorIgnoreLevel = kError;
  TGaxis::SetMaxDigits(4);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  //gROOT->SetStyle("Plain");
  //gROOT->ForceStyle();
  
  TFile* intree = new TFile(output_folder + "/hist_" + var_nm + ".root");
  cout << intree -> GetName() << endl;
  
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

  TH1D* hist_isr3pi_sc = (TH1D *)intree -> Get("hist_isr3pi_sc");
  TH1D* hist_bkgsum_sc = (TH1D *)intree -> Get("hist_bkgsum_sc");
  TH1D* hist_etagam_sc = (TH1D *)intree -> Get("hist_etagam_sc");
  TH1D* hist_ksl_sc = (TH1D *)intree -> Get("hist_ksl_sc");
  TH1D* hist_omegapi_sc = (TH1D *)intree -> Get("hist_omegapi_sc");
  TH1D* hist_mcrest_sc = (TH1D *)intree -> Get("hist_mcrest_sc");
  TH1D* hist_eeg_sc = (TH1D *)intree -> Get("hist_eeg_sc");
  TH1D* hist_data = (TH1D *)intree -> Get("hist_data");
  
  // mcsum
  hist_isr3pi_sc -> Scale(sfw1d_isr3pi);
  //hist_isr3pi_sc -> Draw("Histo");

  TH1D* hist_mcsum_sc = (TH1D*) hist_bkgsum_sc -> Clone();
  hist_mcsum_sc -> Add(hist_isr3pi_sc, 1.);
  hist_mcsum_sc -> SetName("hist_mcsum_sc");
  format_h(hist_mcsum_sc, 2, 2);

  // plot distri.
  const double binwidth = getbinwidth(hist_data);
  const double ymax = hist_data -> GetMaximum();
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
  hist_data -> SetMarkerStyle(20);
  hist_data -> SetMarkerSize(0.7);

  //hist_data -> GetYaxis() -> SetNdivisions(505);
  //hist_data -> GetYaxis() -> SetLabelSize(20);
  //hist_data -> GetYaxis() -> SetTitleFont(43);
  //hist_data -> GetYaxis() -> SetLabelFont(43); // Absolute font size in pixel (precision 3)

  TPaveText *pt34 = new TPaveText(0.2, 0.8, 0.4, 0.85, "NDC");
  PteAttr(pt34); pt34 -> SetTextSize(0.06);
  //pt34 -> AddText("(a)");
  //pt34 -> AddText("(b)");
  //pt34 -> AddText("(c)");
  //pt34 -> AddText("(d)");

  // Fit data
  double xmin = hist_data -> GetXaxis() -> GetXmin();
  double xmax = hist_data -> GetXaxis() -> GetXmax();
  const int npar = 3;
  double rms = hist_data -> GetRMS();
  double mean = hist_data -> GetMean();
  double peak = hist_data -> GetMaximum();
  //double bin_width = getbinwidth(hist_data);
  //double fitpara[5] = {1, -1, peak, mean, rms};
  double fitpara[3] = {peak, mean, rms};
  
  //TF1 *fitFcn = new TF1("fitFcn", "[0] + [1]*x + [2]*exp(-0.5*((x-[3])/[4])^2)", 540, 560);
  const double fit_min = 545.;
  const double fit_max = 551.;
  TF1 *fitFcn = new TF1("fitFcn", fitFunction, fit_min, fit_max, 3);
  fitFcn -> SetParameters(fitpara);
  fitFcn -> SetLineColor(kRed);
  fitFcn -> SetLineWidth(1);
  fitFcn -> SetNpx(5000);
    
  TFitResultPtr r = hist_data -> Fit("fitFcn", "LM0", "", fit_min, fit_max);
  TF1 *fit_fun = hist_data -> GetFunction("fitFcn");
  fit_fun -> SetLineColor(kRed);
  fit_fun -> SetLineStyle(1);
  fit_fun -> SetLineWidth(1);
  fit_fun -> SetNpx(5000);

  // Calculate significance
  const double etamass_pdg = 547.862;
  const double etamass_pdg_err = 0.017;

  const double etamass_thesis = fit_fun -> GetParameter(1);
  const double etamass_thesis_err = fit_fun -> GetParError(1);
  double etamass_diff_thesis = TMath::Abs(etamass_thesis - etamass_pdg) * 1e3;
  
  const double etamass_kloe2007 = 547.873;
  const double etamass_kloe2007_err = TMath::Sqrt(0.007 * 0.007 + 0.031 *0.031);
  double etamass_diff_kloe2007 = TMath::Abs(etamass_kloe2007 - etamass_pdg) * 1e3;
  
  // https://arxiv.org/pdf/0707.4616
  
  cout << "eta mass pdg: " << etamass_pdg << "+/-" << etamass_pdg_err << "\n"
       << "thesis: " << etamass_thesis << "+/-" << etamass_thesis_err << ", mass - pdg: " << etamass_diff_thesis << " keV/c^2\n"
       << "kloe2007: " << etamass_kloe2007 << "+/-" << etamass_kloe2007_err << ", mass - pdg: " << etamass_diff_kloe2007 << " keV/c^2\n";
  
  // Draw
  TCanvas *cv = new TCanvas("cv", " ", 800, 700);
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

  //hist_data -> GetXaxis() -> SetTitle("M_{3#pi} [MeV/c^{2}]");
  hist_data -> GetXaxis() -> SetNdivisions(505);
  hist_data -> GetXaxis() -> SetTitle(var_symb + " " + unit);
  hist_data -> GetXaxis() -> CenterTitle();
  hist_data -> GetXaxis() -> SetTitleSize(0.06);
  hist_data -> GetXaxis() -> SetTitleOffset(1.0);
  hist_data -> GetXaxis() -> SetLabelOffset(0.01);
  hist_data -> GetXaxis() -> SetLabelSize(0.05);//0.03
  //hist_data -> GetXaxis() -> SetRangeUser(0., 40.); //chi2
  //hist_data -> GetXaxis() -> SetRangeUser(20., 140.); //pi0angle
  //hist_data -> GetXaxis() -> SetRangeUser(0.5, 1.);
  //hist_data -> GetXaxis() -> SetRangeUser(650., 950.); // 3pi mass omega region
  
  //hist_data -> GetYaxis() -> SetLabelFont(43); // Absolute front size in pixel (precision 3)
  //hist_data -> GetYaxis() -> SetLabelSize(0.05);
  hist_data -> GetYaxis() -> SetRangeUser(0.01, ymax * 1.2); 
  //hist_data -> GetYaxis() -> SetRangeUser(0.01, ymax * 1e2); // 3pi mass full region
  //hist_data -> GetYaxis() -> SetRangeUser(0., 5e4); // 3pi mass omega region  
  //hist_data -> GetYaxis() -> SetTitle(TString::Format("Entries/%0.2f", binwidth) + " " + unit);
  hist_data -> GetYaxis() -> SetTitle("Events");
  hist_data -> GetYaxis() -> CenterTitle();
  hist_data -> GetYaxis() -> SetTitleSize(0.06);
  hist_data -> GetYaxis() -> SetTitleOffset(1.2);
  
  //hist_data -> SetStats(0);

  format_h(hist_mcsum_sc, 1, 2);

  hist_isr3pi_sc -> SetLineStyle(2); // dashed line
  hist_etagam_sc -> SetLineStyle(3); // dotted line
  hist_ksl_sc -> SetLineStyle(4); // short-dotted-dashed line
  hist_omegapi_sc -> SetLineStyle(5); // dotted-dashed line
  hist_mcrest_sc -> SetLineStyle(6); // triple-dotted-dashed line
  hist_eeg_sc -> SetLineStyle(7); // long-dashed line
  
  hist_data -> Draw();
  fit_fun -> Draw("Same");

  /*
  hist_mcsum_sc -> Draw("SameHist");
  hist_isr3pi_sc -> Draw("SameHist");
  hist_omegapi_sc -> Draw("SameHist");
  hist_etagam_sc -> Draw("SameHist");
  hist_ksl_sc -> Draw("SameHist");
  hist_mcrest_sc -> Draw("SameHist");
  hist_eeg_sc -> Draw("SameHist");
  */
  
  pt34 -> AddText(Form("M_{#eta}=%0.3f#pm%0.3f [MeV/c^{2}]", fit_fun -> GetParameter(1), fit_fun -> GetParError(1)));
  pt34 -> Draw("Same");
  
  //gPad->SetLogy();

  TLegend *legd_cv = new TLegend(0.75, 0.35, 0.9, 0.9);
  
  legd_cv -> SetTextFont(132);
  legd_cv -> SetFillStyle(0);
  legd_cv -> SetBorderSize(0);
  legd_cv -> SetNColumns(1);
  
  legd_cv -> AddEntry(hist_data, "Data", "lep");
  legd_cv -> AddEntry(hist_mcsum_sc, "MC sum", "l");
  legd_cv -> AddEntry(hist_isr3pi_sc, "#pi^{+}#pi^{-}#pi^{0}#gamma", "l");
  legd_cv -> AddEntry(hist_omegapi_sc, "#omega#pi^{0}", "l");
  //legd_cv -> AddEntry(hist_etagam_sc, "#varphi#rightarrow#eta#gamma#rightarrow#pi^{+}#pi^{-}#pi^{0}#gamma", "l");
  legd_cv -> AddEntry(hist_etagam_sc, "#eta#gamma", "l");
  legd_cv -> AddEntry(hist_ksl_sc, "K_{L}K_{S}", "l");
  legd_cv -> AddEntry(hist_eeg_sc, "e^{+}e^{-}#gamma", "l");
  legd_cv -> AddEntry(hist_mcrest_sc, "Others", "l");
  
  //legd_cv -> Draw("Same");
  
  legtextsize(legd_cv, 0.04);

  cv -> SaveAs("etamass.pdf");
 
  return 0;
  
}
