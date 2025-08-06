#include "../hist.h"
#include "./plot.h"
#include "../fitfun.h"

int IM3pi() {
  
  gErrorIgnoreLevel = kError;
  TGaxis::SetMaxDigits(3);
  //gStyle->SetOptStat(0);

  cout << "Ploting" << endl;

  TFile *f_input = new TFile("../crx3pi/output8373_finebin/crx3pi0.root");

  TIter next_tree(f_input -> GetListOfKeys());

  TString objnm_tree, classnm_tree;
  
  int i = 0;

  while ( (key = (TKey *) next_tree() ) ) {// start tree while lop

    i ++;

    objnm_tree   =  key -> GetName();
    classnm_tree = key -> GetClassName();
    key -> GetSeekKey();

    //cout << " tree" << i << ": classnm = " << classnm_tree << ", objnm = " << objnm_tree << endl;
      
  }

  // MC sum
  TH1D* h1d_isr3pi = (TH1D*)f_input -> Get("h1d_IM3pi_TISR3PI_" + sig_type);
  
  TH1D* h1d_IM3pi_TMCSUM = (TH1D*) h1d_IM3pi_TEEG_Scaled -> Clone();
  h1d_IM3pi_TMCSUM -> Add(h1d_IM3pi_TOMEGAPI_Scaled, 1.);
  h1d_IM3pi_TMCSUM -> Add(h1d_IM3pi_TETAGAM_Scaled, 1.);
  h1d_IM3pi_TMCSUM -> Add(h1d_IM3pi_TKSL_Scaled, 1.);
  h1d_IM3pi_TMCSUM -> Add(h1d_isr3pi, 1.);
  h1d_IM3pi_TMCSUM -> Add(h1d_IM3pi_TMCREST_Scaled, 1.);
  h1d_IM3pi_TMCSUM -> SetName("h1d_IM3pi_TMCSUM");
  format_h(h1d_IM3pi_TMCSUM, 2, 2);

  // residual and disr.
  const double bins = h1d_IM3pi_TDATA -> GetNbinsX();
  const double xmin = h1d_IM3pi_TDATA -> GetXaxis() -> GetXmin();
  const double xmax = h1d_IM3pi_TDATA -> GetXaxis() -> GetXmax();

  cout << "bins = " << bins << endl;
  
  TH1D * hresidul = new TH1D("hresidul", "", bins, xmin, xmax);
  TH1D * hresidul_distr = new TH1D("hresidul_distr", "", 150, -50., 50.);

  double nb_data = 0., evnt_err = 0., nb_mcsum = 0., residul = 0.;

  for (int j = 1; j <= h1d_IM3pi_TDATA -> GetNbinsX(); j ++ ) {

    nb_data =  h1d_IM3pi_TDATA -> GetBinContent(j);
    evnt_err = TMath::Sqrt(nb_data);

    nb_mcsum = h1d_IM3pi_TMCSUM -> GetBinContent(j);  

    // residual
    residul = (nb_data - nb_mcsum) / evnt_err;

    //cout << j << ": nb_data = " << nb_data << " +/- " << evnt_err << ", mcsum = " << nb_mcsum << ", residul_sfw = " << residul_sfw << "\n";

    if (nb_data > 0. && nb_mcsum > 0.) {
      
     hresidul -> SetBinContent(j, residul);
     hresidul_distr -> Fill(residul);

    }

  }

  // fit residual distr.
  const int residul_npar = 3;
 
  double residul_fitpar[residul_npar] = {hresidul_distr -> GetMaximum(), hresidul_distr -> GetMean(), hresidul_distr -> GetRMS()};

  TF1 *fitfun = new TF1("fitfun", gauss1d, residul_min, residul_max, residul_npar);
  fitfun -> SetParNames("peak","mean","sigma");
  fitfun -> SetParameters(residul_fitpar);
  fitfun -> SetLineWidth(2);
  fitfun -> SetNpx(5000);

  TFitResultPtr r = hresidul_distr -> Fit(fitfun, "L0", " ", -10., 10.);

  residul_fitpar[0]  = fitfun -> GetParameter(0); // peak
  residul_fitpar[1]  = fitfun -> GetParameter(1); // mean
  residul_fitpar[2]  = fitfun -> GetParameter(2); // width

  // get gaussians 
  gfratio_fit = new TF1("gfratio_fit", gauss1d, residul_min, residul_max, residul_npar);
  
  gfratio_fit -> SetParameters(residul_fitpar);
  gfratio_fit -> SetLineColor(kRed); 
  gfratio_fit -> SetLineStyle(1);
  gfratio_fit -> SetNpx(5000);

  // plot
  TCanvas *cv = new TCanvas("cv", "", 0, 0, 700, 700);

  TPad *p2 = new TPad("p2", "p2", 0., 0., 1., 0.25);
  p2 -> Draw();
  p2 -> SetBottomMargin(0.3);
  p2 -> SetGrid();
  
  TPad *p1 = new TPad("p1", "p1", 0., 0.24, 1., 1.);
  p1 -> Draw();
  p1 -> SetBottomMargin(0.03);//0.007
  p1 -> cd();

  const double binwidth = getbinwidth(h1d_IM3pi_TDATA);
  const double ymax = h1d_IM3pi_TDATA -> GetMaximum();

  h1d_IM3pi_TDATA -> SetMarkerStyle(21);
  h1d_IM3pi_TDATA -> SetMarkerSize(0.7);

  //h1d_IM3pi_TDATA -> GetYaxis() -> SetNdivisions(505);
  //h1d_IM3pi_TDATA -> GetYaxis() -> SetLabelSize(20);
  //h1d_IM3pi_TDATA -> GetYaxis() -> SetTitleFont(43);
  //h1d_IM3pi_TDATA -> GetYaxis() -> SetLabelFont(43); // Absolute font size in pixel (precision 3)
  
  //h1d_IM3pi_TDATA -> GetYaxis() -> SetRangeUser(0.1, ymax * 1.2);
  h1d_IM3pi_TDATA -> GetYaxis() -> SetRangeUser(0., 3.5e4);
  h1d_IM3pi_TDATA -> GetYaxis() -> SetTitle(TString::Format("Entries/%0.2f", binwidth) + " " + unit);
  h1d_IM3pi_TDATA -> GetYaxis() -> CenterTitle();
  h1d_IM3pi_TDATA -> GetYaxis() -> SetTitleSize(0.04);
  h1d_IM3pi_TDATA -> GetYaxis() -> SetTitleOffset(1.2);
  
  h1d_IM3pi_TDATA -> GetXaxis() -> SetRangeUser(600., 950.); 
  h1d_IM3pi_TDATA -> GetXaxis() -> SetTitle(xtit + " " + unit);
  h1d_IM3pi_TDATA -> GetXaxis() -> CenterTitle();
  h1d_IM3pi_TDATA -> GetXaxis() -> SetTitleSize(0.04);
  h1d_IM3pi_TDATA -> GetXaxis() -> SetTitleOffset(1.1);
  h1d_IM3pi_TDATA -> GetXaxis() -> SetLabelSize(0.03);
  h1d_IM3pi_TDATA -> SetStats(0);
  
  h1d_IM3pi_TDATA -> Draw("E0");
  h1d_IM3pi_TEEG_Scaled -> Draw("SameHist");
  h1d_IM3pi_TOMEGAPI_Scaled -> Draw("SameHist");
  h1d_IM3pi_TETAGAM_Scaled -> Draw("SameHist");
  h1d_IM3pi_TKSL_Scaled -> Draw("SameHist");
  h1d_IM3pi_TISR3PI_SIG_Scaled -> Draw("SameHist");
  h1d_IM3pi_TMCREST_Scaled -> Draw("SameHist");
  h1d_IM3pi_TMCSUM -> Draw("SameHist");
    
  //gPad -> SetLogy();

  TLegend *legd_cv = new TLegend(0.6, 0.55, 1., 0.99);
  
  legd_cv -> SetTextFont(132);
  legd_cv -> SetFillStyle(0);
  legd_cv -> SetBorderSize(0);
  legd_cv -> SetNColumns(1);
  
  legd_cv -> AddEntry(h1d_IM3pi_TDATA,  "Data", "lep");
  legd_cv -> AddEntry(h1d_IM3pi_TMCSUM, "MC sum", "l");
  legd_cv -> AddEntry(h1d_IM3pi_TISR3PI_SIG_Scaled,  "#pi^{+}#pi^{-}#pi^{0}#gamma (phok5)", "l");
  legd_cv -> AddEntry(h1d_IM3pi_TOMEGAPI_Scaled,     "#phi#rightarrow#omega#pi^{0}", "l");
  legd_cv -> AddEntry(h1d_IM3pi_TETAGAM_Scaled,      "#phi#rightarrow#eta#gamma#rightarrow#pi^{+}#pi^{-}#pi^{0}#gamma", "l");
  legd_cv -> AddEntry(h1d_IM3pi_TKSL_Scaled,         "#phi#rightarrowK_{L}K_{S}", "l");
  legd_cv -> AddEntry(h1d_IM3pi_TEEG_Scaled,         "e^{+}e^{-}#gamma", "l");
  legd_cv -> AddEntry(h1d_IM3pi_TMCREST_Scaled,      "MC rest", "l");
  
  legd_cv -> Draw("Same");
  
  legtextsize(legd_cv, 0.05);

  p2 -> cd();

  hresidul -> GetYaxis() -> SetNdivisions(505);
  hresidul -> GetYaxis() -> SetTitleSize(25);
  hresidul -> GetYaxis() -> SetTitleFont(43);
  hresidul -> GetYaxis() -> SetTitleOffset(1.2);
  hresidul -> GetYaxis() -> SetLabelFont(43); // Absolute front size in pixel (precision 3)
  hresidul -> GetYaxis() -> SetLabelSize(15);
  hresidul -> GetYaxis() -> CenterTitle();
  hresidul -> GetYaxis() -> SetTitle("Residual");
  hresidul -> GetYaxis() -> SetRangeUser(residul_min, residul_max);
  //hresidul -> GetYaxis() -> SetTitle("(N_{d} - N_{mc})/#sqrt{N_{d}}");
  hresidul -> GetXaxis() -> SetLabelFont(43); // Absolute font size in pixel (precision 3)
  //hresidul -> GetXaxis() -> SetTitleOffset(1.2);
  //hresidul -> GetXaxis() -> SetTitleSize(45);
  hresidul -> GetXaxis() -> SetRangeUser(600., 950.);
  hresidul -> GetXaxis() -> SetLabelSize(15);
  hresidul -> GetXaxis() -> SetLabelOffset(0.03);
  hresidul -> GetXaxis() -> SetTitleSize(0.13);
  hresidul -> GetXaxis() -> SetTitle(xtit + " " + unit);
  hresidul -> GetXaxis() -> CenterTitle();
  hresidul -> SetStats(0);      // No statistics on lower plot
  hresidul -> SetMarkerStyle(21);
  hresidul -> SetMarkerSize(0.5);
  hresidul -> Draw("E0");

  // plot residul distr.
  const double binwidth_residul = getbinwidth(hresidul_distr);

  hresidul_distr -> SetMarkerStyle(21);
  hresidul_distr -> SetMarkerSize(0.7);
  hresidul_distr -> GetYaxis() -> SetTitle(TString::Format("Entries/%0.2f", binwidth_residul));
  hresidul_distr -> GetYaxis() -> CenterTitle();
  //hresidul_distr -> GetYaxis() -> SetRangeUser(0.1, 1200.);
  hresidul_distr -> GetYaxis() -> SetNdivisions(505);
  hresidul_distr -> GetYaxis() -> SetTitleSize(25);
  hresidul_distr -> GetYaxis() -> SetLabelSize(20);
  hresidul_distr -> GetYaxis() -> SetTitleFont(43);
  //hresidul_distr -> GetYaxis() -> SetTitleOffset(1.2);
  hresidul_distr -> GetYaxis() -> SetLabelFont(43); // Absolute font size in pixel (precision 3)
  hresidul_distr -> GetYaxis() -> SetTitleOffset(1.3);
  
  hresidul_distr -> GetXaxis() -> SetTitle(xtit + " Residual");
  //hresidul_distr -> GetXaxis() -> SetRangeUser(-30., 30.);
  hresidul_distr -> GetXaxis() -> SetTitleOffset(1.2);
  hresidul_distr -> GetXaxis() -> SetLabelSize(0.03);
  hresidul_distr -> GetXaxis() -> CenterTitle();

  TCanvas *cv_residul = new TCanvas("cv_residul", "cv_residul", 0, 0, 700, 700);

  hresidul_distr -> Draw();
  gfratio_fit -> Draw("Same");

  //pt1 -> Draw("Same");
  //pt2 -> Draw("Same");
  //pt3 -> Draw("Same");
  //pt4 -> Draw("Same");

  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(11);
  gStyle->SetStatX(0.85);
  gStyle->SetStatY(0.85);

  TPaveStats *pp = (TPaveStats*)hresidul_distr -> GetListOfFunctions() -> FindObject("stats");
  pp -> SetName("mystats");

  format_h(hresidul_distr, 4, 2);

  hresidul_distr -> SetStats(0);
  //hresidul_distr_sfw -> GetLineWith("mean") -> SetTextColor(kBlue);
  pp -> Draw("Same");

  

  return 0;

}
