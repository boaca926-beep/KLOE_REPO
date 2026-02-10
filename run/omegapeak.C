#include "../header/plot.h"
#include "../header/sm_para.h"

TRandom *generator = new TRandom();
TRandom *rnd=0;
rnd=new TRandom3();

//
double DetectorEvent_fcn(double m, double *para_tmp) {
  // m[0]: mTrue
  // para[0] = frac
  // para[1] = smallBias
  // para[2] = smallSigma
  // para[3] = wideBias
  // para[4] = wideSigma

  // smear by double-gaussian
  if(rnd->Rndm()>para_tmp[0]) {
    return rnd->Gaus(m+para_tmp[1],para_tmp[2]);
  } else {
    return rnd->Gaus(m+para_tmp[3],para_tmp[4]);
  }
}

// Linear background function
Double_t background(Double_t *x, Double_t *par) {
  return par[0] + par[1]*x[0];
}

// Gaussian peak function with 3 parameters
double gaussianPeak(double *x, double *par) {
  double r1 = (x[0] - par[1]) / par[2];
  double fitval = par[0] * TMath::Exp(-0.5 * r1 * r1);
  //cout << par[0] << ", " << par[1] << ", " << par[2] << endl;
  return fitval;
}

// Sum of background and peak function
Double_t fitFunction(Double_t *x, Double_t *par) {
  //double *p1 = &par[0];
  //double *p2 = &par[2];
  //cout << par[0] << endl;
  //return gaussianPeak(x, par);
  return background(x, &par[3]) + gaussianPeak(x, par);
}

int omegapeak() {

  //gROOT->SetBatch(kTRUE);  
  gErrorIgnoreLevel = kError;
  TGaxis::SetMaxDigits(4);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  //gROOT->SetStyle("Plain");
  //gROOT->ForceStyle();
  
  // get root input files from: ~/Desktop/analysis_root_v6/sf_2g/plot_hist.C (smeared signal 3pi invariant mass distribution)
  // /home/bo/Desktop/analysis_root_v6/sf_1g/plot_hist.C (reconstructed signal 3pi invariant mass distribution)
  
  // input smeared
  TFile* intree_smear = new TFile("/home/bo/Desktop/analysis_root_v6/sf_2g/plot_hist.root");
  cout << intree_smear -> GetName() << endl;
  
  TIter next_tree(intree_smear -> GetListOfKeys());

  TString objnm_tree, classnm_tree;

  int i = 0;
  TKey *key;
  
  while ( (key = (TKey *) next_tree() ) ) {
    
    i ++;
    
    objnm_tree   =  key -> GetName();
    classnm_tree = key -> GetClassName();
    key -> GetSeekKey();
    
    cout << "tree" << i << ": classnm = " << classnm_tree << ", objnm = " << objnm_tree << endl;
    
  }

  TH1D* hIM3pi_smeared = (TH1D *)intree_smear -> Get("hIM3pi_isr3pi_sc_corr");

  // Fit smeared
  double xmin = hIM3pi_smeared -> GetXaxis() -> GetXmin();
  double xmax = hIM3pi_smeared -> GetXaxis() -> GetXmax();
  const int npar = 5;
  double rms = hIM3pi_smeared -> GetRMS();
  double mean = hIM3pi_smeared -> GetMean();
  double peak = hIM3pi_smeared -> GetMaximum();
  double fitpara[npar] = {peak, mean, rms, 1, 1};
  
  cout << "smear mass range = [" << xmin << ", " << xmax << "] MeV/c^2 \n"
       << "rms = " << rms << ", mean = " << mean << ", peak = " << peak << "\n";

  const double fit_min = 774, fit_max = 792;
  const double binsize = hIM3pi_smeared -> GetNbinsX(); //cout<<"binsize = " << binsize << endl;
  double binwidth=(xmax-xmin)/binsize; //cout<<width<<endl;

  TF1 *fitFcn = new TF1("fitFcn", fitFunction, fit_min, fit_max, 5);
  fitFcn -> SetParameters(fitpara);

  
  //fitFcn -> SetLineColor(kBlue);
  //fitFcn -> SetLineWidth(2);
  //fitFcn -> SetNpx(10000);
  
  TFitResultPtr r = hIM3pi_smeared -> Fit("fitFcn", "LM0", "", fit_min, fit_max);

  TF1 *fit_fun = hIM3pi_smeared -> GetFunction("fitFcn");
  fit_fun -> SetLineColor(kRed);
  fit_fun -> SetLineStyle(1);
  fit_fun -> SetLineWidth(2);
  fit_fun -> SetNpx(10000);
  
  // input reconstructed
  TFile* intree_rec = new TFile("/home/bo/Desktop/analysis_root_v6/sf_1g/plot_hist.root");

  cout << intree_rec -> GetName() << endl;
  
  TIter next_tree_rec(intree_rec -> GetListOfKeys());

  while ( (key = (TKey *) next_tree_rec() ) ) {
    
    i ++;
    
    objnm_tree   =  key -> GetName();
    classnm_tree = key -> GetClassName();
    key -> GetSeekKey();
    
    cout << "tree" << i << ": classnm = " << classnm_tree << ", objnm = " << objnm_tree << endl;
    
  }

  TH1D* hIM3pi_rec = (TH1D *)intree_rec -> Get("hIM3pi_isr3pi_sc_rec");
  
  // Fit reconstructed
  double xmin_rec = hIM3pi_rec -> GetXaxis() -> GetXmin();
  double xmax_rec = hIM3pi_rec -> GetXaxis() -> GetXmax();
  double rms_rec = hIM3pi_rec -> GetRMS();
  double mean_rec = hIM3pi_rec -> GetMean();
  double peak_rec = hIM3pi_rec -> GetMaximum();
  double fitpara_rec[npar] = {peak_rec, mean_rec, rms_rec, 1, 1};

  cout << "rec mass range = [" << xmin_rec << ", " << xmax_rec << "] MeV/c^2 \n"
       << "rms_rec = " << rms_rec << ", mean_rec = " << mean_rec << ", peak_rec = " << peak_rec << "\n";

  const double fit_min_rec = 775, fit_max_rec = 793;

  TF1 *fitFcn_rec = new TF1("fitFcn_rec", fitFunction, fit_min_rec, fit_max_rec, 5);
  fitFcn_rec -> SetParameters(fitpara_rec);

  TFitResultPtr r_rec = hIM3pi_rec -> Fit("fitFcn_rec", "LM0", "", fit_min_rec, fit_max_rec);
  TF1 *fit_fun_rec = hIM3pi_rec -> GetFunction("fitFcn_rec");
  fit_fun_rec -> SetLineColor(kGreen);
  fit_fun_rec -> SetLineStyle(1);
  fit_fun_rec -> SetLineWidth(2);
  fit_fun_rec -> SetNpx(10000);

  
  double rms_diff = TMath::Abs(rms - rms_rec);
  cout << "rms difference: (rms - rms_rec) = " << rms_diff << "%\n";

  double width_smear = fitFcn -> GetParameter(2);
  double width_smear_err = fitFcn -> GetParError(2);
  double width_rec = fitFcn_rec -> GetParameter(2);
  double width_rec_err = fitFcn_rec -> GetParError(2);
  double width_diff = TMath::Abs(width_smear - width_rec);
  double Z_width = width_diff / TMath::Sqrt(width_rec_err * width_rec_err + width_smear * width_smear);
  
  cout << "width_smear = " << width_smear << "+/-" << width_smear_err << "\n"
       << "width_rec = " << width_rec << "+/-" << width_rec_err << "\n"
       << "width_diff = " << width_diff << ", rel. diff. = " << width_diff / width_smear * 100. << "%\n"
       << "Z_width = " << Z_width << endl;

  
  
  /// Plots
  TCanvas *cv = new TCanvas("cv", "", 800, 700);
  //cv -> SetBottomMargin(0.2);//0.007
  cv -> SetLeftMargin(0.15);
  cv -> SetLeftMargin(0.15);

  TPaveText *pt1 = new TPaveText(0.2, 0.85, 0.33, 0.87, "NDC");
  PteAttr(pt1);
  pt1 -> SetTextSize(0.04);
  pt1 -> SetTextColor(kRed);
  pt1 -> AddText(Form("#Gamma=%0.2f#pm%0.2f MeV", width_smear, width_smear_err));

  TPaveText *pt2 = new TPaveText(0.2, 0.81, 0.33, 0.83, "NDC");
  PteAttr(pt2);
  pt2 -> SetTextSize(0.04);
  pt2 -> SetTextColor(kGreen);
  pt2 -> AddText(Form("#Gamma=%0.2f#pm%0.2f MeV", width_rec, width_rec_err));

  TPaveText *pt3 = new TPaveText(0.2, 0.77, 0.33, 0.79, "NDC");
  PteAttr(pt3);
  pt3 -> SetTextSize(0.04);
  //pt3 -> SetTextColor(kGreen);
  //pt3 -> AddText(Form("Z=%0.2f", Z_width));
  
  
  //hIM3pi_smeared -> GetXaxis() -> SetTitle("M_{3#pi} [MeV/c^{2}]");
  //hIM3pi_smeared -> GetXaxis() -> SetNdivisions(505);
  hIM3pi_smeared -> GetXaxis() -> SetTitle("M_{3#pi} [MeV/c^{2}]");
  hIM3pi_smeared -> GetXaxis() -> CenterTitle();
  hIM3pi_smeared -> GetXaxis() -> SetTitleSize(0.05);
  hIM3pi_smeared -> GetXaxis() -> SetTitleOffset(0.9);
  //hIM3pi_smeared -> GetXaxis() -> SetLabelOffset(0.01);
  hIM3pi_smeared -> GetXaxis() -> SetLabelSize(0.04);//0.03
  //hIM3pi_smeared -> GetXaxis() -> SetRangeUser(0., 40.); //chi2
  //hIM3pi_smeared -> GetXaxis() -> SetRangeUser(20., 140.); //pi0angle
  //hIM3pi_smeared -> GetXaxis() -> SetRangeUser(0.5, 1.);
  //hIM3pi_smeared -> GetXaxis() -> SetRangeUser(650., 950.); // 3pi mass omega region
  
  //hIM3pi_smeared -> GetYaxis() -> SetLabelFont(43); // Absolute front size in pixel (precision 3)
  //hIM3pi_smeared -> GetYaxis() -> SetLabelSize(0.05);
  hIM3pi_smeared -> GetYaxis() -> SetRangeUser(0., peak * 1.5); // h_max * 1e3
  //hIM3pi_smeared -> GetYaxis() -> SetRangeUser(0.01, ymax * 1e2); // 3pi mass full region
  //hIM3pi_smeared -> GetYaxis() -> SetRangeUser(0., 5e4); // 3pi mass omega region  
  //hIM3pi_smeared -> GetYaxis() -> SetTitle(TString::Format("Entries/%0.2f", binwidth) + " " + unit);
  hIM3pi_smeared -> GetYaxis() -> SetTitle("Events");
  hIM3pi_smeared -> GetYaxis() -> CenterTitle();
  hIM3pi_smeared -> GetYaxis() -> SetTitleSize(0.06);
  hIM3pi_smeared -> GetYaxis() -> SetTitleOffset(1.2);
  
  format_h(hIM3pi_rec, kBlack, 2);
  
  hIM3pi_smeared -> Draw();
  hIM3pi_rec -> Draw("Same");
  fit_fun_rec -> Draw("Same");
  fit_fun -> Draw("Same");

  pt1 -> Draw("Same");
  pt2 -> Draw("Same");
  pt3 -> Draw("Same");
  
  TLegend *legd_cv = new TLegend(0.6, 0.65, 0.9, 0.9);
  
  legd_cv -> SetTextFont(132);
  legd_cv -> SetFillStyle(0);
  legd_cv -> SetBorderSize(0);
  legd_cv -> SetNColumns(1);
  
  legd_cv -> AddEntry(hIM3pi_smeared, "Smeared", "lep");
  legd_cv -> AddEntry(hIM3pi_rec, "Reconstructed", "l");
  
  legd_cv -> Draw("Same");
  
  legtextsize(legd_cv, 0.04);

  cv -> SaveAs("omegapeak.pdf");
  
  return 0;
}
