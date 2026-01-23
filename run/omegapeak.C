
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

  // get root input files from: ~/Desktop/analysis_root_v6/sf_2g/plot_hist.C (smeared signal 3pi invariant mass distribution)
  // 

  TFile* intree = new TFile("/home/bo/Desktop/analysis_root_v6/sf_2g/plot_hist.root");
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
    
    cout << "tree" << i << ": classnm = " << classnm_tree << ", objnm = " << objnm_tree << endl;
    
  }

  TH1D* hIM3pi_smeared = (TH1D *)intree -> Get("hIM3pi_isr3pi_sc_corr");

  // Fit data
  double xmin = hIM3pi_smeared -> GetXaxis() -> GetXmin();
  double xmax = hIM3pi_smeared -> GetXaxis() -> GetXmax();
  const int npar = 5;
  double rms = hIM3pi_smeared -> GetRMS();
  double mean = hIM3pi_smeared -> GetMean();
  double peak = hIM3pi_smeared -> GetMaximum();
  double fitpara[npar] = {peak, mean, rms, 1, 1};
  
  cout << "mass range = [" << xmin << ", " << xmax << "] MeV/c^2 \n"
       << "rms = " << rms << ", mean = " << mean << ", peak = " << peak << "\n";

  const double fit_min = 774, fit_max = 792;
  const double binsize = hIM3pi_smeared -> GetNbinsX(); //cout<<"binsize = " << binsize << endl;
  double binwidth=(xmax-xmin)/binsize; //cout<<width<<endl;

  
  TF1 *fitFcn = new TF1("fitFcn", fitFunction, fit_min, fit_max, 5);
  fitFcn -> SetParameters(fitpara);
  fitFcn -> SetLineColor(kBlue);
  fitFcn -> SetLineWidth(2);
  fitFcn -> SetNpx(10000);
    
  TFitResultPtr r = hIM3pi_smeared -> Fit("fitFcn", "LM0", "", fit_min, fit_max);
  TF1 *fit_fun = hIM3pi_smeared -> GetFunction("fitFcn");
  fit_fun -> SetLineColor(kRed);
  fit_fun -> SetLineStyle(1);
  fit_fun -> SetLineWidth(2);
  fit_fun -> SetNpx(10000);

  /// Plots
  TCanvas *cv = new TCanvas("cv", "", 700, 700);

  //cv -> SetBottomMargin(0.2);//0.007
  cv -> SetLeftMargin(0.15);

  //hIM3pi_isr3pi_true -> GetYaxis() -> SetTitleOffset(1.2);
  hIM3pi_smeared -> GetYaxis() -> SetLabelSize(0.04);
  hIM3pi_smeared -> GetYaxis() -> SetTitleSize(0.05);
  hIM3pi_smeared -> GetYaxis() -> SetTitle(TString::Format("Entries/%0.2f", binwidth) + " [MeV/c^{2}]");
  hIM3pi_smeared -> GetYaxis() -> CenterTitle();
  
  hIM3pi_smeared -> GetXaxis() -> SetTitle("M_{3#pi} [MeV/c^{2}]");
  hIM3pi_smeared -> GetXaxis() -> CenterTitle();
  hIM3pi_smeared -> GetXaxis() -> SetTitleSize(0.04);
  hIM3pi_smeared -> GetXaxis() -> SetTitleOffset(1.2);
  
  hIM3pi_smeared -> Draw();
  fit_fun -> Draw("Same");
  
  return 0;
}
