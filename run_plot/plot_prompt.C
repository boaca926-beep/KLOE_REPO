#include "../hist.h"

int plot_prompt() {

  gErrorIgnoreLevel = kError;
  TGaxis::SetMaxDigits(3);
  gStyle->SetOptTitle(0);
  gStyle->SetStatBorderSize(0);
  gStyle->SetOptStat(0);

  // get histo
  TFile *intree = new TFile("sig_pre.root");

  TH1D *prompt_distr = (TH1D *)intree -> Get("prompt_distr");
  
  formatfill_h(prompt_distr, 4, 3001);

  
  TIter next_tree(intree -> GetListOfKeys());
  
  TString objnm_tree, classnm_tree;
  
  int i = 0;

  TKey *key;
  
  while ( (key = (TKey *) next_tree() ) ) {// start tree while lop
    
    i ++;
    
    objnm_tree   =  key -> GetName();
    classnm_tree = key -> GetClassName();
    key -> GetSeekKey();
    
    cout << " tree" << i << ": classnm = " << classnm_tree << ", objnm = " << objnm_tree << endl;
    
  }

  // normalize histo
  const double bins = prompt_distr -> GetNbinsX(); 
  double entries  = prompt_distr -> Integral(1, bins);
  prompt_distr -> Scale(1. / entries);

  cout << "bins = " << bins << ", entries = " << entries << endl;

  // plots
  TCanvas *cv = new TCanvas("prompt_distr", "", 0, 0, 1000, 800);
  cv -> SetBottomMargin(0.15);//0.007
  cv -> SetLeftMargin(0.15);
  
  prompt_distr -> SetMarkerStyle(21);
  prompt_distr -> SetMarkerSize(0.7);
  prompt_distr -> GetYaxis() -> SetTitle("Probabilities");
  prompt_distr -> GetYaxis() -> SetTitleOffset(1.);
  prompt_distr -> GetYaxis() -> SetTitleSize(0.07);
  prompt_distr -> GetYaxis() -> SetLabelSize(0.05);
  prompt_distr -> GetYaxis() -> CenterTitle();
  
  //prompt_distr -> GetYaxis() -> SetRangeUser(0.1, ymax);
  prompt_distr -> GetXaxis() -> SetTitle("Number of prompt photons n_{#gamma}");
  prompt_distr -> GetXaxis() -> SetRangeUser(0, 8);
  prompt_distr -> GetXaxis() -> SetTitleOffset(0.9);
  prompt_distr -> GetXaxis() -> SetTitleSize(0.07);
  prompt_distr -> GetXaxis() -> SetLabelSize(0.05);
  prompt_distr -> GetXaxis() -> CenterTitle();

  prompt_distr -> Draw("hist");

  gPad -> SetLogy();

  cv -> SaveAs("prompt.pdf");
  

  return 0;
  
}


void plot_distr() {

  
}
