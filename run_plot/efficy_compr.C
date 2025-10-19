// Compare efficiencies between the standard analysis and ECL correction with additional two background suppression critera: KSL and etagamm suppresion
#include "../header/plot.h"

int efficy_compr(){

  //gROOT->SetBatch(kTRUE);
  gErrorIgnoreLevel = kError;
  TGaxis::SetMaxDigits(4);
  gStyle->SetOptStat(0);

  // get histos
  TFile *intree1 = new TFile("/home/bo/Desktop/input_norm_TDATA/hist/hist.root"); // nominal
  TFile *intree2 = new TFile("/home/bo/Desktop/input_vertex_TDATA/hist/hist.root"); // ECL efficiency correction
  
  //cout << intree1 << endl;

  TIter next_tree(intree1 -> GetListOfKeys());

  TString objnm_tree, classnm_tree;
  
  int i = 0;
  TKey *key;
  
  while ( (key = (TKey *) next_tree() ) ) {// start tree while lop
    
    i ++;
    
    objnm_tree   =  key -> GetName();
    classnm_tree = key -> GetClassName();
    key -> GetSeekKey();
    
    cout << "tree" << i << ": classnm = " << classnm_tree << ", objnm = " << objnm_tree << endl;
    
  }

  TList *HSIG1 = (TList *) intree1 -> Get("HSIG");
  TH1D *hsig_true1 = (TH1D *) HSIG1 -> FindObject("hsig_true");

  TList *HSIG2 = (TList *) intree2 -> Get("HSIG");
  TH1D *hsig_true2 = (TH1D *) HSIG2 -> FindObject("hsig_true");

  hsig_true1->Sumw2();
  hsig_true2->Sumw2();
  
  //hsig_true2 -> Draw();
  //hsig_true1 -> Draw("Same");

  // Divide two histos
  TH1D *hratio = (TH1D*)hsig_true2 -> Clone();
  hratio -> Divide(hsig_true1);
  //hratio -> SetTitle("hratio");
  format_h(hratio, 4, 2);// color 45
  hratio -> SetMarkerStyle(21);
  hratio -> SetMarkerSize(0.5);
  hratio -> SetMarkerColor(4);

  // plot
  TCanvas *cv_ratio = new TCanvas("cv_ratio", "", 0, 0, 1000, 800);

  //cv_ratio -> SetBottomMargin(0.15);
  cv_ratio -> SetLeftMargin(0.1);
  //cv_ratio -> SetRightMargin(0.15);

  hratio -> GetYaxis() -> SetNdivisions(512);
  hratio -> GetYaxis() -> SetRangeUser(0.5, 1.0); 
  hratio -> GetYaxis() -> SetTitle("Efficiency ratio");
  hratio -> GetYaxis() -> SetTitleSize(0.05);
  hratio -> GetYaxis() -> SetTitleOffset(1.);
  hratio -> GetYaxis() -> SetLabelSize(0.035);
  //hratio -> GetYaxis() -> SetTitleFont(132);
  hratio -> GetYaxis() -> CenterTitle();
      
  hratio -> GetXaxis() -> SetTitle("M_{3#pi} [MeV/c^{2}]");
  hratio -> GetXaxis() -> SetTitleOffset(.8);
  hratio -> GetXaxis() -> SetTitleSize(0.05);
  hratio -> GetXaxis() -> SetLabelSize(0.04);
  //hratio -> GetXaxis() -> SetTitleFont(132);
  hratio -> GetXaxis() -> SetRangeUser(760., 800.);
  hratio -> GetXaxis() -> CenterTitle();
  
  hratio -> Draw();

  cv_ratio -> SaveAs("efficy_compr.pdf");
  
  return 0;
}

  
