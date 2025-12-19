#include "../header/plot.h"

const TString hist_nm = "_corred";
const TString infile_nm = "../../../analysis/crx3pi/output_norm/crx3pi0.root";

const double xrange1 = 740., xrange2 = 820.; // mc true
const double yrange1 = 760., yrange2 = 800.; // mc recon


int smearmatr_ppt() {

  gErrorIgnoreLevel = kError;
  TGaxis::SetMaxDigits(4);
  gStyle->SetOptStat(0);
  //gROOT->SetBatch(kTRUE);
  gStyle->SetPaintTextFormat("4.4f");

  cout << "Plotting smearing matrix ... " << endl;

  TFile *infile = new TFile(infile_nm);

  TIter next_tree(infile -> GetListOfKeys());

  TString objnm_tree, classnm_tree;
  
  int i = 0;
  TKey *key;
  
  while ( (key = (TKey *) next_tree() ) ) {// start tree while lop
    
    i ++;
    
    objnm_tree   =  key -> GetName();
    classnm_tree = key -> GetClassName();
    key -> GetSeekKey();
    
    //cout << "tree" << i << ": classnm = " << classnm_tree << ", objnm = " << objnm_tree << endl;
    
  }

  TH1D * hprobsum = (TH1D *) infile -> Get("hprobsum");
  
  TH2D * hsmearmatr = (TH2D *) infile -> Get("hsmearmatr");
  TH2D * h2d_scatter_corred = (TH2D *) infile -> Get("h2d_scatter_corred");

  TLine *linex = new TLine(740., 800., 822., 800.); // horiz
  linex -> SetLineColor(kBlack);
  linex -> SetLineWidth(3);

  TLine *liney = new TLine(740., 760., 822., 760.); // horiz
  liney -> SetLineColor(kBlack);
  liney -> SetLineWidth(3);

  //
  TCanvas * cv1 = new TCanvas("cv1", "smearing matrix PPT (signal sfw corrected)", 0, 0, 1300, 800);

  cv1 -> SetBottomMargin(0.16);
  cv1 -> SetLeftMargin(0.15);
  cv1 -> SetRightMargin(0.15);
  
  hsmearmatr -> GetXaxis() -> SetTitle("M^{true}_{3#pi} [MeV/c^{2}]"); //SetTitle("M^{TRUE}_{3#pi} [MeV/c^{2}]");
  hsmearmatr -> GetXaxis() -> SetTitleSize(0.07);
  
  hsmearmatr -> GetYaxis() -> SetTitle("M^{smear}_{3#pi} [MeV/c^{2}]"); //-> SetTitle("M^{REC}_{3#pi} [MeV/c^{2}]");
  hsmearmatr -> GetYaxis() -> SetTitleSize(0.07);
  
  hsmearmatr -> GetXaxis() -> SetTitleOffset(1.);
  hsmearmatr -> GetYaxis() -> SetTitleOffset(1.);

  hsmearmatr -> GetXaxis() -> SetRangeUser(xrange1, xrange2);
  hsmearmatr -> GetYaxis() -> SetRangeUser(xrange1, xrange2);

  hsmearmatr -> GetXaxis() -> SetLabelSize(0.05); //20, 0.03
  hsmearmatr -> GetYaxis() -> SetLabelSize(0.05);
  hsmearmatr -> GetZaxis() -> SetLabelSize(0.05);

  hsmearmatr -> GetXaxis() -> CenterTitle();
  hsmearmatr -> GetYaxis() -> CenterTitle();

  //hsmearmatr -> GetXaxis() -> SetLabelOffset(0.1);
  
  hsmearmatr -> SetStats(0);  

  hsmearmatr -> Draw("COLZ");
  linex -> Draw("Same");
  liney -> Draw("Same");
  
  gPad -> SetLogz();

  cv1 -> SaveAs("smearmatr_corred_PPT.pdf");
  
  return 0;

}
