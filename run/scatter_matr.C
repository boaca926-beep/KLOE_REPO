//#include "smearmatr.h"
//#include "crx3pi.h"
#include "../header/sm_para.h"

const double fit_min = 732;
const double fit_max = 832;

const double IM3pi_min = 380; //
const double IM3pi_max = 1020; //

const TString hist_nm = "_corred";
//const TString infile_nm = "/home/bo/Desktop/analysis/crx3pi/output_norm/crx3pi0.root";
const TString infile_nm = "/home/bo/Desktop/analysis_root_v6/input_vertex_TDATA_bkgrej";
TFile *infile = new TFile(infile_nm + "/cut/tree_pre.root");
TFile *inputfile = new TFile(infile_nm + "/input/sig.root"); 
TFile *f_hist = new TFile(infile_nm + "/hist/hist.root"); 
TList *HIM3pi_fit = (TList *)f_hist -> Get("HIM3pi_fit");
TH1D *hdata = (TH1D *)HIM3pi_fit -> FindObject("h1d_IM3pi_TDATA") -> Clone();
//hdata -> Draw();

TH2D *hcorrmatrix_tuned;
TH2D *hcorrmatrix_good;

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

//
double DetectorEvent(double mTrue) {
  // smear by double-gaussian
  if(rnd->Rndm()>frac) {
    return rnd->Gaus(mTrue+smallBias,smallSigma);
  } else {
    return rnd->Gaus(mTrue+wideBias,wideSigma);
  }
}

//
TH2D *getCorrMatrix_tuned(int bins, double xmin, double xmax) {

  TTree *TISR3PI_SIG = (TTree*)infile -> Get("TISR3PI_SIG");

  // create corrlation matrix
  TH2D *h2d_corr = new TH2D("", "", bins, xmin, xmax, bins, xmin, xmax);
  
  //cout << "correlation matrix is created!!! " << TISR3PI_SIG -> GetName() << endl;

  double IM3pi_true = 0.;
  double IM3pi_corred = 0.;
  double PARA[5] = {frac, smallBias, smallSigma, wideBias, wideSigma};
  
  for (Int_t irow = 0; irow < TISR3PI_SIG -> GetEntries(); irow++) {// loop chain

    TISR3PI_SIG -> GetEntry(irow);
    
    IM3pi_true = TISR3PI_SIG -> GetLeaf("Br_IM3pi_true") -> GetValue(0);
    
    //IM3pi_corred = DetectorEvent(TMath::Abs(IM3pi_true));
    IM3pi_corred = DetectorEvent_fcn(TMath::Abs(IM3pi_true), PARA);
    //cout << IM3pi_corred << ", " << DetectorEvent(TMath::Abs(IM3pi_true), para) << endl;
    
    h2d_corr -> Fill(IM3pi_true, IM3pi_corred);
    
  }
  
  return h2d_corr;
  
}

//
TH2D *getCorrMatrix_good(int bins, double xmin, double xmax) {// good signal events

  TTree *ALLCHAIN_CUT = (TTree*)inputfile -> Get("ALLCHAIN_CUT");

  // create corrlation matrix
  TH2D *h2d_corr = new TH2D("", "", bins, xmin, xmax, bins, xmin, xmax);
  
  //cout << "correlation matrix is created!!! " << TISR3PI_SIG -> GetName() << endl;

  int phid = 9999;
  int bkg_indx = 9999;
  int sig_type = 9999;
  int recon_indx = 9999;

  double IM3pi = 0., IM3pi_true = 0.;
  double IM3pi_corred = 0.;
  double PARA[5] = {frac, smallBias, smallSigma, wideBias, wideSigma};
  
  for (Int_t irow = 0; irow < ALLCHAIN_CUT -> GetEntries(); irow++) {// loop chain

    ALLCHAIN_CUT -> GetEntry(irow);

    phid = ALLCHAIN_CUT -> GetLeaf("Br_phid") -> GetValue(0);
    sig_type = ALLCHAIN_CUT -> GetLeaf("Br_sig_type") -> GetValue(0);
    bkg_indx = ALLCHAIN_CUT -> GetLeaf("Br_bkg_indx") -> GetValue(0);
    recon_indx = ALLCHAIN_CUT -> GetLeaf("Br_recon_indx") -> GetValue(0);

    IM3pi_true = ALLCHAIN_CUT -> GetLeaf("Br_IM3pi_true") -> GetValue(0);
    IM3pi = ALLCHAIN_CUT -> GetLeaf("Br_IM3pi_7C") -> GetValue(0);
    
    TH2D *h2d_corr = new TH2D("", "", bins, xmin, xmax, bins, xmin, xmax);

    if (irow > 1e4) break;
    
    if (phid == 0) {// e+ e- -> omega pi0
      
      //hist -> Fill(IM3pi);
      //h2d_scatter -> Fill(IM3pi_true, IM3pi);
      
      if (recon_indx == 2 && bkg_indx == 1) {
	h2d_corr -> Fill(IM3pi_true, IM3pi);
	cout << IM3pi_true << ", " << IM3pi << endl;
	
	//h1d_resol_Eisr -> Fill(Eisr_resol);
      }
      else {
	//hist_wrong -> Fill(IM3pi);
	//h2d_scatter_wrong -> Fill(IM3pi_true, IM3pi);
      }

    }
    //cout << "phid: " << phid << ", sig_type: " << sig_type << endl;
    //h2d_corr -> Fill(IM3pi_true, IM3pi_corred);
    
  }
  
  return h2d_corr;
  
}

int scatter_matr() {

  gErrorIgnoreLevel = kError;
  TGaxis::SetMaxDigits(4);
  gStyle->SetOptStat(0);
  //gROOT->SetBatch(kTRUE);
  gStyle->SetPaintTextFormat("4.0f");

  cout << "Plotting IM3pi MC true v.s. recon ... infile: " << infile_nm << endl;

  TIter next_tree(infile -> GetListOfKeys());

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

  int binsize = hdata -> GetNbinsX();
  hcorrmatrix_tuned = getCorrMatrix_tuned(binsize, IM3pi_min, IM3pi_max);
  hcorrmatrix_good = getCorrMatrix_good(binsize, IM3pi_min, IM3pi_max);
  hcorrmatrix_good -> Draw();
  //hcorrmatrix_tuned -> Draw();
  /*
  //TH2D * h2d = (TH2D*)infile -> Get("h2d_scatter" + hist_nm); // smeared matrix
  TH2D * h2d = (TH2D*)infile -> Get("h2d_scatter");  // reconstructed matrix

  double norm_factor = h2d -> Integral();
  //h2d -> Scale(1. / norm_factor);
  cout << norm_factor << ", h2d:" << h2d -> GetName() << endl;
  
  // plot

  TLine *line1 = new TLine(fit_min, IM3pi_max, fit_max, IM3pi_max); //upper horiz
  line1 -> SetLineColor(kBlack);
  line1 -> SetLineWidth(3);

  TLine *line2 = new TLine(fit_min, IM3pi_min, fit_max, IM3pi_min); //lower horiz
  line2 -> SetLineColor(kBlack);
  line2 -> SetLineWidth(3);

  TLine *line3 = new TLine(fit_min, IM3pi_min, fit_min, IM3pi_max); //left vertical
  line3 -> SetLineColor(kBlack);
  line3 -> SetLineWidth(3);

  TLine *line4 = new TLine(fit_max, IM3pi_min, fit_max, IM3pi_max); //right vertical
  line4 -> SetLineColor(kBlack);
  line4 -> SetLineWidth(3);

  double entries = h2d -> GetEntries();
  
  TCanvas * cv = new TCanvas("cv", "scatter" + hist_nm, 0, 0, 700, 700);
  cv -> SetBottomMargin(0.15);
  cv -> SetLeftMargin(0.15);
  cv -> SetRightMargin(0.15);
  
  char display1[50];

  TPaveText *pt1 = new TPaveText(0.15, 0.82, 0.25, 0.83, "NDC");
  
  pt1 -> SetTextSize(0.05);
  pt1 -> SetFillColor(0);
  pt1 -> SetTextAlign(12);

  //sprintf(display, test);
  sprintf(display1,"Entries %0.0f ", entries);
  pt1 -> AddText(display1);

  //const double axis_min = 550., axis_max = 900.;
  
  h2d -> GetZaxis() -> SetLabelSize(0.045);

  h2d -> GetXaxis() -> CenterTitle();
  h2d -> GetXaxis() -> SetTitle("M^{true}_{3#pi} [MeV/c^{2}]");
  h2d -> GetXaxis() -> SetTitleOffset(1.0);
  h2d -> GetXaxis() -> SetTitleSize(0.06);
  h2d -> GetXaxis() -> SetLabelSize(0.045);
  //h2d -> GetXaxis() -> SetRangeUser(axis_min, axis_max);
  
  h2d -> GetYaxis() -> CenterTitle();
  h2d -> GetYaxis() -> SetTitle("M^{rec}_{3#pi} [MeV/c^{2}]");
  //h2d -> GetYaxis() -> SetTitle("M^{smear}_{3#pi} [MeV/c^{2}]");
  h2d -> GetYaxis() -> SetTitleOffset(1.2);
  h2d -> GetYaxis() -> SetTitleSize(0.06);
  h2d -> GetYaxis() -> SetLabelSize(0.045);
  //h2d -> GetYaxis() -> SetRangeUser(axis_min, axis_max);
  
  //h2d -> Draw("TEXT0COL");
  h2d -> Draw("COLZ");
  //pt1 -> Draw("Same");
  //line1 -> Draw("Same");
  //line2 -> Draw("Same");
  line3 -> Draw("Same");
  line4 -> Draw("Same");

  gPad -> SetLogz();
  // save
  //TFile *fout = new TFile("./plots/smearmatr.root", "recreate");

  //cout << cv_nm << endl;
  //cv -> SaveAs("./plots/scatter" + hist_nm +  ".pdf");
  cv -> SaveAs("./plots/scatter.pdf");
  
  //fout -> Close();

  */  
  return 0;

}
