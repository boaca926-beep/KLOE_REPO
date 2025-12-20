#include "../header/plot.h"
//#include "crx3piCompr.h"

//#include "Riostream.h"

// Compare KLOE/KLOE-2, SND (Phys.Rev.D 68 (2003) 052006, 2003.), BESIII (arXiv:1912.11208), and Babar (arXiv:2110.00520), CMD-2 (Phys.Lett.B 476 (2000) 33-39, 2000) 

TVectorD GetDiffVect(TGraphErrors *gf_kloe2, double m3pi_target, double crx3pi_target, double crx3piErr_target) {

  
  TVectorD vect(2);

  double vNout = gf_kloe2 -> GetN();
  double m3piKloe2, crx3piKloe2, crx3piKloe2Err;
  double m3piDiff_min = 999999., m3piDiff = 0.;
  double crx3piDiff_min = 0., crx3piDiff = 0.;
  double crx3piDiffErr_min = 0., crx3piDiffErr = 0.;
  int kloe2BinIndx_min = 0;

  for (Int_t k=0;k<vNout;k++) {
    gf_kloe2 -> GetPoint(k, m3piKloe2, crx3piKloe2);
    crx3piKloe2Err = gf_kloe2 -> GetErrorY(k);
    m3piDiff = TMath::Abs(m3piKloe2 - m3pi_target);
    crx3piDiff = crx3piKloe2 - crx3pi_target;
    crx3piDiffErr = TMath::Sqrt(crx3piKloe2Err * crx3piKloe2Err + crx3piErr_target * crx3piErr_target);
    //crx3piDiffErr = crx3piErr_target;
    
    if (m3piDiff < m3piDiff_min) {
      m3piDiff_min = m3piDiff;
      crx3piDiff_min = crx3piDiff;
      crx3piDiffErr_min = crx3piDiffErr;
      kloe2BinIndx_min = k + 1;
    }

    /*
    cout << "kloe2 bin " << k + 1 << "  m3pi = " << m3piKloe2
	 << " MeV/c^2  crx3pi = " << crx3piKloe2 << "+/-" << crx3piKloe2Err << ", target m3pi = " << m3pi_target << " MeV/c^2" << ", |m3piDiff| = " << m3piDiff << ", min |m3piDiff| = " << m3piDiff_min << ", kloe2 m3pi bin indx for m3piDiff_min = " << kloe2BinIndx_min << endl;
    */
    
  }

  //cout << "min m3piDiff at mass bin " << kloe2BinIndx_min << ", crx3piDiff_min = " << crx3piDiff_min << "+/-" << crx3piDiffErr_min << endl;

  /*
  vect(0) = crx3piDiff_min;
  vect(1) = crx3piDiffErr_min;
  */

  vect(0) = crx3piDiff_min;
  vect(1) = crx3piDiffErr_min;
  
  return vect;

}

int crx3piCompr() {

  //gROOT->SetBatch(kTRUE);  
  gErrorIgnoreLevel = kError;
  TGaxis::SetMaxDigits(4);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetStatBorderSize(0);
  gStyle->SetFitFormat("6.2g");

  /// input root files
  TFile* fr_cmd2 = new TFile("./plotOutput/cmd2Crx3pi.root", "READ");
  TFile* fr_snd = new TFile("./plotOutput/sndCrx3pi.root", "READ");
  TFile* fr_kloe2 = new TFile("./plotOutput/crx3piObs.root", "READ");
  TFile* fr_babar = new TFile("./plotOutput/Babar.root", "READ");
  TFile* fr_bes = new TFile("./plotOutput/bes.root", "READ");
  

  TIter next_tree( fr_snd -> GetListOfKeys());

  TString objnm_tree, classnm_tree;
  TKey *key;
  
  int i = 0;
  
  while ( (key = (TKey *) next_tree() ) ) {
    
    i ++;
    
    objnm_tree   =  key -> GetName();
    classnm_tree = key -> GetClassName();
    //key -> GetSeekKey();
    
    cout << "tree" << i << ": classnm = " << classnm_tree << ", objnm = " << objnm_tree << endl;
    
  }

  /// getting gfaphs
  // CMD2
  TGraphErrors *gf_cmd2 = (TGraphErrors*) fr_cmd2 -> Get("Table 1/Graph1D_y1");
  gf_cmd2 -> SetMarkerStyle(27);
  gf_cmd2 -> SetMarkerColor(3);
  gf_cmd2 -> SetLineColor(3);
  gf_cmd2 -> SetMarkerSize(0.7);

  // SND
  TGraphErrors *gf_snd = (TGraphErrors*) fr_snd -> Get("Table 1/Graph1D_y1");
  gf_snd -> SetMarkerStyle(25);
  gf_snd -> SetMarkerColor(2);
  gf_snd -> SetLineColor(2);
  gf_snd -> SetMarkerSize(0.7);

  // Kloe2
  TCanvas *cnv = (TCanvas*)fr_kloe2 -> Get("cv");
  TGraphErrors *gf_kloe2 = (TGraphErrors*) cnv -> FindObject("Graph0");
  gf_kloe2 -> SetMarkerStyle(20);
  gf_kloe2 -> SetMarkerColor(4);
  gf_kloe2 -> SetLineColor(4);
  gf_kloe2 -> SetMarkerSize(0.7);

  // save Kloe2 visible crs3pi data in root files

  // Babar
  TGraphErrors *gf_babar = (TGraphErrors*) fr_babar -> Get("Graph");
  gf_babar -> SetMarkerStyle(30);
  gf_babar -> SetMarkerColor(46);
  gf_babar -> SetLineColor(46);
  gf_babar -> SetMarkerSize(0.7);
  
  // BesIII
  TGraphErrors *gf_bes = (TGraphErrors*) fr_bes -> Get("Graph");
  gf_bes -> SetMarkerStyle(32);
  gf_bes -> SetMarkerColor(28);
  gf_bes -> SetLineColor(28);
  gf_bes -> SetMarkerSize(0.7);

  // pulls, diviation between KLOE-2 reult and results of other experiments
  // BaBar
  int nbins = gf_babar -> GetN();
  double CRX3PIDIFF_BABAR[1000], CRX3PIDIFFERR_BABAR[1000];
  double M3PI_BABAR[1000], M3PIERR_BABAR[1000];
  double m3pi = 0., crx3piBaBar = 0., crx3piBaBarErr = 0.;
  TVectorD diffVect(2);
  
  cout << "loop over BaBar data set, number of bins = " << nbins << endl;
  
  for (Int_t k = 0; k < nbins; k++) {

    //if (k > 3) break;

    gf_babar -> GetPoint(k, m3pi, crx3piBaBar);
    crx3piBaBarErr = gf_babar -> GetErrorY(k);
   
    diffVect = GetDiffVect(gf_kloe2, m3pi, crx3piBaBar, crx3piBaBarErr);
    M3PI_BABAR[k] = m3pi;
    M3PIERR_BABAR[k] = 0.;
    
    CRX3PIDIFF_BABAR[k] = diffVect(0);
    CRX3PIDIFFERR_BABAR[k] = diffVect(1);

    cout << "BaBar bin  " << k + 1 << "  m3pi = " << m3pi
	 << " MeV/c^2  crx3pi = " << crx3piBaBar << "+/-" << crx3piBaBarErr
	 << " crx3pi diff. = " << CRX3PIDIFF_BABAR[k] << "+/-" << CRX3PIDIFFERR_BABAR[k] << endl;
     
  }

  TGraphErrors * gf_babar_diff = new TGraphErrors(nbins, M3PI_BABAR, CRX3PIDIFF_BABAR, M3PIERR_BABAR, CRX3PIDIFFERR_BABAR);
  gf_babar_diff -> SetMarkerStyle(30);
  gf_babar_diff -> SetMarkerColor(46);
  gf_babar_diff -> SetLineColor(46);
  gf_babar_diff -> SetMarkerSize(0.7);
  
  // cmd-2
  nbins = gf_cmd2 -> GetN();
  double CRX3PIDIFF_CMD2[1000], CRX3PIDIFFERR_CMD2[1000];
  double M3PI_CMD2[1000], M3PIERR_CMD2[1000];
  double m3pi_cmd2 = 0.;
  double crx3piCMD2 = 0., crx3piCMD2Err = 0.;
  TVectorD diffVect_CMD2(2);
  
  cout << "loop over CMD-2 data set, number of bins = " << nbins << endl;

  for (Int_t k = 0; k < nbins; k++) {

    //if (k > 0) break;

    gf_cmd2 -> GetPoint(k, m3pi_cmd2, crx3piCMD2);
    crx3piCMD2Err = gf_cmd2 -> GetErrorY(k);
   
    diffVect_CMD2 = GetDiffVect(gf_kloe2, m3pi_cmd2, crx3piCMD2, crx3piCMD2Err);
    M3PI_CMD2[k] = m3pi_cmd2;
    M3PIERR_CMD2[k] = 0.;
    
    CRX3PIDIFF_CMD2[k] = diffVect_CMD2(0);
    CRX3PIDIFFERR_CMD2[k] = diffVect_CMD2(1);

    cout << "CMD2 bin  " << k + 1 << "  m3pi = " << m3pi_cmd2
	 << " MeV/c^2  crx3pi = " << crx3piCMD2 << "+/-" << crx3piCMD2Err
	 << " crx3pi diff. = " << CRX3PIDIFF_CMD2[k] << "+/-" << CRX3PIDIFFERR_CMD2[k] << endl;
     
  }

  TGraphErrors * gf_cmd2_diff = new TGraphErrors(nbins, M3PI_CMD2, CRX3PIDIFF_CMD2, M3PIERR_CMD2, CRX3PIDIFFERR_CMD2);
  gf_cmd2_diff -> SetMarkerStyle(27);
  gf_cmd2_diff -> SetMarkerColor(3);
  gf_cmd2_diff -> SetLineColor(3);
  gf_cmd2_diff -> SetMarkerSize(0.7);

  // besIII
  nbins = gf_bes -> GetN();
  double CRX3PIDIFF_BES[1000], CRX3PIDIFFERR_BES[1000];
  double M3PI_BES[1000], M3PIERR_BES[1000];
  double m3pi_bes = 0.;
  double crx3piBES = 0., crx3piBESErr = 0.;
  TVectorD diffVect_BES(2);
  
  cout << "loop over BESIII data set, number of bins = " << nbins << endl;

  for (Int_t k = 0; k < nbins; k++) {

    //if (k > 0) break;

    gf_bes -> GetPoint(k, m3pi_bes, crx3piBES);
    crx3piBESErr = gf_bes -> GetErrorY(k);
   
    diffVect_BES = GetDiffVect(gf_kloe2, m3pi_bes, crx3piBES, crx3piBESErr);
    M3PI_BES[k] = m3pi_bes;
    M3PIERR_BES[k] = 0.;
    
    CRX3PIDIFF_BES[k] = diffVect_BES(0);
    CRX3PIDIFFERR_BES[k] = diffVect_BES(1);

    cout << "BES bin  " << k + 1 << "  m3pi = " << m3pi_bes
	 << " MeV/c^2  crx3pi = " << crx3piBES << "+/-" << crx3piBESErr
	 << " crx3pi diff. = " << CRX3PIDIFF_BES[k] << "+/-" << CRX3PIDIFFERR_BES[k] << endl;
     
  }

  TGraphErrors * gf_bes_diff = new TGraphErrors(nbins, M3PI_BES, CRX3PIDIFF_BES, M3PIERR_BES, CRX3PIDIFFERR_BES);
  gf_bes_diff -> SetMarkerStyle(32);
  gf_bes_diff -> SetMarkerColor(28);
  gf_bes_diff -> SetLineColor(28);
  gf_bes_diff -> SetMarkerSize(0.7);

  // SND
  nbins = gf_snd -> GetN();
  double CRX3PIDIFF_SND[1000], CRX3PIDIFFERR_SND[1000];
  double M3PI_SND[1000], M3PIERR_SND[1000];
  double m3pi_snd = 0.;
  double crx3piSND = 0., crx3piSNDErr = 0.;
  TVectorD diffVect_SND(2);
  
  cout << "loop over SND data set, number of bins = " << nbins << endl;

  for (Int_t k = 0; k < nbins; k++) {

    //if (k > 0) break;

    gf_snd -> GetPoint(k, m3pi_snd, crx3piSND);
    crx3piSNDErr = gf_snd -> GetErrorY(k);
   
    diffVect_SND = GetDiffVect(gf_kloe2, m3pi_snd, crx3piSND, crx3piSNDErr);
    M3PI_SND[k] = m3pi_snd;
    M3PIERR_SND[k] = 0.;
    
    CRX3PIDIFF_SND[k] = diffVect_SND(0);
    CRX3PIDIFFERR_SND[k] = diffVect_SND(1);

    cout << "SND bin  " << k + 1 << "  m3pi = " << m3pi_bes
	 << " MeV/c^2  crx3pi = " << crx3piSND << "+/-" << crx3piSNDErr
	 << " crx3pi diff. = " << CRX3PIDIFF_SND[k] << "+/-" << CRX3PIDIFFERR_SND[k] << endl;
     
  }

  TGraphErrors * gf_snd_diff = new TGraphErrors(nbins, M3PI_SND, CRX3PIDIFF_SND, M3PIERR_SND, CRX3PIDIFFERR_SND);
  gf_snd_diff -> SetMarkerStyle(25);
  gf_snd_diff -> SetMarkerColor(2);
  gf_snd_diff -> SetLineColor(2);
  gf_snd_diff -> SetMarkerSize(0.7);

  /// plot
  const double xmin = 740.;
  const double xmax = 820.;

  //
  TCanvas *cv = new TCanvas("cv", "crx3pi comparison", 1200, 600);
  cv -> SetBottomMargin(0.2);//0.007
  cv -> SetLeftMargin(0.15);

  gf_kloe2 -> GetXaxis() -> CenterTitle();
  gf_kloe2 -> GetXaxis() -> SetTitle("#sqrt{s} [MeV]");
  gf_kloe2 -> GetXaxis() -> SetTitleOffset(1.0);
  gf_kloe2 -> GetXaxis() -> SetTitleSize(0.06);
  gf_kloe2 -> GetXaxis() -> SetLabelSize(0.045);
  gf_kloe2 -> GetXaxis() -> SetRangeUser(xmin, xmax); 
  gf_kloe2 -> GetXaxis() -> SetTitleFont(132);
  
  gf_kloe2 -> GetYaxis() -> CenterTitle();
  gf_kloe2 -> GetYaxis() -> SetTitle("#sigma_{3#pi} [nb]");
  gf_kloe2 -> GetYaxis() -> SetTitleOffset(1.2);
  gf_kloe2 -> GetYaxis() -> SetTitleSize(0.06);
  gf_kloe2 -> GetYaxis() -> SetLabelSize(0.045);
  gf_kloe2 -> GetYaxis() -> SetTitleFont(132);
  gf_kloe2 -> GetYaxis() -> SetRangeUser(0., 1700.); 
  
  gf_kloe2 -> Draw("ap");
  gf_snd -> Draw("p");
  gf_cmd2 -> Draw("p");
  gf_babar -> Draw("p");
  gf_bes -> Draw("p");
  
  
  TLegend *legd_cv = new TLegend(0.2, 0.5, 0.6, 0.9);
  
  legd_cv -> SetTextFont(132);
  legd_cv -> SetFillStyle(0);
  legd_cv -> SetBorderSize(0);
  legd_cv -> SetNColumns(1);

  legd_cv -> AddEntry(gf_kloe2, "This work", "lep");
  legd_cv -> AddEntry(gf_snd, "SND", "lep");
  legd_cv -> AddEntry(gf_cmd2, "CMD-2", "lep");
  legd_cv -> AddEntry(gf_babar, "BaBar", "lep");
  legd_cv -> AddEntry(gf_bes, "BESIII", "lep");
  
  legd_cv -> Draw("Same");
  
  legtextsize(legd_cv, 0.04);

  /*
  //
  TCanvas *cv_diff = new TCanvas("cv_diff", "crx3pi diviation comparison", 1000, 800);
  cv_diff -> SetBottomMargin(0.2);//0.007
  cv_diff -> SetLeftMargin(0.15);

  gf_babar_diff -> GetXaxis() -> CenterTitle();
  gf_babar_diff -> GetXaxis() -> SetTitle("#sqrt{s} [MeV]");
  gf_babar_diff -> GetXaxis() -> SetTitleOffset(1.0);
  gf_babar_diff -> GetXaxis() -> SetTitleSize(0.06);
  gf_babar_diff -> GetXaxis() -> SetLabelSize(0.045);
  gf_babar_diff -> GetXaxis() -> SetRangeUser(xmin, xmax); 
  gf_babar_diff -> GetXaxis() -> SetTitleFont(132);
  
  gf_babar_diff -> GetYaxis() -> CenterTitle();
  gf_babar_diff -> GetYaxis() -> SetTitle("Residuals [nb]");
  gf_babar_diff -> GetYaxis() -> SetTitleOffset(1.2);
  gf_babar_diff -> GetYaxis() -> SetTitleSize(0.06);
  gf_babar_diff -> GetYaxis() -> SetLabelSize(0.045);
  gf_babar_diff -> GetYaxis() -> SetTitleFont(132);
  gf_babar_diff -> GetYaxis() -> SetRangeUser(-400., 400.); 
  
  //gf_babar_diff -> Draw("ap");
  gf_cmd2_diff -> Draw("ap");
  //gf_bes_diff -> Draw("p");
  //gf_snd_diff -> Draw("p");
  
  TLegend *legd_cv_diff = new TLegend(0.2, 0.7, 0.9, 0.9);
  
  legd_cv_diff -> SetTextFont(132);
  legd_cv_diff -> SetFillStyle(0);
  legd_cv_diff -> SetBorderSize(0);
  legd_cv_diff -> SetNColumns(4);

  // legd_cv_diff -> AddEntry(gf_kloe2, "KLOE-2", "lep");
  legd_cv_diff -> AddEntry(gf_snd_diff, "SND", "lep");
  legd_cv_diff -> AddEntry(gf_cmd2_diff, "CMD-2", "lep");
  legd_cv_diff -> AddEntry(gf_babar_diff, "BaBar", "lep");
  legd_cv_diff -> AddEntry(gf_bes_diff, "BESIII", "lep");
  
  legd_cv_diff -> Draw("Same");
  
  legtextsize(legd_cv_diff, 0.04);
  */
  
  // save
  cv -> SaveAs("./plotOutput/crx3piCompr_W.pdf");
  //cv_diff -> SaveAs("./plotOutput/crx3piDiff_W.pdf");
 

  return 0;
  
}

