//#include "para.h"
#include "../header/sm_para.h"
#include "../header/graph.h"
#include "../header/plot.h"
#include "../header/method.h"

const int list_size = 1000;
double m3pi_list[list_size], isrlumi_list[list_size], W0_list[list_size];
double isrlumi_apprx_list[list_size];
double isrlumi_diff[list_size];

//const double sqrtS = M_phi *1e-3;
//const double sample_size = 8373;
//const double Lumi_int = 1724470 * (sample_size / 8373);

int isrlumi() {

  gErrorIgnoreLevel = kError;
  TGaxis::SetMaxDigits(4);
  //gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);  

  //const TString infile = "/home/bo/Desktop/analysis/crx3pi/output_norm/crx3pi0.root";
  const TString infile = "/home/bo/Desktop/analysis_root_v6/input_chain_TDATA/omega_fit/omega_fit.root";
   
  cout << "\nCreate histograms from tree " << infile << "\n";

  TFile *intree = new TFile(infile);

  //getObj(intree);

  TIter next_tree(intree -> GetListOfKeys());

  TString objnm_tree, classnm_tree;
  
  int i = 0;

  const double fit_min = 760., fit_max = 800.;
  
  TKey *key;
  
  while ( (key = (TKey *) next_tree() ) ) {// start tree while lop
    
    i ++;
    
    objnm_tree   =  key -> GetName();
    classnm_tree = key -> GetClassName();
    key -> GetSeekKey();
    
    cout << "tree" << i << ": classnm = " << classnm_tree << ", objnm = " << objnm_tree << endl;
    
  }

  TTree* TRESULT = (TTree*)intree -> Get("TRESULT");

  //double isrlumi = 0., m3pi = 0., W0_full = 0.;
  int counter = 0;
  for (Int_t irow = 0; irow < TRESULT -> GetEntries(); irow++) {// start fill histos

      TRESULT -> GetEntry(irow); 

      isrlumi_list[irow] = TRESULT -> GetLeaf("Br_isrlumi") -> GetValue(0);

      isrlumi_apprx_list[irow] = TRESULT -> GetLeaf("Br_isrlumi_apprx") -> GetValue(0);

      isrlumi_diff[irow] = TMath::Abs((isrlumi_apprx_list[irow] - isrlumi_list[irow]) / isrlumi_list[irow]) * 1e2; 
		   
      if (isrlumi_list[irow] > 0.) {
	
	m3pi_list[irow] = TRESULT -> GetLeaf("Br_m3pi") -> GetValue(0);
	W0_list[irow] = TRESULT -> GetLeaf("Br_W0_full") -> GetValue(0);
	
	counter ++;
   
      }

      if (m3pi_list[irow] <= fit_max && m3pi_list[irow] >= fit_min) {

	cout << "m3pi = " << m3pi_list[irow] << ", ISR Lumi. Residual = " << isrlumi_diff[irow] << " [%]" << endl;

      }
      
      //cout << irow << ", m3pi = " << m3pi_list[irow] << ", W0_full = " << W0_list[irow] << ", isrlumi (apprx) = " << isrlumi_list[irow] << " (" << isrlumi_apprx_list[irow] << ")" << endl;

  }

  // plot graphs
  const int gf_size = counter;
  const double Delta_m3pi = m3pi_list[1] - m3pi_list[0];
  
  cout << "gf_size = " << gf_size << endl;
  //TGraph *gf_chi2, *gf_chi2_sub, *gf_smearing_fact, *gf_chi2_min;
  TGraph *gf_isrlumi, *gf_W0, *gf_isrlumi_apprx, *gf_isrlumi_diff;

  gf_isrlumi  = new TGraph(gf_size, m3pi_list, isrlumi_list);
  gf_isrlumi_apprx  = new TGraph(gf_size, m3pi_list, isrlumi_apprx_list);
  gf_W0  = new TGraph(gf_size, m3pi_list, W0_list);
  gf_isrlumi_diff = new TGraph(gf_size, m3pi_list, isrlumi_diff);
    
  //SetGFAttr(gf_isrlumi_apprx, "M_{3#pi} [MeV/c^{2}]", "ISR Luminosity [nb^{-1}]");
  SetGFAttr(gf_isrlumi_apprx, "M_{3#pi} [MeV/c^{2}]", "L_{ISR}^{apprx} [nb^{-1}]");
  SetGFAttr(gf_isrlumi, "M_{3#pi} [MeV/c^{2}]", "ISR Luminosity [nb^{-1}]");
  SetGFAttr(gf_isrlumi_diff, "M_{3#pi} [MeV/c^{2}]", "[%]");
  SetGFAttr(gf_W0, "M_{3#pi} [MeV/c^{2}]", "W_{0}(s,x)");

  //const double fit_min = 400., fit_max = 1000.;
  
  //gf_isrlumi -> GetXaxis() -> SetRangeUser(fit_min, fit_max);
  gf_isrlumi -> SetMarkerColor(kRed);
  gf_isrlumi -> SetMarkerSize(0.8);

  //gf_isrlumi -> GetYaxis() -> SetNdivisions(512);
  //gf_isrlumi -> GetXaxis() -> SetRangeUser(xrange0, xrange1);
  //gf_isrlumi -> GetYaxis() -> SetRangeUser(-50., h1d_IM3pi_TISR3PI_SIG -> GetMaximum() * 1.2);
  gf_isrlumi -> GetYaxis() -> SetTitleSize(0.07);
  gf_isrlumi -> GetYaxis() -> SetTitleOffset(0.7);
  gf_isrlumi -> GetYaxis() -> SetLabelSize(0.05);
  
  gf_isrlumi -> GetXaxis() -> SetLabelSize(0.03);
  //gf_isrlumi -> SetLineColor(kBlack);
  //gf_isrlumi -> GetYaxis() -> SetTitle(TString::Format("Events/[%0.2f MeV]", Delta_m3pi));
  //gf_isrlumi -> GetYaxis() -> SetTitle("#sigma(e^{+}e^{-}#rightarrow#pi^{+}#pi^{-}#pi^{0}) [nb]" );
  gf_isrlumi -> GetXaxis() -> SetTitleOffset(1.2);
  gf_isrlumi -> GetXaxis() -> SetLabelSize(0);
  
  char display1[50];
  char display2[50];
  char display3[50];
  char display4[50];
  
  TPaveText *pt1 = new TPaveText(0.12, 0.87, 0.2, 0.88, "NDC");
  TPaveText *pt2 = new TPaveText(0.12, 0.82, 0.2, 0.83, "NDC");
  TPaveText *pt3 = new TPaveText(0.12, 0.77, 0.2, 0.78, "NDC");

  TPaveText *pt4 = new TPaveText(0.12, 0.87, 0.2, 0.88, "NDC");
  TPaveText *pt5 = new TPaveText(0.12, 0.82, 0.2, 0.83, "NDC");
  
  pt1 -> SetTextSize(0.03);
  pt1 -> SetFillColor(0);
  pt1 -> SetTextAlign(12);

  pt2 -> SetTextSize(0.03);
  pt2 -> SetFillColor(0);
  pt2 -> SetTextAlign(12);

  pt3 -> SetTextSize(0.03);
  pt3 -> SetFillColor(0);
  pt3 -> SetTextAlign(12);

  pt4 -> SetTextSize(0.03);
  pt4 -> SetFillColor(0);
  pt4 -> SetTextAlign(12);

  pt5 -> SetTextSize(0.03);
  pt5 -> SetFillColor(0);
  pt5 -> SetTextAlign(12);

  //sprintf(display, test);
  sprintf(display1,"c.m. E=%0.2f [GeV]", sqrtS);
  sprintf(display2,"#DeltaM_{3#pi}=%0.2f [MeV/c^{2}]", Delta_m3pi);
  sprintf(display3, "%0.0f^{#circ} < #theta_{0} < %0.0f^{#circ}", 0., 180.);
  sprintf(display4,"L_{int}=%0.0f [nb^-1]", Lumi_int);
    
  pt1 -> AddText(display1); 
  pt2 -> AddText(display2);
  pt5 -> AddText(display3);
  pt4 -> AddText(display4);

  TCanvas *cv_isrlumi_apprx = new TCanvas("cv_isrlumi_apprx", "ISR Lumi. apprx", 0, 0, 1000, 700);

  cv_isrlumi_apprx -> SetBottomMargin(0.15);//0.007
  cv_isrlumi_apprx -> SetLeftMargin(0.13);

  gf_isrlumi_apprx -> SetMarkerSize(0.8);

  //gf_isrlumi_apprx -> GetYaxis() -> SetNdivisions(512);
  //gf_isrlumi_apprx -> GetXaxis() -> SetRangeUser(xrange0, xrange1);
  //gf_isrlumi_apprx -> GetYaxis() -> SetRangeUser(-50., h1d_IM3pi_TISR3PI_SIG -> GetMaximum() * 1.2);
  gf_isrlumi_apprx -> GetYaxis() -> SetTitleSize(0.07);
  gf_isrlumi_apprx -> GetYaxis() -> SetTitleOffset(0.8);
  gf_isrlumi_apprx -> GetYaxis() -> SetLabelSize(0.05);
  
  gf_isrlumi_apprx -> GetXaxis() -> SetLabelSize(0.05);
  gf_isrlumi_apprx -> GetXaxis() -> SetTitleSize(0.07);
  //gf_isrlumi_apprx -> SetLineColor(kBlack);
  //gf_isrlumi_apprx -> GetYaxis() -> SetTitle(TString::Format("Events/[%0.2f MeV]", Delta_m3pi));
  //gf_isrlumi_apprx -> GetYaxis() -> SetTitle("#sigma(e^{+}e^{-}#rightarrow#pi^{+}#pi^{-}#pi^{0}) [nb]" );
  gf_isrlumi_apprx -> GetXaxis() -> SetTitleOffset(0.9);
  
  gf_isrlumi_apprx -> Draw("AP");
  gPad -> SetLogy();
  
  //
  TCanvas *cv_isrlumi_compr = new TCanvas("cv_isrlumi", "ISR Lumi. comparison ", 0, 0, 1000, 700);

  TPad *p2 = new TPad("p2", "p2", 0., 0., 1., 0.25);
  p2 -> Draw();
  p2 -> SetBottomMargin(0.4);
  p2 -> SetLeftMargin(0.1);
  p2 -> SetGrid();
  
  TPad *p1 = new TPad("p1", "p1", 0., 0.25, 1., 1.);
  p1 -> Draw();
  p1 -> SetBottomMargin(0.02);//0.007
  p1 -> SetLeftMargin(0.1);
  p1 -> cd();

  gf_isrlumi -> Draw("AP");
  gf_isrlumi_apprx -> Draw("PSame");
  //pt2 -> Draw("Same");
  //pt4 -> Draw("Same");
  gPad -> SetLogy();
  
  TLegend *legd_cv = new TLegend(0.15, 0.7, 0.6, 0.8);
  
  legd_cv -> SetTextFont(132);
  legd_cv -> SetFillStyle(0);
  legd_cv -> SetBorderSize(0);
  legd_cv -> SetNColumns(2);

  legd_cv -> AddEntry(gf_isrlumi,  "L^{exact}_{ISR}", "lep");
  legd_cv -> AddEntry(gf_isrlumi_apprx, "L^{apprx}_{ISR}", "lep");  
  legd_cv -> Draw("Same");

  //TLegendEntry *header = (TLegendEntry*)legd_cv->GetListOfPrimitives()->First();
  //header->SetTextAlign(22);
  //header->SetTextColor(2);
  //header->SetTextSize(.09);
  
  legtextsize(legd_cv, 0.08);

  p2 -> cd();

  //gf_isrlumi_diff -> GetXaxis() -> SetRangeUser(fit_min, fit_max);
  gf_isrlumi_diff -> SetLineColor(0);
  gf_isrlumi_diff -> GetYaxis() -> SetNdivisions(505);
  gf_isrlumi_diff -> GetYaxis() -> SetTitleSize(20);
  gf_isrlumi_diff -> GetYaxis() -> SetTitleFont(43);
  gf_isrlumi_diff -> GetYaxis() -> SetLabelFont(43); // Absolute front size in pixel (precision 3)
  gf_isrlumi_diff -> GetYaxis() -> SetLabelSize(25);
  gf_isrlumi_diff -> GetYaxis() -> CenterTitle();
  gf_isrlumi_diff -> GetYaxis() -> SetTitle("Difference [%]");
  gf_isrlumi_diff -> GetYaxis() -> SetRangeUser(-1e-1, 1e-1);
  gf_isrlumi_diff -> GetYaxis() -> SetTitleOffset(2.5);
  //gf_isrlumi_diff -> GetYaxis() -> SetTitle("L^{apprx}_{ISR} - L^{exact}_{ISR})/L^{exact}_{ISR}");
  //gf_isrlumi_diff -> GetXaxis() -> SetLabelFont(43); // Absolute font size in pixel (precision 3)
  gf_isrlumi_diff -> GetXaxis() -> SetTitleOffset(0.9);
  //gf_isrlumi_diff -> GetXaxis() -> SetTitleSize(45);
  gf_isrlumi_diff -> GetXaxis() -> SetLabelSize(0.15);
  //gf_isrlumi_diff -> GetXaxis() -> SetLabelOffset(0.2);
  gf_isrlumi_diff -> GetXaxis() -> SetTitleSize(0.2);
  gf_isrlumi_diff -> GetXaxis() -> SetTitle("M_{3#pi} [MeV/c^{2}]");
  gf_isrlumi_diff -> GetXaxis() -> CenterTitle();
  gf_isrlumi_diff -> SetMarkerStyle(21);
  gf_isrlumi_diff -> SetMarkerSize(0.7);
  gf_isrlumi_diff -> Draw("E0");

  gf_isrlumi_diff -> Draw();
    
  TCanvas *cv_W0 = new TCanvas("cv_W0", " ", 0, 0, 1000, 700);

  cv_W0 -> SetBottomMargin(0.15);//0.007
  cv_W0 -> SetLeftMargin(0.13);

  gf_W0 -> SetMarkerSize(0.8);

  //gf_W0 -> GetYaxis() -> SetNdivisions(512);
  //gf_W0 -> GetXaxis() -> SetRangeUser(xrange0, xrange1);
  //gf_W0 -> GetYaxis() -> SetRangeUser(-50., h1d_IM3pi_TISR3PI_SIG -> GetMaximum() * 1.2);
  gf_W0 -> GetYaxis() -> SetTitleSize(0.07);
  gf_W0 -> GetYaxis() -> SetTitleOffset(0.8);
  gf_W0 -> GetYaxis() -> SetLabelSize(0.05);
  
  gf_W0 -> GetXaxis() -> SetLabelSize(0.05);
  gf_W0 -> GetXaxis() -> SetTitleSize(0.07);
  //gf_W0 -> SetLineColor(kBlack);
  //gf_W0 -> GetYaxis() -> SetTitle(TString::Format("Events/[%0.2f MeV]", Delta_m3pi));
  //gf_W0 -> GetYaxis() -> SetTitle("#sigma(e^{+}e^{-}#rightarrow#pi^{+}#pi^{-}#pi^{0}) [nb]" );
  gf_W0 -> GetXaxis() -> SetTitleOffset(0.9);
  
  gf_W0 -> Draw("AP");
  //pt1 -> Draw("Same");
  //pt3 -> Draw("Same");

  gPad -> SetLogy();
  
  cv_isrlumi_compr -> SaveAs("isrlumi_compr.pdf");
  cv_isrlumi_apprx -> SaveAs("isrlumi_apprx.pdf");
  cv_W0 -> SaveAs("W0.pdf");
  
  return 0;
  
}

