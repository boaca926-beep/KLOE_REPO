#include "../header/plot.h"

const TString hist_nm = "_corred";
const TString infile_nm = "../../../analysis/crx3pi/output_norm/crx3pi0.root";

int smearmatr() {

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
    
  // print smearing matrix
  double p_value = 0.;
  double m3pi_lower = 0., m3pi_upper = 0.;

  cout << "Total number of bins = " << hsmearmatr -> ProjectionX() -> GetNbinsX() << "\n";
  TH1D* hproject_true = hsmearmatr -> ProjectionX();
  
  for (int j = 1; j <= hsmearmatr -> ProjectionX() -> GetNbinsX(); j ++ ) {

    double p_sum = 0.;
    
    for (int i = 1; i <= hsmearmatr -> ProjectionY() -> GetNbinsX(); i ++ ) {

      p_value = hsmearmatr -> GetBinContent(j, i);
      p_sum += p_value;

      if (p_value != 0.) {
	//cout << "(j, i) = (" << j << ", " << i << "), p_value = " << p_value << endl;
      }
  
    }

    m3pi_lower = hproject_true -> GetBinLowEdge(j);
    m3pi_upper = hproject_true -> GetBinLowEdge(j + 1);
    
    //cout << "bin_indx = " << j << ", m3pi = " << hproject_true -> GetBinCenter(j) << ", mass range [MeV/c²] = [" << m3pi_lower << ", " << m3pi_upper << "], mass interval = " << m3pi_upper - m3pi_lower << ", p_sum = " << p_sum << "\n\n";

  }

  //
  double entries = h2d_scatter_corred -> GetEntries();

  // plot
  //TCanvas * cv1 = new TCanvas("cv1", "", 0, 0, 1000, 800);
  TH1D* hprobsum_norm = (TH1D*) hprobsum -> Clone();
  
  hprobsum_norm -> Scale(0.01);
  //hprobsum_norm -> Draw("Hist");
  //gPad -> SetLogy();
  
  TCanvas * cv = new TCanvas("cv", "smearing matrix (signal sfw corrected)", 0, 0, 1300, 800);

  TPad *p2 = new TPad("p2", "p2", 0., 0., 1., 0.25);
  p2 -> Draw();
  p2 -> SetBottomMargin(0.4);
  p2 -> SetLeftMargin(0.15);
  p2 -> SetRightMargin(0.15);
  p2 -> SetGrid();
  
  TPad *p1 = new TPad("p1", "p1", 0., 0.24, 1., 1.);
  p1 -> Draw();
  p1 -> SetBottomMargin(0.025);//0.007
  p1 -> SetLeftMargin(0.15);
  p1 -> SetRightMargin(0.15);
  p1 -> cd();

  char display1[50];

  TPaveText *pt1 = new TPaveText(0.15, 0.82, 0.25, 0.83, "NDC");
  
  pt1 -> SetTextSize(0.05);
  pt1 -> SetFillColor(0);
  pt1 -> SetTextAlign(12);

  //sprintf(display, test);
  sprintf(display1,"Entries %0.0f ", entries);
  pt1 -> AddText(display1);

  const double xrange1 = 500., xrange2 = 900.; // mc true
  const double yrange1 = 760., yrange2 = 810.; // mc recon
  
  TLine *line = new TLine(760., xrange1, 760., xrange2);
  line -> SetLineColor(kBlack);
  line -> SetLineWidth(2);

  TLine *line1 = new TLine(800., xrange1, 800., xrange2);
  line1 -> SetLineColor(kBlack);
  line1 -> SetLineWidth(2);

  TLine *linex = new TLine(500., 800., 900., 800.); // horiz
  linex -> SetLineColor(kBlack);
  linex -> SetLineWidth(3);

  TLine *liney = new TLine(500., 760., 900., 760.); // horiz
  liney -> SetLineColor(kBlack);
  liney -> SetLineWidth(3);

  hsmearmatr -> GetXaxis() -> SetTitle(TString::Format("M^{true}_{3#pi} Events/%0.2f", m3pi_upper - m3pi_lower) + " [MeV/c^{2}]"); //SetTitle("M^{TRUE}_{3#pi} [MeV/c^{2}]");
  hsmearmatr -> GetYaxis() -> SetTitleSize(0.07);
  
  hsmearmatr -> GetYaxis() -> SetTitle(TString::Format("M^{smear}_{3#pi} Events/[%0.2f", m3pi_upper - m3pi_lower) + " MeV/c^{2}]"); //-> SetTitle("M^{REC}_{3#pi} [MeV/c^{2}]");

  hsmearmatr -> GetXaxis() -> SetTitleOffset(1.2);
  hsmearmatr -> GetYaxis() -> SetTitleOffset(1.);

  hsmearmatr -> GetXaxis() -> SetRangeUser(xrange1, xrange2);
  hsmearmatr -> GetYaxis() -> SetRangeUser(xrange1, xrange2);

  hsmearmatr -> GetXaxis() -> SetLabelSize(0.03); //20, 0.03
  hsmearmatr -> GetYaxis() -> SetLabelSize(0.05);
  hsmearmatr -> GetZaxis() -> SetLabelSize(0.05);

  hsmearmatr -> GetXaxis() -> CenterTitle();
  hsmearmatr -> GetYaxis() -> CenterTitle();

  hsmearmatr -> GetXaxis() -> SetLabelOffset(0.1);
  
  hsmearmatr -> SetStats(0);  

  //hsmearmatr -> Draw("TEXT0COL");
  hsmearmatr -> Draw("COLZ");
  //pt1 -> Draw("Same");

  // draw vertical lines indicating the omega region
  linex -> Draw("Same");
  liney -> Draw("Same");
  
  gPad -> SetLogz();

  p2 -> cd();

  // plot hprobsum_norm
  hprobsum_norm -> SetMarkerStyle(21);
  hprobsum_norm -> SetMarkerSize(0.5);
  
  hprobsum_norm -> GetYaxis() -> SetTitle("Probability");
  //hprobsum_norm -> GetXaxis() -> SetTitle("M^{TRUE}_{3#pi} [MeV/c^{2}]");
  hprobsum_norm -> GetXaxis() -> SetTitle(TString::Format("M^{true}_{3#pi} Events/[%0.2f", m3pi_upper - m3pi_lower) + " MeV/c^{2}]");  
  
  hprobsum_norm -> GetYaxis() -> CenterTitle();
  hprobsum_norm -> GetXaxis() -> CenterTitle();

  //hprobsum_norm -> GetYaxis() -> SetRangeUser(50., 150.);
  hprobsum_norm -> GetXaxis() -> SetRangeUser(xrange1, xrange2);

  hprobsum_norm -> GetYaxis() -> SetNdivisions(502);
  
  hprobsum_norm -> GetYaxis() -> SetTitleSize(0.14);
  hprobsum_norm -> GetXaxis() -> SetTitleSize(0.2);

  hprobsum_norm -> GetYaxis() -> SetLabelSize(0.15);
  hprobsum_norm -> GetXaxis() -> SetLabelSize(0.15);
  
  //hprobsum_norm -> GetYaxis() -> SetTitleFont(43);
  //hprobsum_norm -> GetXaxis() -> SetLabelFont(43); // Absolute font size in pixel (precision 3)
  
  hprobsum_norm -> GetYaxis() -> SetTitleOffset(.5);
  hprobsum_norm -> GetXaxis() -> SetTitleOffset(.86);

  hprobsum_norm -> GetXaxis() -> SetLabelOffset(0.02);
  //hprobsum_norm -> GetYaxis() -> SetLabelOffset(0.01);
  
  //hprobsum_norm -> GetYaxis() -> SetLabelFont(43); // Absolute font size in pixel (precision 3)

  hprobsum_norm -> SetStats(0);

  hprobsum_norm -> Draw("E0");
  //gPad -> SetLogy();
  
  format_h(hprobsum_norm, 4, 2);
  
  // get prob_avrge
  const double sigma_nb = 1.5;
  const double xmax = hprobsum -> GetXaxis() -> GetXmax(); //cout<<xmax<<endl;
  const double xmin = hprobsum -> GetXaxis() -> GetXmin();
  //double prob_theo = (1- 2 * (1. - TMath::Freq(sigma_nb - 1))) * 100.;
  double prob_theo = TMath::Erf(sigma_nb /TMath::Sqrt(2)) * 100.;
 
  double counter = 0., prob_sum_tmp = 0. , prob_sumsum = 0., prob_sum_avrge = 0.;

  for (int i = 1; i <= hprobsum -> GetNbinsX(); i ++) {

    prob_sum_tmp = hprobsum -> GetBinContent(i);
    prob_sumsum = 0.;

    m3pi_lower = hprobsum -> GetBinLowEdge(i);
    m3pi_upper = hprobsum -> GetBinLowEdge(i + 1);
    
    if (prob_sum_tmp != 0.){

      counter ++;

      prob_sumsum += prob_sum_tmp;
      
      prob_sum_avrge = prob_sumsum / counter;

      //cout << "bin_indx = " << i << ", mass range [MeV/c²] = (" << m3pi_lower << ", " << m3pi_upper << "), prob_sum = " << prob_sum_tmp << ", average = " << prob_sum_avrge << endl;
    
    }
    
  }//

  cout << "sigma band = " << sigma_nb << " sigma interval = " << prob_theo << "[%] \n"
       << "prob_sum_avrge = " << prob_sum_avrge << "\n";

  TLine *line2 = new TLine(xrange1, prob_theo * 1e-2, xrange2, prob_theo * 1e-2);
  line2 -> SetLineColor(kRed);
  line2 -> SetLineWidth(2);
  
  line2 -> Draw("Same");

  // save
  cv -> SaveAs("smearmatr_corred.pdf");
  
  return 0;

}
