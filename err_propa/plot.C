#include "sm_para.h"
#include "method.h"
#include "plot.h"

int plot(){

  gErrorIgnoreLevel = kError;
  TGaxis::SetMaxDigits(4);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetFitFormat("6.4g");

  // Files
  TFile *f_omega = new TFile("omega_fit.root"); 
  //getObj(f_omega);

  // Histos
  TH1D *hsig_true = (TH1D *) f_omega -> Get("hsig_true");
  TH1D *hsig_gen = (TH1D *) f_omega -> Get("hsig_gen");
  TH1D *hefficy = (TH1D *) f_omega -> Get("hefficy");
  
  TH2D *hsmearmatr = (TH2D *) f_omega -> Get("hsmearmatr");
  TH2D *hsmearmatr_err = (TH2D *) f_omega -> Get("hsmearmatr_err");
  TH2D *hcorrmatrix = (TH2D *) f_omega -> Get("hcorrmatrix");

  
  // Projections
  TH1D* hsmearmatr_projX = hsmearmatr -> ProjectionX("hsmearmatr_projX");
  hsmearmatr_projX -> SetName("hsmearmatr_projX");

  TH1D* hsmearmatr_projY = hsmearmatr -> ProjectionY("hsmearmatr_projY");
  hsmearmatr_projY -> SetName("hsmearmatr_projY");

  TH1D* hcorrmatrix_projX = hcorrmatrix -> ProjectionX("hcorrmatrix_projX");
  hcorrmatrix_projX -> SetName("hcorrmatrix_projX");

  TH1D* hcorrmatrix_projY = hcorrmatrix -> ProjectionY("hcorrmatrix_projY");
  hcorrmatrix_projY -> SetName("hcorrmatrix_projY");

  // Make projection

  // loop over histos comparing with omega_fit.h fit function
  
  // Create canvas and draw
  TCanvas *c0 = new TCanvas("c0", "efficy, hsig_true, hsig_gen", 800, 800);

  c0 -> Divide(1,2);
  c0 -> cd(1);

  gPad -> SetRightMargin(0.12);
  gPad -> SetLeftMargin(0.12);

  formatfill_h(hsig_true, 2, 3001);
  format_h(hsig_gen, 4, 2);
  
  hsig_gen -> Draw("Hist");
  hsig_true -> Draw("HistSame");
  gPad -> SetLogy();

  TLegend * legd_c0_1 = new TLegend(0.15, 0.75, 0.5, 0.9);
  SetLegend(legd_c0_1);
  legd_c0_1 -> SetNColumns(1);
  legd_c0_1 -> AddEntry(hsig_gen, "hsig_gen", "l");
  legd_c0_1 -> AddEntry(hsig_true, "hsig_true", "f");
  legd_c0_1 -> Draw("same");
  legtextsize(legd_c0_1, 0.06);

  c0 -> cd(2);

  format_h(hefficy, 4, 2);
  
  gPad -> SetRightMargin(0.12);
  gPad -> SetLeftMargin(0.12);

  hefficy -> Draw("E0");

  TLegend * legd_c0_2 = new TLegend(0.15, 0.75, 0.5, 0.9);
  SetLegend(legd_c0_2);
  legd_c0_2 -> SetNColumns(1);
  legd_c0_2 -> AddEntry(hsig_gen, "hefficy", "lep");
  legd_c0_2 -> Draw();
  legtextsize(legd_c0_2, 0.06);

  //
  TCanvas *c1 = new TCanvas("c1", "matrices", 800, 800);

  // Draw with colors and bin content
  c1 -> Divide(2,3);
  c1 -> cd(1);

  gPad -> SetRightMargin(0.12);
  gPad -> SetLeftMargin(0.12);
  
  hsmearmatr -> Draw("COLZ TEXT");

  // Formatting options
  hsmearmatr -> SetMarkerSize(1.2);
  hsmearmatr -> SetMarkerColor(kBlack);
  hsmearmatr -> GetYaxis() -> CenterTitle();
  hsmearmatr -> GetYaxis() -> SetTitle("M^{rec}_{3#pi}");
  hsmearmatr -> GetXaxis() -> CenterTitle();
  hsmearmatr -> GetXaxis() -> SetTitle("M^{true}_{3#pi}");
  
  c1->cd(2);

  gPad -> SetRightMargin(0.12);
  gPad -> SetLeftMargin(0.12);
  
  hcorrmatrix -> Draw("COLZ TEXT");
  
  hcorrmatrix -> SetMarkerSize(1.2);
  hcorrmatrix -> SetMarkerColor(kBlack);
  hcorrmatrix -> GetYaxis() -> CenterTitle();
  hcorrmatrix -> GetYaxis() -> SetTitle("M^{rec}_{3#pi}");
  hcorrmatrix -> GetXaxis() -> CenterTitle();
  hcorrmatrix -> GetXaxis() -> SetTitle("M^{true}_{3#pi}");
  
  gStyle->SetPaintTextFormat(".1f"); // Integer format
  
  // Make sure text is visible regardless of bin color
  hsmearmatr->SetContour(100);
  gPad->Update();
  
  // Add statistics box if needed
  //gStyle->SetOptStat(1111);

  c1->cd(3);

  gPad -> SetRightMargin(0.12);
  gPad -> SetLeftMargin(0.12);

  format_h(hsmearmatr_projX, 4, 2);
  
  hsmearmatr_projX -> Draw("Hist");
  hsmearmatr_projX -> GetXaxis() -> CenterTitle();
  hsmearmatr_projX -> GetXaxis() -> SetTitle("M^{true}_{3#pi}");
  
  TLegend * legd1 = new TLegend(0.15, 0.75, 0.5, 0.9);
  SetLegend(legd1);
  legd1 -> SetNColumns(1);
  legd1 -> AddEntry(hsmearmatr_projX, "hsmearmatr_projX", "lep");
  legd1 -> Draw();
  legtextsize(legd1, 0.06);

  c1->cd(4);

  gPad -> SetRightMargin(0.12);
  gPad -> SetLeftMargin(0.12);

  format_h(hcorrmatrix_projX, 4, 2);
  
  hcorrmatrix_projX -> Draw("Hist");
  hsig_true -> Draw("HistSame");
  hcorrmatrix_projX -> GetXaxis() -> CenterTitle();
  hcorrmatrix_projX -> GetXaxis() -> SetTitle("M^{true}_{3#pi}");
  
  TLegend * legd2 = new TLegend(0.15, 0.75, 0.5, 0.9);
  SetLegend(legd2);
  legd2 -> SetNColumns(1);
  legd2 -> AddEntry(hcorrmatrix_projX, "hcorrmatrix_projX", "lep");
  legd2 -> AddEntry(hsig_true, "hsig_true", "f");
  legd2 -> Draw();
  legtextsize(legd2, 0.06);
  
  c1->cd(5);

  gPad -> SetRightMargin(0.12);
  gPad -> SetLeftMargin(0.12);

  format_h(hsmearmatr_projY, 4, 2);
  
  hsmearmatr_projY -> Draw("Hist");
  hsmearmatr_projY -> GetXaxis() -> CenterTitle();
  hsmearmatr_projY -> GetXaxis() -> SetTitle("M^{true}_{3#pi}");
  
  TLegend * legd3 = new TLegend(0.15, 0.75, 0.5, 0.9);
  SetLegend(legd3);
  legd3 -> SetNColumns(1);
  legd3 -> AddEntry(hsmearmatr_projY, "hsmearmatr_projY", "lep");
  legd3 -> Draw();
  legtextsize(legd3, 0.06);

  c1->cd(6);

  gPad -> SetRightMargin(0.12);
  gPad -> SetLeftMargin(0.12);

  format_h(hcorrmatrix_projY, 4, 2);
  
  hcorrmatrix_projY -> Draw("Hist");
  hcorrmatrix_projY -> GetXaxis() -> CenterTitle();
  hcorrmatrix_projY -> GetXaxis() -> SetTitle("M^{true}_{3#pi}");
  
  TLegend * legd4 = new TLegend(0.15, 0.75, 0.5, 0.9);
  SetLegend(legd4);
  legd4 -> SetNColumns(1);
  legd4 -> AddEntry(hcorrmatrix_projY, "hcorrmatrix_projY", "lep");
  legd4 -> Draw();
  legtextsize(legd4, 0.06);

  //Since the uncertainties on εiεi​ are independent, the variance of the sum is the sum of the variances of each term. For a product Ni×εi, where Ni is exact (no uncertainty), the variance is:
  //Var(Ni^εi)=Ni^2×Var(εi)=Ni^2σi^2
  
  const int binsize = hsmearmatr -> GetNbinsX();
  cout << "binsize = " << binsize << endl;

  int total_bins = 0;

  for (int j = 1; j <= binsize; j ++ ) {

    double nb_rec = 0.;
    double nb_smeared = 0., nb_smeared_err = 0.;
    double nb_check = 0., nb_check_err = 0.;
    double nb_true = 0., nb_true_err = 0.;
    double nb_gen = 0., nb_gen_err = 0.;
    double efficy = 0., efficy_err = 0.;
  
    double prob_sum = 0.;
    double prob_tmp = 0., prob_tmp_err = 0.;

    double N2sigma2 = 0.;
    double N2sigma2_sum = 0.;

    double sigma2_y = 0.;
    double sigma2_y_sum = 0.;
    
    for (int i = 1; i <= binsize; i ++ ) {

      total_bins ++;
      
      prob_tmp = hsmearmatr -> GetBinContent(i, j);
      prob_tmp_err = hsmearmatr_err -> GetBinContent(i, j);
      prob_sum += prob_tmp;

      nb_true = hsig_true -> GetBinContent(i);
      nb_true_err = hsig_true -> GetBinError(i);
      nb_smeared += nb_true * prob_tmp;

      nb_gen = hsig_gen -> GetBinContent(i);
      nb_gen_err = hsig_gen -> GetBinError(i);
      efficy = hefficy -> GetBinContent(i);
      efficy_err = hefficy -> GetBinError(i);
      nb_check += nb_gen * efficy * prob_tmp; 
    
      if (prob_tmp != 0.) {
	N2sigma2 = nb_true * nb_true * prob_tmp_err * prob_tmp_err + nb_true_err * nb_true_err * prob_tmp * prob_tmp;
	sigma2_y = get_sigma2_y(nb_gen, efficy, prob_tmp, nb_gen_err, efficy_err, prob_tmp_err);
      }
      else {
	N2sigma2 = 0.;
	sigma2_y = 0.;
      }
      
      N2sigma2_sum += N2sigma2;
      sigma2_y_sum += sigma2_y;
      
      // check nb_rec = prob_tmp * efficy * nb_gen == corr_tmp_norm * nb_gen 

      if (j == 6) {

	cout << "a = " << nb_gen << "+/-" << nb_gen_err << ", b = " << efficy << "+/-" << efficy_err << ", c = " << prob_tmp << "+/-" << prob_tmp_err << ", sigma2_y = " << sigma2_y << endl;
	//cout << "bin " << total_bins << ", (rec_indx, true_indx) = (" << j << ", " << i << "), nb_true = " << nb_true << "+/-" << nb_true_err << ", prob_tmp = " << prob_tmp << "+/-" << prob_tmp_err << ", N2sigma2 = " << N2sigma2 << ", N2sigma2_sum = " << N2sigma2_sum << "\n";
      }
      
    }

    nb_rec = hcorrmatrix_projY -> GetBinContent(j);
    nb_smeared_err = TMath::Sqrt(N2sigma2_sum);
    nb_check_err = TMath::Sqrt(sigma2_y_sum);

    if (j == 6) {
      cout << "rec. indx = " << j << ", prob_sum = " << prob_sum << ", nb_rec = " << nb_rec << ", nb_smeared (check) = " << nb_smeared << "+/-" << nb_smeared_err << " (" << nb_check << "+/-" << nb_check_err << ")\n\n";
    }
    
  }
  
    
  TFile *f_out = new TFile("plot.root", "recreate");
  hcorrmatrix_projY -> Write();
  hcorrmatrix_projX -> Write();
  
  return 0;
  
}
