#include "../header/sm_para.h"
//#include "../header/hist.h"
#include "../header/plot.h"
//#include "../crx3pi/crx3pi_dg.C"
#include "crx3pi.h"
#include "compr.h"

int fit_result() {
  
  gErrorIgnoreLevel = kError;
  TGaxis::SetMaxDigits(4);
  //gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetFitFormat("6.4g");
  
  //output8373_finebinned
  
  const double sample_frac = Lumi_int / 1724470 * 100.;
  
  
  cout << "Getting input file from " << infile_tmp << "\n"
       << "Lumi. int [%] = " << sample_frac << "\n"
       << "input file = " << infile_tmp << "\n";
  
  TFile *intree = new TFile(infile_tmp);
  
  TIter next_tree(intree -> GetListOfKeys());
  
  TString objnm_tree, classnm_tree;
  
  int i = 0;
  TKey *key;
  
  while ( (key = (TKey *) next_tree() ) ) {// start tree while lop
    
    i ++;
    
    objnm_tree   =  key -> GetName();
    classnm_tree = key -> GetClassName();
    key -> GetSeekKey();
    
    //cout << " tree/histo" << i << ": classnm = " << classnm_tree << ", objnm = " << objnm_tree << endl;
    
  }

  // corrected signal sfw, get from ../crx3pi_sfw1d
  cout << "isr3pi_sfw1d = " << isr3pi_sfw1d << endl;
  //TH1D * h1d_IM3pi_TISR3PI_SIG_Scaled_corr = (TH1D*) h1d_IM3pi_TISR3PI_SIG_CORRED -> Clone();
  //TH1D * h1d_IM3pi_TISR3PI_SIG_Scaled_corr = (TH1D *) intree -> FindObject("h1d_IM3pi_TISR3PI_SIG_Scaled_corr") -> Clone();
  
  //format_h(h1d_IM3pi_TISR3PI_SIG_Scaled_corr, 4, 2);
      
  //h1d_IM3pi_TISR3PI_SIG_Scaled_corr -> Scale(isr3pi_sfw1d);
  //h1d_IM3pi_TISR3PI_SIG_Scaled_corr -> Draw();

  TH1D * h1d_IM3pi_TISR3PI_SIG = (TH1D *) intree -> Get("h1d_IM3pi_TISR3PI_SIG");
  
  
  // fill lists
  const int list_size = 1000;

  double M3PI_FIT[list_size], M3PI_FIT_ERR[list_size];
  double N3PI_OBS_FIT[list_size], N3PI_OBS_FIT_ERR[list_size];
  double N3PI_FIT[list_size], N3PI_FIT_ERR[list_size];
  double N3PI_DIFF[list_size], N3PI_DIFF_ERR[list_size];
  double ISRLUMI[list_size], CRX3PI_BW[list_size];
  double EFFICY[list_size], EFFICY_ERR[list_size];
  
  double Mass_omega_fit = 0., Mass_omega_err_fit = 0.;
  double Gam_omega_fit = 0., Gam_omega_err_fit = 0.;
  double BB_fit = 0., BB_err_fit = 0.;
  double chi2_sum = 0.;
  int ndf = 0;
  int bin_indx = 0;
  
  int fit_indx = 0;

  TTree *TCRX3PI = (TTree*)intree -> Get("TCRX3PI");
    
  for (Int_t irow = 0; irow < TCRX3PI -> GetEntries(); irow++) {// start fill
    
    TCRX3PI -> GetEntry(irow); //cout << irow << endl;

    //
    bin_indx = TCRX3PI -> GetLeaf("Br_bin_indx") -> GetValue(0);

    chi2_sum = TCRX3PI -> GetLeaf("Br_chi2_sum_crx3pi") -> GetValue(0);
    ndf = TCRX3PI -> GetLeaf("Br_ndf") -> GetValue(0);

    Mass_omega_fit = TCRX3PI -> GetLeaf("Br_Mass_omega_fit") -> GetValue(0);
    Mass_omega_err_fit = TCRX3PI -> GetLeaf("Br_Mass_omega_err_fit") -> GetValue(0);

    Gam_omega_fit = TCRX3PI -> GetLeaf("Br_Gam_omega_fit") -> GetValue(0);
    Gam_omega_err_fit = TCRX3PI -> GetLeaf("Br_Gam_omega_err_fit") -> GetValue(0);
    
    BB_fit = TCRX3PI -> GetLeaf("Br_BB_fit") -> GetValue(0);
    BB_err_fit = TCRX3PI -> GetLeaf("Br_BB_err_fit") -> GetValue(0);
    
    //
    M3PI_FIT[irow] = TCRX3PI -> GetLeaf("Br_M3PI") -> GetValue(bin_indx);
    M3PI_FIT_ERR[irow] = 0.;

    EFFICY[irow] = TCRX3PI -> GetLeaf("Br_EFFICY") -> GetValue(bin_indx);
    //EFFICY_ERR[irow] = TCRX3PI -> GetLeaf("Br_EFFICY_ERR") -> GetValue(bin_indx);

    N3PI_OBS_FIT[irow] = TCRX3PI -> GetLeaf("Br_N3PI_OBS") -> GetValue(bin_indx);
    N3PI_OBS_FIT_ERR[irow] = TCRX3PI -> GetLeaf("Br_N3PI_OBS_ERR") -> GetValue(bin_indx);

    N3PI_FIT[irow] = TCRX3PI -> GetLeaf("Br_N3PI_SMEAR") -> GetValue(bin_indx);
    N3PI_FIT_ERR[irow] = TCRX3PI -> GetLeaf("Br_N3PI_SMEAR_ERR") -> GetValue(bin_indx);

    N3PI_DIFF[irow] = TCRX3PI -> GetLeaf("Br_N3PI_DIFF") -> GetValue(bin_indx);
    N3PI_DIFF_ERR[irow] = TCRX3PI -> GetLeaf("Br_N3PI_DIFF_ERR") -> GetValue(bin_indx);
    //N3PI_DIFF[irow] = TCRX3PI -> GetLeaf("Br_N3PI_DIFF") -> GetValue(bin_indx) / N3PI_DIFF_ERR[irow];
    
    //
    ISRLUMI[irow] = TCRX3PI -> GetLeaf("Br_ISRLUMI") -> GetValue(bin_indx);
    CRX3PI_BW[irow] = TCRX3PI -> GetLeaf("Br_CRX3PI_BW") -> GetValue(bin_indx);
    
    fit_indx ++;

    cout << irow + 1 << ": n3pi obs. = " << N3PI_OBS_FIT[irow] << "+/-" << N3PI_OBS_FIT_ERR[irow] << ", n3pi diff. = " << N3PI_DIFF[irow] << "+/-" << N3PI_DIFF_ERR[irow] << endl;

    cout << "list indx = " << irow + 1 << "\n"
	 << "\tbin_indx = " << bin_indx << "\n"
	 << "\tm3pi = " << M3PI_FIT[irow] << " +/- " << M3PI_FIT_ERR[irow] << "\n"
	 << "\tn3pi_obs = " << N3PI_OBS_FIT[irow] << " +/- " << N3PI_OBS_FIT_ERR[irow] << "\n"
	 << "\tn3pi_fit = " << N3PI_FIT[irow] << " +/- " << N3PI_FIT_ERR[irow] << "\n"
	 << "\tn3pi_diff = " << N3PI_DIFF[irow] << " +/- " << N3PI_DIFF_ERR[irow] << "\n"
	 << "\tisr_lumi = " << ISRLUMI[irow] << "\n"
	 << "\tcrx3pi_bw = " << CRX3PI_BW[irow] << "\n"
	 << "\tefficy = " << EFFICY[irow] * 100. << " +/- " << EFFICY_ERR[irow] * 100. << " [%]\n"
	 << "\chi2_sum = " << chi2_sum << "\n\n";
    
    
  }

  //intree -> Close();
  double prob = TMath::Prob(chi2_sum, ndf);
  
  cout << "chi2_sum = " << chi2_sum << ", ndf = " << ndf << ", chi2_sum / ndf = " << chi2_sum / ndf << "\n"
       << "fit_indx = " << fit_indx << "\n"
       << "prob = " << prob << "\n";

  // check chi2-sum
  chi2_sum = 0.;
  double chi2_tmp = 0.;
  double diff_max = 0., diff_max_indx = 0;
  
  for (int i = 0; i < fit_indx; i ++) {

    if (chi2_tmp > diff_max) {
      diff_max = chi2_tmp;
      diff_max_indx = i;
    }
    
    chi2_tmp = ((N3PI_OBS_FIT[i] - N3PI_FIT[i]) / N3PI_OBS_FIT_ERR[i]) * ((N3PI_OBS_FIT[i] - N3PI_FIT[i]) / N3PI_OBS_FIT_ERR[i]);
    chi2_sum += chi2_tmp;

  }

  //chi2_sum = chi2_sum - diff_max;
  //cout << "diff_max = " << diff_max << ", n3pi_obs = " << N3PI_OBS_FIT[diff_max_indx] << endl;

  // Scale the signal histo
  double isr3pi_sfw = 0.;
  
  for (Int_t irow = 0; irow < TRESULT -> GetEntries(); irow++) {
    
    TRESULT -> GetEntry(irow); //cout << irow << endl;

    isr3pi_sfw = TRESULT -> GetLeaf("Br_isr3pi_sfw") -> GetValue(0);
    
  }

  h1d_IM3pi_TISR3PI_SIG -> Scale(isr3pi_sfw);
  

  cout << isr3pi_sfw << endl;

  // graphs
  const double Delta_m3pi = M3PI_FIT[1] - M3PI_FIT[0];
  //TGraph *gf_n3pi_diff = new TGraph (fit_indx, M3PI_FIT, N3PI_DIFF);
  TGraphErrors * gf_n3pi_diff = new TGraphErrors(fit_indx, M3PI_FIT, N3PI_DIFF, M3PI_FIT_ERR, N3PI_DIFF_ERR);
  SetGFAttr(gf_n3pi_diff, "M_{3#pi} [MeV/c^{2}]", "Residual");
  
  //
  //cout << N3PI_OBS_FIT_ERR[5] << endl;
  TGraphErrors * gf_n3pi_obs = new TGraphErrors(fit_indx, M3PI_FIT, N3PI_OBS_FIT, M3PI_FIT_ERR, N3PI_OBS_FIT_ERR);
  SetGFAttr(gf_n3pi_obs, "M_{3#pi} [MeV/c^{2}]", " ");
  gf_n3pi_obs -> SetMarkerColor(kBlack);
  gf_n3pi_obs -> SetMarkerSize(0.8);

  gf_n3pi_obs -> GetYaxis() -> SetNdivisions(512);
  //gf_n3pi_obs -> GetXaxis() -> SetRangeUser(xrange0, xrange1);
  gf_n3pi_obs -> GetYaxis() -> SetRangeUser(-50., h1d_IM3pi_TISR3PI_SIG -> GetMaximum() * 1.2);
  gf_n3pi_obs -> GetYaxis() -> SetTitleSize(0.07);
  gf_n3pi_obs -> GetYaxis() -> SetTitleOffset(0.6);
  gf_n3pi_obs -> GetYaxis() -> SetLabelSize(0.05);
  
  gf_n3pi_obs -> GetXaxis() -> SetLabelSize(0.03);
  //gf_n3pi_obs -> SetLineColor(kBlack);
  gf_n3pi_obs -> GetYaxis() -> SetTitle(TString::Format("Events/[%0.2f MeV/c^{2}]", Delta_m3pi));
  //gf_n3pi_obs -> GetYaxis() -> SetTitle("#sigma(e^{+}e^{-}#rightarrow#pi^{+}#pi^{-}#pi^{0}) [nb]" );
  gf_n3pi_obs -> GetXaxis() -> SetTitleOffset(1.2);
  gf_n3pi_obs -> GetXaxis() -> SetLabelSize(0);
  
  //
  TGraph* gf_n3pi_fit = new TGraphErrors(fit_indx, M3PI_FIT, N3PI_FIT);
  SetGFAttr(gf_n3pi_fit, "M_{3#pi} [MeV/c^{2}]", " ");
  gf_n3pi_fit -> SetMarkerStyle(21);
  //gf_n3pi_fit -> SetMarkerColor(kBlue);
  gf_n3pi_fit -> SetLineColor(kBlack);
  gf_n3pi_fit -> SetLineWidth(2);
  
  // plot
  TCanvas *cv_n3pi_fit = new TCanvas("cv_n3pi_fit", "Nobs 3pi", 0, 0, 1000, 800);

  char display1[50], display2[50], display3[50], display4[50], display5[50];

  TPaveText *pt1 = new TPaveText(0.14, 0.68, 0.43, 0.70, "NDC");
  TPaveText *pt2 = new TPaveText(0.14, 0.61, 0.43, 0.63, "NDC");
  TPaveText *pt3 = new TPaveText(0.14, 0.54, 0.43, 0.57, "NDC");
  
  TPaveText *pt4 = new TPaveText(0.67, 0.66, 0.85, 0.68, "NDC");
  TPaveText *pt5 = new TPaveText(0.67, 0.59, 0.85, 0.61, "NDC");
  
  PteAttr(pt1); pt1 -> SetTextSize(0.05);
  PteAttr(pt2); pt2 -> SetTextSize(0.05);
  PteAttr(pt3); pt3 -> SetTextSize(0.05);
  PteAttr(pt4); pt4 -> SetTextSize(0.05);
  PteAttr(pt5); pt5 -> SetTextSize(0.05);

  TLatex t(0.1, 0.2,"M_{#omega}");
  
  //sprintf(display1,"#chi^{2}/ndf=%0.2f", chi2_sum / ndf);
  pt1 -> AddText(Form("#Gamma_{#omega}=%0.2f#pm%0.2f MeV", Gam_omega_fit, Gam_omega_err_fit)); 
  pt2 -> AddText(Form("B_{ee}B_{3#pi}=(%0.2f#pm%0.2f)#times10^{-5}", BB_fit * 1e5, BB_err_fit * 1e5));
  pt3 -> AddText(Form("M_{#omega}=%0.2f#pm%0.2f MeV/c^{2}", Mass_omega_fit, Mass_omega_err_fit)); 
  pt4 -> AddText(Form("#chi^{2}/ndf=%0.2f", chi2_sum / ndf)); 
  pt5 -> AddText(Form("L_{int}=%0.2f fb^{-1}", Lumi_int * 1e-6)); 
  
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

  
  gf_n3pi_obs -> Draw("PA");
  gf_n3pi_fit -> Draw("cSame"); //P
  //h1d_IM3pi_TISR3PI_SIG_Scaled_corr -> Draw("HistSame");
  pt1 -> Draw("Same");
  pt2 -> Draw("Same");
  pt3 -> Draw("Same");
  pt4 -> Draw("Same");
  pt5 -> Draw("Same");
  //t.SetNDC(kTRUE);
  //t.Draw("Same");
  
  
  
  //
  TLegend * legd_cv = new TLegend(0.15, 0.75, 0.5, 0.85);
  
  SetLegend(legd_cv);
  legd_cv -> SetNColumns(2);
  
  //
  //legd_cv -> AddEntry(gf_n3pi_obs, "N^{obs}_{3#pi}", "lep");
  legd_cv -> AddEntry(gf_n3pi_obs, "Data", "lep");
  legd_cv -> AddEntry(gf_n3pi_fit, model_type, "l");
  //legd_cv -> AddEntry(h1d_IM3pi_TISR3PI_SIG,  "#pi^{+}#pi^{-}#pi^{0}", "lep");
  
  legd_cv -> Draw("Same");
  
  legtextsize(legd_cv, 0.07);
  
  p2 -> cd();
  
  gf_n3pi_diff -> SetLineColor(0);
  gf_n3pi_diff -> GetYaxis() -> SetNdivisions(505);
  gf_n3pi_diff -> GetYaxis() -> SetTitleSize(33);
  gf_n3pi_diff -> GetYaxis() -> SetTitleFont(43);
  gf_n3pi_diff -> GetYaxis() -> SetTitleOffset(1.2);
  gf_n3pi_diff -> GetYaxis() -> SetLabelFont(43); // Absolute front size in pixel (precision 3)
  gf_n3pi_diff -> GetYaxis() -> SetLabelSize(25);
  gf_n3pi_diff -> GetYaxis() -> CenterTitle();
  gf_n3pi_diff -> GetYaxis() -> SetTitle("Residual");
  gf_n3pi_diff -> GetYaxis() -> SetRangeUser(-1000., 1000.);
  //gf_n3pi_diff -> GetYaxis() -> SetTitle("(N_{d} - N_{mc})/#sqrt{N_{d}}");
  //gf_n3pi_diff -> GetXaxis() -> SetLabelFont(43); // Absolute font size in pixel (precision 3)
  gf_n3pi_diff -> GetXaxis() -> SetTitleOffset(0.85);
  //gf_n3pi_diff -> GetXaxis() -> SetTitleSize(45);
  gf_n3pi_diff -> GetXaxis() -> SetLabelSize(0.15);
  //gf_n3pi_diff -> GetXaxis() -> SetLabelOffset(0.2);
  gf_n3pi_diff -> GetXaxis() -> SetTitleSize(0.2);
  gf_n3pi_diff -> GetXaxis() -> SetTitle("M_{3#pi} [MeV/c^{2}]");
  gf_n3pi_diff -> GetXaxis() -> CenterTitle();
  gf_n3pi_diff -> SetMarkerStyle(21);
  gf_n3pi_diff -> SetMarkerSize(0.7);
  //gf_n3pi_diff -> Draw("E0");

  //TMultiGraph *mg = new TMultiGraph();
  //mg -> SetTitle("Exclusion graphs");
  
  //mg -> Add(gf_n3pi_diff);
  gf_n3pi_diff -> SetFillColor(1);
  gf_n3pi_diff -> SetLineStyle(9);
  gf_n3pi_diff -> SetFillStyle(3001);
  
  gf_n3pi_diff -> Draw("a3");
  
  
  //mg -> Draw("PC3");

  TLine *line = new TLine(740., 0., 845., 0.);
  line -> SetLineColor(kRed);
  line -> SetLineWidth(2);

  //line -> Draw("Same");

  TFile *f_output = new TFile("./result_" + file_type + ".root", "recreate");

  TTree* TRESULT = new TTree("TRESULT", "recreate");

  TRESULT -> SetAutoSave(0);

  double Mass_V[2] = {Mass_omega_fit, Mass_omega_err_fit};
  double Gamma_V[2] = {Gam_omega_fit, Gam_omega_err_fit};
  double BB[2] = {BB_fit, BB_err_fit};
  
  //TRESULT -> Branch("Br_Mass", &Mass_V, "Br_Mass[2]/D");
  //TRESULT -> Branch("Br_Gamma", &Gamma_V, "Br_Gamma[2]/D");
  //TRESULT -> Branch("Br_BB", &BB, "Br_BB[2]/D");
  //TRESULT -> Branch("Br_chi2_sum", &chi2_sum, "Br_chi2_sum/D");
  //TRESULT -> Branch("Br_ndf", &ndf, "Br_ndf/I");
  //TRESULT -> Branch("Br_Delta_m3pi", &Delta_m3pi, "Br_Delta_m3pi/D");
  
  //TRESULT -> Fill();
  //TRESULT -> Write();
  
  gf_n3pi_obs -> SetName("n3pi_obs");
  gf_n3pi_obs -> Write();

  gf_n3pi_fit -> SetName("n3pi_fit");
  gf_n3pi_fit -> Write();

  gf_n3pi_diff -> SetName("n3pi_diff");
  gf_n3pi_diff -> Write();
    
  f_output -> Close();
  
  // save
  cv_n3pi_fit -> SaveAs("./plots/n3pi_fit_" + file_type + ".pdf");
  
  return 0;


}


