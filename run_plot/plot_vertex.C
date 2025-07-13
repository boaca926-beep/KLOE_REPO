#include "../header/plot.h"

int plot_vertex(){

  gErrorIgnoreLevel = kError;
  TGaxis::SetMaxDigits(4);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetFitFormat("6.4g");

  const TString mc_type = "sig"; // sig ksl
  
  TString file_path = "../../" + mc_type + "_pre.root";
  
  TFile* file_input = new TFile(file_path);

  TIter next_tree(file_input -> GetListOfKeys());
  TString objnm_tree, classnm_tree;
  
  //int k = 0;
  TKey *key;
  
  while ( (key = (TKey *) next_tree() ) ) {// start cutlist_tree while loop
    
    objnm_tree   = key -> GetName();
    classnm_tree = key -> GetClassName();
    key -> GetSeekKey();
	
    cout << "classnm = " << classnm_tree << ", objnm = " << objnm_tree << endl;
    
  }

  TH1D *h_iv_ip = (TH1D *)file_input->Get("h_iv_ip");
  TH1D *h_iv_indx = (TH1D *)file_input->Get("h_iv_indx");

  TH1D *h_iv_ip_1 = (TH1D *)file_input->Get("h_iv_ip_1");
  TH1D *h_iv_indx_1 = (TH1D *)file_input->Get("h_iv_indx_1");

  TH2D *h_nv_ip = (TH2D *)file_input->Get("h_nv_ip");
  TH2D *h_nv_ip_1 = (TH2D *)file_input->Get("h_nv_ip_1");

  TH2D *h_xpca = (TH2D *)file_input->Get("h_xpca");
  TH2D *h_ypca = (TH2D *)file_input->Get("h_ypca");
  TH2D *h_zpca = (TH2D *)file_input->Get("h_zpca");

  // normalization
  h_nv_ip -> Scale(1./h_nv_ip->Integral(), "width");
  h_nv_ip_1 -> Scale(1./h_nv_ip_1->Integral(), "width");

  h_xpca -> Scale(1./h_xpca->Integral(), "width");
  h_ypca -> Scale(1./h_ypca->Integral(), "width");
  h_zpca -> Scale(1./h_zpca->Integral(), "width");
  
  
  // attribute
  h_iv_ip -> SetLineWidth(2);
  h_iv_indx -> SetLineWidth(2);
  
  h_iv_ip_1 -> SetFillStyle(3005);
  h_iv_ip_1 -> SetFillColor(kBlue);
  
  h_iv_indx_1 -> SetFillStyle(3005);
  h_iv_indx_1 -> SetFillColor(kBlue);
  
  TPaveText *pt_1 = new TPaveText(0.16, 0.82, 0.7, 0.86, "NDC");
  TString pt_1_text = mc_type + " (After trigger, filfo, streaming)"; 
  SetPte(pt_1, pt_1_text);

  TPaveText *pt_2 = new TPaveText(0.16, 0.82, 0.7, 0.86, "NDC");
  TString pt_2_text = mc_type + " (After track and photon selection)"; 
  SetPte(pt_2, pt_2_text);
  
  TCanvas *cv_2d = new TCanvas("cv_2d", "nv vs. iv_ip", 0, 0, 1600, 800);
  cv_2d -> Divide(2, 1);
  
  cv_2d -> cd(1);
  gPad -> SetBottomMargin(0.15);//0.007
  gPad -> SetLeftMargin(0.15);
  gPad -> SetRightMargin(0.17);
  
  h_nv_ip -> GetXaxis() -> SetNdivisions(5);
  h_nv_ip -> GetXaxis() -> SetTitle("n_{v}");
  h_nv_ip -> GetXaxis() -> SetTitleOffset(1.2);
  h_nv_ip -> GetXaxis() -> SetTitleSize(0.06);
  h_nv_ip -> GetXaxis() -> CenterTitle();
  h_nv_ip -> GetXaxis() -> SetLabelSize(0.06);
  h_nv_ip -> GetXaxis() -> SetLabelOffset(0.01);
  //h_nv_ip -> GetXaxis() -> SetRangeUser(0.2, 0.6);
  
  h_nv_ip -> GetYaxis() -> SetTitle("iv_ip");
  h_nv_ip -> GetYaxis() -> SetLabelOffset(0.01);
  h_nv_ip -> GetYaxis() -> SetTitleOffset(1.2);
  h_nv_ip -> GetYaxis() -> SetLabelSize(0.06);
  h_nv_ip -> GetYaxis() -> SetTitleSize(0.06);
  h_nv_ip -> GetYaxis() -> CenterTitle();

  h_nv_ip -> GetZaxis() -> SetLabelSize(0.06);

  h_nv_ip -> Draw("COLZ");
  pt_1 -> Draw("Same");
  gPad->SetLogz(1);

  cv_2d -> cd(2);
  gPad -> SetBottomMargin(0.15);//0.007
  gPad -> SetLeftMargin(0.15);
  gPad -> SetRightMargin(0.17);

  h_nv_ip_1 -> GetXaxis() -> SetNdivisions(5);
  h_nv_ip_1 -> GetXaxis() -> SetTitle("n_{v}");
  h_nv_ip_1 -> GetXaxis() -> SetTitleOffset(1.2);
  h_nv_ip_1 -> GetXaxis() -> SetTitleSize(0.06);
  h_nv_ip_1 -> GetXaxis() -> CenterTitle();
  h_nv_ip_1 -> GetXaxis() -> SetLabelSize(0.06);
  h_nv_ip_1 -> GetXaxis() -> SetLabelOffset(0.01);
  //h_nv_ip_1 -> GetXaxis() -> SetRangeUser(0.2, 0.6);
  
  h_nv_ip_1 -> GetYaxis() -> SetTitle("iv_ip");
  h_nv_ip_1 -> GetYaxis() -> SetLabelOffset(0.01);
  h_nv_ip_1 -> GetYaxis() -> SetTitleOffset(1.2);
  h_nv_ip_1 -> GetYaxis() -> SetLabelSize(0.06);
  h_nv_ip_1 -> GetYaxis() -> SetTitleSize(0.06);
  h_nv_ip_1 -> GetYaxis() -> CenterTitle();

  h_nv_ip_1 -> GetZaxis() -> SetLabelSize(0.06);

  h_nv_ip_1 -> Draw("COLZ");
  pt_2 -> Draw("Same");
  gPad->SetLogz(1);


  TCanvas *cv_pca = new TCanvas("cv_pca", "pca track1 vs. track2 " + mc_type, 0, 0, 1600, 800);
  cv_pca -> Divide(3, 1);
  
  cv_pca -> cd(1);
  gPad -> SetBottomMargin(0.15);//0.007
  gPad -> SetLeftMargin(0.15);
  gPad -> SetRightMargin(0.2);
  
  h_xpca -> GetXaxis() -> SetNdivisions(5);
  h_xpca -> GetXaxis() -> SetTitle("xpca[1]");
  h_xpca -> GetXaxis() -> SetTitleOffset(1.2);
  h_xpca -> GetXaxis() -> SetTitleSize(0.06);
  h_xpca -> GetXaxis() -> CenterTitle();
  h_xpca -> GetXaxis() -> SetLabelSize(0.06);
  h_xpca -> GetXaxis() -> SetLabelOffset(0.01);
  //h_xpca -> GetXaxis() -> SetRangeUser(0.2, 0.6);
  
  h_xpca -> GetYaxis() -> SetTitle("xpca[2]");
  h_xpca -> GetYaxis() -> SetLabelOffset(0.01);
  h_xpca -> GetYaxis() -> SetTitleOffset(1.2);
  h_xpca -> GetYaxis() -> SetLabelSize(0.06);
  h_xpca -> GetYaxis() -> SetTitleSize(0.06);
  h_xpca -> GetYaxis() -> CenterTitle();

  h_xpca -> GetZaxis() -> SetLabelSize(0.06);

  h_xpca -> Draw("COLZ");
  gPad->SetLogz(1);

  cv_pca -> cd(2);
  gPad -> SetBottomMargin(0.15);//0.007
  gPad -> SetLeftMargin(0.15);
  gPad -> SetRightMargin(0.2);

  h_ypca -> GetXaxis() -> SetNdivisions(5);
  h_ypca -> GetXaxis() -> SetTitle("ypca[1]");
  h_ypca -> GetXaxis() -> SetTitleOffset(1.2);
  h_ypca -> GetXaxis() -> SetTitleSize(0.06);
  h_ypca -> GetXaxis() -> CenterTitle();
  h_ypca -> GetXaxis() -> SetLabelSize(0.06);
  h_ypca -> GetXaxis() -> SetLabelOffset(0.01);
  //h_ypca -> GetXaxis() -> SetRangeUser(0.2, 0.6);
  
  h_ypca -> GetYaxis() -> SetTitle("ypca[2]");
  h_ypca -> GetYaxis() -> SetLabelOffset(0.01);
  h_ypca -> GetYaxis() -> SetTitleOffset(1.2);
  h_ypca -> GetYaxis() -> SetLabelSize(0.06);
  h_ypca -> GetYaxis() -> SetTitleSize(0.06);
  h_ypca -> GetYaxis() -> CenterTitle();

  h_ypca -> GetZaxis() -> SetLabelSize(0.06);

  h_ypca -> Draw("COLZ");
  gPad->SetLogz(1);

  cv_pca -> cd(3);
  gPad -> SetBottomMargin(0.15);//0.007
  gPad -> SetLeftMargin(0.15);
  gPad -> SetRightMargin(0.2);

  h_zpca -> GetXaxis() -> SetNdivisions(5);
  h_zpca -> GetXaxis() -> SetTitle("zpca[1]");
  h_zpca -> GetXaxis() -> SetTitleOffset(1.2);
  h_zpca -> GetXaxis() -> SetTitleSize(0.06);
  h_zpca -> GetXaxis() -> CenterTitle();
  h_zpca -> GetXaxis() -> SetLabelSize(0.06);
  h_zpca -> GetXaxis() -> SetLabelOffset(0.01);
  //h_zpca -> GetXaxis() -> SetRangeUser(0.2, 0.6);
  
  h_zpca -> GetYaxis() -> SetTitle("zpca[2]");
  h_zpca -> GetYaxis() -> SetLabelOffset(0.01);
  h_zpca -> GetYaxis() -> SetTitleOffset(1.2);
  h_zpca -> GetYaxis() -> SetLabelSize(0.06);
  h_zpca -> GetYaxis() -> SetTitleSize(0.06);
  h_zpca -> GetYaxis() -> CenterTitle();

  h_zpca -> GetZaxis() -> SetLabelSize(0.06);

  h_zpca -> Draw("COLZ");
  gPad->SetLogz(1);

  /*
  TCanvas *cv = new TCanvas("cv", " ", 0, 0, 1600, 800);
  cv -> Divide(2,1);
  cv -> SetBottomMargin(0.1);//0.007
  
  cv -> cd(1);
  gPad -> SetLeftMargin(0.15);

  h_iv_ip -> GetXaxis() -> CenterTitle();
  h_iv_ip -> GetXaxis() -> SetTitle("n_{v}");
  h_iv_ip -> GetXaxis() -> SetTitleSize(0.06);
  h_iv_ip -> GetXaxis() -> SetLabelSize(0.04);  
  h_iv_ip -> GetXaxis() -> SetTitleOffset(.7);
  
  h_iv_ip -> Draw();
  h_iv_ip_1 -> Draw("Same");
  pt -> Draw("Same");
  
  gPad->SetLogy(1);

  TLegend * legd_cv = new TLegend(0.45, 0.8, 0.6, 0.9);
  
  SetLegend(legd_cv);
  legd_cv -> SetNColumns(3);
  legd_cv -> AddEntry(h_iv_ip_1, "n_{v}=1", "f");
  legd_cv -> Draw("Same");
  legtextsize(legd_cv, 0.04);
  
  cv -> cd(2);
  gPad -> SetLeftMargin(0.15);

  h_iv_indx -> GetXaxis() -> CenterTitle();
  h_iv_indx -> GetXaxis() -> SetTitle("v_{indx}");
  h_iv_indx -> GetXaxis() -> SetTitleSize(0.06);
  h_iv_indx -> GetXaxis() -> SetLabelSize(0.04);  
  h_iv_indx -> GetXaxis() -> SetTitleOffset(.7);
  
  h_iv_indx -> Draw();
  h_iv_indx_1 -> Draw("Same");
  gPad->SetLogy(1);
  pt -> Draw("Same");

  TLegend * legd_cv1 = new TLegend(0.45, 0.8, 0.6, 0.9);
  
  SetLegend(legd_cv1);
  legd_cv1 -> SetNColumns(3);
  legd_cv1 -> AddEntry(h_iv_indx_1, "n_{v}=1", "f");
  legd_cv1 -> Draw("Same");
  legtextsize(legd_cv1, 0.04);
  */
  
  //cv -> SaveAs("../code_ref/vertex_" + mc_type + ".pdf");
  cv_2d -> SaveAs("../code_ref/vertex2D_" + mc_type + ".pdf");
  cv_pca -> SaveAs("../code_ref/pca_" + mc_type + ".pdf");
  
  return 0;
  
}
