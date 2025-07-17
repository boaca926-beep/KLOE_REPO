#include "../header/sm_para.h"
#include "../header/path.h"
#include "../header/psf_fun.h"
#include "../header/amp_fun.h"
#include "../header/model.h"
#include "../header/method.h"

int plot_vmd_bw(){

  //gROOT->SetBatch(kTRUE);  
  gErrorIgnoreLevel = kError;
  TGaxis::SetMaxDigits(3);
  gStyle->SetOptStat(1110);
  gStyle->SetOptTitle(0);
  gStyle->SetStatBorderSize(0);
  //gStyle->SetFitFormat("6.2g");
  //gROOT->SetStyle("Plain");
  //gROOT->ForceStyle();
  
  const double M3pi_min = 710., M3pi_max = 1050.; // M3pi range in MeV/c^2
  double sigma_max = 1606.53; //nb
  
  double BB = GetBB_new(sigma_max, mass_omega); //sigma_max * 1e-9 * (mass_omega * mass_omega) / 12. / pi;

  TF1 *fun_bw = new TF1("func_bw", Func_bw, M3pi_min, M3pi_max, 3);
  fun_bw -> SetParameters(mass_omega, Gam_omega, BB);
  fun_bw -> SetNpx(2000);
  fun_bw -> SetLineColor(kBlack);

  TF1 *fun_vmd = new TF1("func_vmd", vmd_crx3pi_fcn, M3pi_min, M3pi_max, 3);
  fun_vmd -> SetParameters(mass_omega, Gam_omega, BB);

  TCanvas *cv_bw_vmd = new TCanvas("cv_bw_vmd","BW vs. VMD", 0, 0, 1200, 600);

  fun_vmd -> SetNpx(2000);
  fun_vmd -> SetLineColor(kBlack);
  fun_vmd -> SetLineStyle(kDashed);

  //fun_bw -> SetLineStyle(kDashed);
  fun_bw -> GetYaxis() -> SetTitle("#sigma(e^{+}e^{-}#rightarrow#pi^{+}#pi^{-}#pi^{0}) [nb]");
  fun_bw -> GetYaxis() -> SetLabelSize(0.05);//0.03
  fun_bw -> GetYaxis() -> SetTitleOffset(0.8);
  fun_bw -> GetYaxis() -> SetTitleSize(0.07);
  fun_bw -> GetYaxis() -> CenterTitle();

  // on X-axis
  fun_bw -> GetXaxis() -> SetTitle("M_{3#pi} [MeV/c^{2}]");
  fun_bw -> GetXaxis() -> SetLabelSize(0.05);//0.03
  fun_bw -> GetXaxis() -> SetTitleSize(0.07);
  fun_bw -> GetXaxis() -> SetTitleOffset(.9);
  fun_bw -> GetXaxis() -> CenterTitle();

  // general attri. of canvas
  cv_bw_vmd -> SetBottomMargin(0.15);//0.007
  cv_bw_vmd -> SetLeftMargin(0.13);

  fun_bw -> Draw(); 
  fun_vmd -> Draw("same");

  TLegend *legd_cv = new TLegend(0.6, 0.7, 0.9, 0.9);
  
  legd_cv -> SetTextFont(132);
  legd_cv -> SetFillStyle(0);
  legd_cv -> SetBorderSize(0);
  //legd_cv -> SetNColumns(4);
  
  legd_cv -> AddEntry(fun_vmd, "VMD", "l");
  legd_cv -> AddEntry(fun_bw, "BW", "l");
  legd_cv -> Draw("Same");

  cv_bw_vmd -> SaveAs("./cv_bw_vmd.pdf");
  
  return 0;
  
}
