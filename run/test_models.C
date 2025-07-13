#include "../header/sm_para.h"
#include "../header/plot.h"
#include "../header/graph.h"
#include "../header/psf_fun.h"
#include "../header/amp_fun.h"
#include "../header/model.h"
#include "../header/test_models.h"

int test_models(){
  gErrorIgnoreLevel = kError;
  TGaxis::SetMaxDigits(4);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetFitFormat("6.4g");

  // compare bw and vmd models

  double sigma_max = 1606.53; //nb
  double BB = GetBB_new(sigma_max, mass_omega); //sigma_max * 1e-9 * (mass_omega * mass_omega) / 12. / pi;

  cout << "mass_omega = " << mass_omega << ", Gam_omega = " << Gam_omega << ", BB = " << BB << endl;

  double mass_min = 740.;
  double mass_max = 820.;
  
  TF1 *f_bw = new TF1("f_bw", Func_bw, mass_min, mass_max, 3);
  f_bw -> SetParameters(mass_omega, Gam_omega, BB);
  cout << f_bw -> Eval(mass_omega) << endl;
  
  TF1 *f_vmd = new TF1("f_vmd", vmd_crx3pi_fcn, mass_min, mass_max, 3);
  f_vmd -> SetParameters(mass_omega, Gam_omega, BB);
  cout << f_vmd -> Eval(mass_omega) << endl;
  
  const int np = 200;
  double dm = (mass_max - mass_min) / np;
  double m3pi = mass_min;
  double crx3pi_bw = 0.;
  double crx3pi_vmd = 0.;
  
  //2[i]=a[2]->eval(x[i]); 
  cout << "mass range = [" << mass_min << ", " << mass_max << "], np = " << np << ", dm = " << dm << endl;

  for (Int_t i = 0; i < np; i++) {

    M3PI[i] = m3pi;
    M3PI_ERR[i] = 0.;
    
    CRX3PI_BW[i] = f_bw -> Eval(m3pi);
    CRX3PI_BW_ERR[i] = 0.;

    CRX3PI_VMD[i] = f_vmd -> Eval(m3pi);
    CRX3PI_VMD_ERR[i] = 0.;
    
    m3pi = m3pi + dm;

    //cout << i << ", " << M3PI[i] << ", crx3pi_bw = " << CRX3PI_BW[i] << endl;
  }

  gf_crx3pi_bw = get_graph_syst(M3PI, CRX3PI_BW, M3PI_ERR, CRX3PI_BW_ERR, np);
  gf_crx3pi_bw -> SetMarkerSize(.8);
  gf_crx3pi_bw -> GetYaxis() -> SetTitleOffset(1);
  gf_crx3pi_bw -> GetYaxis() -> SetTitle("#sigma_{3#pi}");
  gf_crx3pi_bw -> GetXaxis() -> SetTitle("M_{3#pi} [MeV/c^{2}]");
  gf_crx3pi_bw -> SetName("gf_crx3pi_bw");
  
  gf_crx3pi_vmd = get_graph_syst(M3PI, CRX3PI_VMD, M3PI_ERR, CRX3PI_VMD_ERR, np);
  gf_crx3pi_vmd -> SetMarkerSize(.8);
  gf_crx3pi_vmd -> SetMarkerColor(kBlue);
  gf_crx3pi_vmd -> SetMarkerStyle(22);
  gf_crx3pi_vmd -> SetTitle("crx3pi vmd");
  gf_crx3pi_vmd -> SetName("gf_crx3pi_vmd");

  TCanvas *cv = new TCanvas("cv", "crx3pi bw vs. vmd", 0, 0, 1000, 800);
  cv -> SetLeftMargin(0.16);
  cv -> SetBottomMargin(0.15);

  gf_crx3pi_bw -> Draw("AP");
  gf_crx3pi_vmd -> Draw("P");
  

  //
  TLegend * legd_cv_compr = new TLegend(0.7, 0.6, 0.9, 0.85);
  
  SetLegend(legd_cv_compr);
  legd_cv_compr -> SetTextSize(0.04);
  legd_cv_compr -> SetNColumns(1);
  
  legd_cv_compr -> AddEntry(gf_crx3pi_bw, "BW", "lep");
  legd_cv_compr -> AddEntry(gf_crx3pi_vmd,  "VMD", "lep");
  
  legd_cv_compr -> Draw("Same");
  
  
  return 0;

}
