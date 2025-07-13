#include "../header/assignment.h"

int assignment() {

  gErrorIgnoreLevel = kError;
  TGaxis::SetMaxDigits(3);
  gStyle->SetOptStat(1110);
  gStyle->SetOptTitle(0);
  gStyle->SetStatBorderSize(0);
  
  // draw kernel function as function of s
  double s_low = 4. * (mass_mu * 1e-3) * (mass_mu * 1e-3);
  double s_upp = 100.;
  const double s_mpi = TMath::Power(4 * mass_pi * 1e-3, 2);
  TF1 *kernel = new TF1("kernel", kernelFcn, s_low, s_upp, 1);
  kernel -> SetParameter(0, mass_mu);
  const int Np = 1000;
  kernel -> SetNpx(Np);
  
  // evaluate hat_kerneal function valule at 4m2pi^2
  double s_2pi = 4. * (mass_pi * 1e-3) * (mass_pi * 1e-3);// c.m energy of 2pi at rest [GeV/c^2]^2
  double kernalM2pi = kernel -> Eval(s_2pi);

  cout << "hat_kernal function at M2pi (GeV/c^2) = " << kernalM2pi << ", lower bound of sqrt(s) = 2Mmu: " << 2 * (mass_mu * 1e-3) << " GeV/c^2" << endl;

  const int np = 1;
  double grX[np] = {s_2pi};
  double grY[np] = {kernalM2pi};
  TGraph *gr = new TGraph(np, grX, grY);
  
  TCanvas *cv_kernel = new TCanvas("CVKernel", "kernel", 1000, 700);
  gPad -> SetBottomMargin(0.12);
  gPad -> SetLeftMargin(0.15);

  gPad -> SetLogx();

  gr -> SetMarkerSize(2);
  gr -> SetLineColor(2);
  gr -> SetLineWidth(4);
  gr -> SetMarkerColor(4);
  gr -> SetMarkerSize(1.5);
  gr -> SetMarkerStyle(21);

  gr -> GetXaxis() -> SetTitleSize(0.05);
  gr -> GetXaxis() -> SetLabelSize(0.04);
  //gr -> GetXaxis() -> SetNdivisions(512);
  gr -> GetXaxis() -> SetLimits(1e-2, s_upp);
  
  gr -> GetXaxis() -> CenterTitle();
  gr -> GetXaxis() -> SetTitle("s [GeV^{2}/c^{4}]");
  gr -> GetYaxis() -> SetTitleSize(0.05);
  gr -> GetYaxis() -> SetLabelSize(0.04);
  //gr -> GetYaxis() -> SetNdivisions(512);
  gr -> GetYaxis() -> CenterTitle();
  gr -> GetYaxis() -> SetTitle("#tilde{K}(s)");
  gr -> GetYaxis() -> SetTitleOffset(1.2);
  gr -> GetYaxis() -> SetNdivisions(10);
  gr -> GetYaxis() -> SetRangeUser(0.4, 1.);
  
  
  gr -> Draw("acp");
  kernel -> Draw("same");
  
  
  cv_kernel -> SetGrid();
  //cv_kernel -> Update();
  
  cv_kernel -> SaveAs("kernel.pdf");
	
  
  // gamma factor in g-2
  const double gamma = 29.3;
  double beta = TMath::Sqrt(1. - 1. / gamma / gamma);

  cout << "gamma = " << gamma << ", beta = " << beta << endl;
  
  
  // area if the normal distribution
  double sigma_size = 1.96; // given number of sigma
  double a = TMath::Erf(sigma_size /TMath::Sqrt(2)) * 100.;
 
  cout << "number of sigma = " << sigma_size << ", area [%] = " << a << endl;
  
  // Calculate probability, Li's Thesis p.82
  // Given chi2=356
  // Number of parameters include Normalization N, nb_para=8+1=9
  // Prob=0.58
  
  const double chi2_sum = 356., bins = 371, nb_para = 9.;
  const double ndf = bins - nb_para;

  double prob = TMath::Prob(chi2_sum, ndf);
    
  cout << "chi2_sum = " << chi2_sum << ", bins = " << bins << ", nb_para = " << nb_para << ", ndf = " << ndf << ", prob = " << prob << "\n";

  /// g-2 results
  // combine g-2 experimental results BNL and Fermi Lab, world average
  const double g2_fermi = 0.00116592040, g2_fermi_err = 0.00000000054; // PHYSICAL REVIEW LETTERS126,141801 (2021)
  const double g2_bnl = 0.00116592080, g2_bnl_stat_err = 0.00000000054, g2_bnl_syst_err = 0.00000000033; //PHYSICAL REVIEW D73,072003 (2006)
  const double g2_bnl_err = TMath::Sqrt(g2_bnl_stat_err * g2_bnl_stat_err + g2_bnl_syst_err * g2_bnl_syst_err);

  double weight_fermi = TMath::Power(1. / g2_fermi_err, 2);
  double weight_bnl = TMath::Power(1. / g2_bnl_err, 2);
  double weight_sum = weight_fermi + weight_bnl;
  //double weight_fermi_ratio = TMath::Power(1. / 5.4, 2) / (TMath::Power(1. / 5.4, 2) + TMath::Power(1. / 6.3, 2)); // get rid of 10^-11
  //double weight_bnl_ratio = TMath::Power(1. / 6.3, 2) / (TMath::Power(1. / 5.4, 2) + TMath::Power(1. / 6.3, 2)); // get rid of 10^-11
  //double g2_comb = g2_fermi * weight_fermi_ratio + g2_bnl * weight_bnl_ratio;
  double g2_comb = (g2_fermi * weight_fermi + g2_bnl * weight_bnl) / (weight_fermi + weight_bnl);
  
  //cout << "weight_fermi_ratio = " << weight_fermi_ratio << endl;
  //cout << "weight_bnl_ratio = " << weight_bnl_ratio << endl;
  printf("g2_comb = %12.10g\n", g2_comb);
  printf("g2_fermi=%12.10g +/- %12.10g\n", g2_fermi, g2_fermi_err);
  printf("g2_bnl = %12.8g +/- %12.8g\n", g2_bnl, g2_bnl_err);
  printf("1/sqrt(weight_sum) = %12.10g\n", 1. / TMath::Sqrt(weight_sum));

  // given theoretical prediction and the world average, calcualte corresponding sigma
  // The accepted theoretical values for the muon are: g-factor: 2.00233183620(86), anomalous magnetic moment: 0.00116591810(43)
  // The new experimental world-average results announced by the Muon g-2 collaboration today are: g-factor: 2.00233184122(82), anomalous magnetic moment: 0.00116592061(41)
  const double g2_theory = 0.00116591810, g2_theory_err = 0.00000000043; // SM prediction, T. Aoyama, N. Asmussen, M. Benayoun, J. Bijnens, T.Blumet al., The anomalous magnetic moment of the muonin the standard model,Phys. Rep.887, 1 (2020).
  const double g2_exp_avrge = 0.00116592061, g2_exp_avrge_err = 0.00000000041; // wold average

  double g2_diff = g2_exp_avrge - g2_theory;
  double g2_err = TMath::Sqrt(g2_exp_avrge_err * g2_exp_avrge_err + g2_theory_err * g2_theory_err);
  double g2_sigma = g2_diff / g2_err;
  double sigma_nb = 4.2;
  prob = GetPValue(sigma_nb);
  
  cout <<"\n";
  printf("g2_theory = %12.10g +/- %12.10g\n", g2_theory, g2_theory_err);
  printf("g2_exp_avrge = %12.10g +/- %12.10g\n", g2_exp_avrge, g2_exp_avrge_err);
  printf("g2_diff = %12.10g\n", g2_diff);
  printf("g2_err = %12.10g\n", g2_err);
  printf("g2_sigma = %12.10g\n", g2_sigma);
  printf("sigma_nb = %12.10g, p-value = %12.10g\n", sigma_nb, prob);

  double NB_SIGMA[5] = {1, 1.5, 2, 2.5, 2.73};
  double p_tmp = 0.;
  
  for (int i = 0; i < 5; i ++) {
    p_tmp = GetPValue(NB_SIGMA[i]);
    cout << NB_SIGMA[i] << " sigma, p-value = " << 1 - 2 * p_tmp << endl;
  }
  
  cout << "\ng2_theory = " << g2_theory << " +/- " << g2_theory_err << "\n"
       << "g2_exp_avrge = " << g2_exp_avrge << " +/- " << g2_exp_avrge_err << "\n"
       << "g2_diff = " << g2_diff << "\n"
       << "g2_err = " << g2_err << "\n"
       << "g2_sigma = " << g2_sigma << "\n";

  // Calculate the number of sigma comparing omega mass results from CMD-2 (2003) and CMD-2 (2004)
  const double momega_a = 783.20; // CMD-2 (2004)
  const double momega_a_err = TMath::Sqrt(0.13 * 0.13 + 0.16 * 0.16);

  const double momega_b = 782.68; // CMD-2 (2003) 
  const double momega_b_err = TMath::Sqrt(0.09 * 0.09 + 0.04 * 0.04);

  double momega_diff = momega_a - momega_b;
  double momega_diff_err = TMath::Sqrt(momega_a_err * momega_a_err + momega_b_err + momega_b_err);
  double momega_diff_sigmanb = momega_diff / momega_diff_err;
  
  cout << "momega_a = " << momega_a << "+/-" << momega_a_err << "\n"
       << "momega_b = " << momega_b << "+/-" << momega_b_err << "\n"
       << "momega_diff = " << momega_diff << "+/-" << momega_diff_err << "\n"
       << "momega_diff_sigmanb = " << momega_diff_sigmanb << endl;

  // calcuate significance for table 7.1
  // omega parameters: mass, width, BB
  const TString PARNM[3] = {"mass", "width", "BB"};
  const double OMEGAPAR_W0[3] = {782.547, 8.91001, 6.34484};
  const double OMEGAPAR_W0_ERR[3] = {0.0210935, 0.0476401, 4.37953e-2};
  
  const double OMEGAPAR_W[3] = {782.395, 9.17892, 6.09512};
  const double OMEGAPAR_W_ERR[3] = {0.0215785, 0.0487003, 4.14672e-2};

  double signifi_tmp = 0.;
  double diff_tmp = 0.;
  double diff_err_tmp = 0.;
  
  for (int i = 0; i < 3; i ++) {

    diff_tmp = TMath::Abs(OMEGAPAR_W0[i] - OMEGAPAR_W[i]);
    diff_err_tmp = TMath::Sqrt(OMEGAPAR_W0_ERR[i] * OMEGAPAR_W0_ERR[i] + OMEGAPAR_W_ERR[i] * OMEGAPAR_W_ERR[i]);
    signifi_tmp = diff_tmp / diff_err_tmp;
    
    cout << PARNM[i] << ": OMEGAPAR_W0 = " << OMEGAPAR_W0[i] << ", OMEGAPAR_W = " << OMEGAPAR_W[i] << ", relative diff. = " << TMath::Abs(OMEGAPAR_W0[i] - OMEGAPAR_W[i]) / OMEGAPAR_W0[i] * 100. << " %" << ", significance = " << signifi_tmp << endl;
      
  }
  
  //break;
  
  // Check relative difference, pi0 width measurement: Precision measurement of the neutral pion lifetime [/home/bo/Desktop/Works/eps2021/gg]
  cout << "\n\tCheck Relative Difference\n";

  const double Gamma_pi0_GBH = 8.10, Gamma_pi0_GBH_err = 0.08; // [eV], GBH, next-to-leading order
  const double Gamma_pi0_ChiAnom = 7.750, Gamma_pi0_ChiAnom_err = 0.016; // [eV]
  double Gamma_pi0_diff = TMath::Abs(Gamma_pi0_ChiAnom - Gamma_pi0_GBH) / Gamma_pi0_ChiAnom * 100.; // relative difference
  
  //double Gamma_pi0_diff

  cout << "Gamma_pi0_GBH = " << Gamma_pi0_GBH << " +/- " << Gamma_pi0_GBH_err << " [eV]\n"
       << "Estimated Uncertainty [%] = " << (Gamma_pi0_GBH_err / Gamma_pi0_GBH) * 100. << "\n"
       << "Gamma_pi0_ChiAnom = " << Gamma_pi0_ChiAnom << " +/- " << Gamma_pi0_ChiAnom_err << " [eV]\n"
       << "Gamma_pi0_diff [%] = " << Gamma_pi0_diff << "\n"; 

  const double Gamma_pi0_ChiAnom_pred = GetGammapi0();
  const double Gamma_pi0_ChiAnom_pred_err = GetGammapi0_err();

  cout << "\n\tCalculate pi0 width chiral anomaly\n"
       << "Gamma_pi0_Chi anomaly [eV] = " << Gamma_pi0_ChiAnom_pred << " +/- " << Gamma_pi0_ChiAnom_pred_err << "\n";

  const double sigma0 = 1457, sigma0_staterr = 23, sigma0_systerr = 19;
  const double Momega = 782.71, Momega_staterr = 0.07, Momega_systerr = 0.04;
  const double Gammaomega = 8.68, Gammaomega_staterr = 0.23, Gammaomega_systerr = 0.10;

  double GammaEEBr3pi = GetGammaEEBr3pi(sigma0, Momega, Gammaomega);
  double GammaEEBr3pi_StatErr = GetGammaEEBr3pi_ErrPropa(sigma0, Momega, Gammaomega, sigma0_staterr, Momega_staterr, Gammaomega_staterr);
  double GammaEEBr3pi_SystErr = GetGammaEEBr3pi_ErrPropa(sigma0, Momega, Gammaomega, sigma0_systerr, Momega_systerr, Gammaomega_systerr);
  
  cout << "\n\tCheck CMD2 results\n"
       << "Fit results: \n"
       << "sigma0 = " << sigma0 << " +/- " << sigma0_staterr << " +/- " << sigma0_systerr << " [nb]\n"
       << "M_omega = " << Momega << " +/- " << Momega_staterr << " +/- " << Momega_systerr << " [MeV]\n"
       << "Gamma_omega = " << Gammaomega << " +/- " << Gammaomega_staterr << " +/- " << Gammaomega_systerr << " [MeV]\n"
       << "Gamma_EE * Br_3pi x 1e3 = " << GammaEEBr3pi * 1e3 << " +/- " << GammaEEBr3pi_StatErr * 1e3 << " +/- " << GammaEEBr3pi_SystErr * 1e3 << endl;

  
  cout << "\n\tStat Err propagation\n";
    
  cout << "\n\tReover systematic errors %\n";
  const int err_list_size = 7;
  double Err_list[err_list_size] = {0.5, 0.1, 0.5, 0.1, 0.2, 0.3, 1.0};
  double errsum2 = 0., errsum_sqrt = 0.;
  
  for (int i = 0; i < err_list_size; i ++) {

    errsum2 += Err_list[i] * Err_list[i];
    errsum_sqrt = TMath::Sqrt(errsum2);
    cout << i + 1 << ": " << Err_list[i] << ", errsum2 = " << errsum2 << "%, errsum_sqrt = " << errsum_sqrt << "%\n";

    
  }

  ///  arXiv:hep-ph/9502406, omega -> 3pi contact term
  cout << "\n\tRecover alpha_K value\n" << endl; 

  double alpha_K = Getalpha_K();
  double F3pi = GetF3pi(mass_omega);
  
  cout << "alpha_K = " << alpha_K << "\n"
       << "F3pi = " << F3pi << "\n"
       << "gamma_euler = " << gamma_euler << "\n"
       << "mK = " << mK << "\n"
       << "mrho = " << mrho << "\n"
       << "f_pi = " << f_pi << "\n";

  // check F3pi function
  TF1 *ff_F3pi=new TF1("ff_F3pi", F3pi_fcn, 500., 800., 0);

  TCanvas *cv_F3pi = new TCanvas("F3pi", "F3pi", 700, 700);
  ff_F3pi -> Draw();
  
  /*
    gamma_euler = 0.5772; // Euler constant
    const double mK = 494.; // MeV Kaon mass
    const double mrho = 768.; // MeV rho mass
    const double f_pi = 93.; // MeV pi0 decay constant
  */

  // syst. error for min. E_clus cut

  const double E1 = 208.; //208
  const double E1err_real = 0.; // error of E1 read from the plot
  
  double Eerr_para = 0.; // cluster energy error evaluated using the parameterization
    
  Eerr_para = GetEerr(E1);
  
  cout << "\nsyst. error for min. E_clus cut\n"
       << "E1 = " << E1 << " MeV, Eerr_para = " << Eerr_para << " MeV, rel. err = " << Eerr_para / E1 << endl;

  // PDG 2022

  // branching ratio phi -> eta gamma
  const double B1 = 1.301e-2;
  const double B1_err = 0.025e-2;
  double B1_relerr2 = TMath::Power(B1_err/B1, 2);
  //cout << B1_relerr2 << endl;
  
  // branching ratio eta -> pi+ pi- pi0
  const double B2 = 23.02e-2;
  const double B2_err = 0.25e-2;
  double B2_relerr2 = TMath::Power(B2_err/B2, 2);
  
  // branching ratio pi0 -> gamma gamma
  const double B3 = 98.823e-2;
  const double B3_err = 0.034e-2;
  double B3_relerr2 = TMath::Power(B3_err/B3, 2);

  double BBB_err = (B1*B2*B3)*TMath::Sqrt(B1_relerr2 + B2_relerr2 + B3_relerr2);
    
  // branching ratio eta -> pi+ pi- gamma
  const double B5 = 4.28e-2;
  const double B5_err = 0.07e-2;
  const double B5_relerr2 = TMath::Power(B5_err/B5, 2);

  double BB_err = (B1*B5)*TMath::Sqrt(B1_relerr2 + B5_relerr2);
  
  cout << "\nProduct of the branching fractions\n"
       << "B(phi -> eta gamma)=" << B1 << " +/- " << B1_err << "\n"
       << "B(eta -> pi+ pi- pi0)=" << B2 << " +/- " << B2_err << "\n"
       << "B(pi0 -> gamma gamma)=" << B3 << " +/- " << B3_err << "\n"
       << "B(phi -> eta gamma) x B(eta -> pi+ pi- pi0) x B(pi0 -> gamma gamma) [10^-3] = " << B1 * B2 * B3 * 1e3 << " +/- " << BBB_err * 1e3 << "\n\n"
       << "B(phi -> eta gamma) x B(eta -> pi+ pi- gamma) [10^-3] = " << B1 * B5 * 1e3 << " +/- " << BB_err * 1e3 << endl;

  // Error propagation, c=axb
  
  cout << "\nErro propagation: 2D histo. bin width, Delta_Eisr X Delta_ppIM" << endl;
  const double Eisr_sigma = 1.24; // Eisr error
  const double ppIM_sigma = 1.15; // ppIM error
  const double Delta_ppIM = 2. * ppIM_sigma; // ppIM bin width
  const double Delta_Eisr = 2. * Eisr_sigma; // Eisr bin width

  double binWidth_2d = Delta_ppIM * Delta_Eisr;
  double binwidth_2d_rel_sigma = TMath::Sqrt(TMath::Power(1./2, 2)+TMath::Power(1./2, 2)); // relative error of 2d histo. bin width  

  // increase 5% of Delta_ppIM and Delta_Eisr at the same time, calcualte increase of binWidth_2d, to reach 1 sigma of binWidth_2d error, i.e. ~70%
  // test increase 50% of Delta_ppIM and Delta_Eisr (under the condition both bin width are twice the corresonding error) so that the 2D bin width increase by 125% 
  const double frac_Delta = 5 * 1e-2; // 0.5
  double Delta_ppIM_tmp = Delta_ppIM;
  double Delta_Eisr_tmp = Delta_Eisr;
  double binWidth_tmp = 0.;
  double binWidth_step_before = binWidth_2d;
  double binWidth_increase = 0., binWidth_increase_sum = 0.;
    
  for (int i=0; i < 7; i++) {

    binWidth_step_before = Delta_ppIM_tmp * Delta_Eisr_tmp;
    Delta_ppIM_tmp = Delta_ppIM_tmp * (1. + frac_Delta);
    Delta_Eisr_tmp = Delta_Eisr_tmp * (1. + frac_Delta);
    
    binWidth_tmp = Delta_ppIM_tmp * Delta_Eisr_tmp;
    binWidth_increase = (binWidth_tmp - binWidth_step_before) / binWidth_step_before * 100.;
    binWidth_increase_sum = binWidth_increase_sum + binWidth_increase;
    
    cout << i + 1 << ":\nDelta_ppIM = " << Delta_ppIM_tmp << ", Delta_Eisr = " << Delta_Eisr_tmp << "\n"
	 << "binWidth_tmp = " << binWidth_tmp << ", step before = " << binWidth_step_before << ", in creased " << binWidth_increase << " [%]" << ", sum increase = " << binWidth_increase_sum << " [%]\n\n";
  }
  
  cout << "Delta_ppIM = " << Delta_ppIM << " +/- " << ppIM_sigma << " [MeV/c^2]\n"
       << "Delta_Eisr = " << Delta_Eisr << " +/- " << Eisr_sigma << " [MeV]\n"
       << "2d histo. bin width = Delta_ppIM X Delta_Eisr = " << binWidth_2d << " [MeV/c^2]X[MeV]\n"
       << "2d histo. bin width absolute error = " << binwidth_2d_rel_sigma *  binWidth_2d << " [MeV/c^2]X[MeV]\n"
       << "2d histo. bin width relative error = " << binwidth_2d_rel_sigma *  100. << " [%]" << endl;
  
  cout << "\nIncrease 5% of Delta_ppIM and Delta_Eisr at the same time, calcualte increase of binWidth_2d, to reach 1 sigma of binWidth_2d error, i.e. ~70%\n"
       << "frac_Delta = " << frac_Delta << endl;
  

  cout << TMath::Power(1/2, 2) << endl;


  // test caluculate systematic error, recover Li's result, Li's thesis P.118 Tab. 4.8
  const int list_size = 9;
  const double SYSTERR_plus[list_size] = {9*1e-4, 1*1e-4, 9*1e-4, 0., 0., 0., 6*1e-4, 1*1e-3, 2*1e-4};
  const double SYSTERR_neg[list_size] = {9*1e-4, 1*1e-4, 9*1e-4, 1*1e-4, 6*1e-4, 0., 0., 1*1e-3, 2*1e-4};
  double systerr_plus=0., systerr_plus_sq_sum = 0.;
  double systerr_neg=0., systerr_neg_sq_sum = 0.;
    
  
  for (int i=0; i < list_size; i++) {

    systerr_plus_sq_sum += SYSTERR_plus[i] * SYSTERR_plus[i];
    systerr_neg_sq_sum += SYSTERR_neg[i] * SYSTERR_neg[i];
    
    cout << "systerr_plus_sq_sum = " << systerr_plus_sq_sum << "\n"
	 << "systerr_neg_sq_sum = " << systerr_neg_sq_sum << endl;
  }

  systerr_plus=TMath::Sqrt(systerr_plus_sq_sum);
  systerr_neg=TMath::Sqrt(systerr_neg_sq_sum);
    
  cout << "\n Syst. Err. = +" << systerr_plus << ", -" << systerr_neg << endl;

  cout << "\n Calculate uncorrelated error" << endl;
  double stat_errA = 0.0629984, stat_errB = 0.0630261, stat_errC = 0.0625479;
    
  double uncorr_err1 = TMath::Sqrt(TMath::Abs(stat_errA * stat_errA - stat_errB * stat_errB));
  double uncorr_err2 = TMath::Sqrt(TMath::Abs(stat_errC * stat_errC - stat_errB * stat_errB));

  cout << "uncorr_err1 = " << uncorr_err1 << ", uncorr_err2 = " << uncorr_err2 << endl;

  // check a following relation
  double m02_max = (Momega  - mass_pi0) * (Momega  - mass_pi0) * 1e-6;

  cout << "m02_max = " << m02_max << endl;
  
  return 0;
    
}

