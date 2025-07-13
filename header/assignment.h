const double mass_mu = 105.6583755; //105.6583755±0.0000023 MeV
const double mass_pi = 139.57039; // MeV, mass of pi+/-
const double mass_pi0 = 134.9768; // MeV, mass of pi0
const double mass_pi0_err = 5e-4; // MeV, mass of pi0 error
const double mK = 494.; // MeV Kaon mass
const double mrho = 768.; // MeV rho mass
const double me = 0.511 * 1e-3; // mass of electron in GeV
const double mass_omega = 782.65; // omega mass in MeV/c^2

const double N_c = 3.; // number of colors in QCD
const double F_pi = 92.277; // MeV, pion decay constant
const double F_pi_err = 0.095; // MeV, error or pion decay constant
const double f_pi = 93.; // MeV pi0 decay constant
const double Gam_omega = 8.49; // omega mass width in MeV/c^2

const double BB_omega = 6.38e-5; // branching ratio product
const double gamma_euler = 0.5772; // Euler constant
  
const double pi = TMath::Pi();
const double alpha = 1. / 137; // fine structure constant 
const double alphapi = alpha / pi;
const double sqrtS = 1.019; // in GeV 
const double theta0_degr = 23.;
const double theta0 = theta0_degr * pi / 180.; // in radius
//const double sample_size = 2000;
//const double Lumi_int = 124529 * (sample_size / 642 ); //nb^-1
//const double Lumi_int = 1724470 * (sample_size / 8373); //nb^-1
const double crx_scale = 0.3894e6; // 1 [GeV]^{-2}=0.3894 mb=0.3894e6nb

const double eeg_lsf = 2., allphys_lsf = 1.;

/// Functions
double kernelFcn(double *s, double *p) {

  double mass_mu_GeV = mass_mu * 1e-3;
  
  double beta_mu = TMath::Sqrt(1. - 4. * mass_mu_GeV * mass_mu_GeV / s[0]);

  //cout << 1. - 4. * mass_mu_GeV * mass_mu_GeV / s[0] << endl;
  
  double x = (1. - beta_mu) / (1. + beta_mu);

  double term1 = x * x / 2. * (2. - x * x);

  double term2 = (1 + x * x) * TMath::Power(1 + x, 2) / x / x;

  double term3 = TMath::Log(1 + x) - x + x * x / 2.;

  double term4 = (1 + x) / (1 - x) * TMath::Power(x, 2) * TMath::Log(x);
  
  double Ks = (3. * s[0]) / mass_mu_GeV / mass_mu_GeV * (term1 + term2 * term3 + term4);
  //double Ks = term1 + term2 * term3 + term4;
  
  //cout << "muon mass: " << p[0] << endl;
  
  return Ks;
}

double GetEerr(double E) {

  double err = E * 5.7 * 1e-2 / TMath::Sqrt(E * 1e-3);
  
  return err;
}

double F3pi_fcn(double *x, double *p) {

  double e_value = TMath::Sqrt(4 * pi * alpha);

  double ratio = TMath::Log((pi * mK * mK) / (x[0] * x[0]));

  double value = -3 * e_value * TMath::Power(mK, 2) / 4. / TMath::Power(2 * pi, 4) / TMath::Power(f_pi, 5) * (ratio + gamma_euler);

  return value;
}

double GetF3pi(double x){

  double e_value = TMath::Sqrt(4 * pi * alpha);

  double ratio = TMath::Log((pi * mK * mK) / (x * x));

  cout << "ratio = " << ratio << endl;

  double value = -3 * e_value * TMath::Power(mK, 2) / 4. / TMath::Power(2 * pi, 4) / TMath::Power(f_pi, 5) * (ratio + gamma_euler);

  return value;
}

double Getalpha_K() {

  double ratio = TMath::Log((pi * mK * mK) / (mrho * mrho));
  
  double value = (mK / (2. *TMath::Sqrt(2) * pi * f_pi)) * TMath::Sqrt(ratio + gamma_euler);

  return value;

}

double GetGammaEEBr3pi_ErrPropa(double c, double d, double e, double c_err, double d_err, double e_err) {
  //c: sigma0, d: Momega, e: Gammaomega
  //x = c * d * d * e
  //a: normalization factor
  const double a = 1e-6 / crx_scale / 12 / pi;
  double x = c * d * d * e;
  //double delta_x = c_err * d_err * d_err * e_err;
  //double c_sumerr2 = TMath::Power(c_err / c, 2);
  //double d_sumerr2 = TMath::Power(d_err / d, 2);
  //double e_sumerr2 = TMath::Power(e_err / e, 2);
  //double sumerr2 = c_sumerr2 + 2 * d_sumerr2 + e_sumerr2;

  double err = a * x * TMath::Sqrt(TMath::Power(c_err / c, 2) + 2 * TMath::Power(2 * d_err / d, 2) + TMath::Power(e_err / e, 2)); 

  //cout << crx_scale << endl;
  /*
  cout << "a = " << a << "\n"
       << "c = " << c << ", d = " << d << ", e = " << e << "\n"
       << "c_err = " << c_err << ", d_err = " << d_err << ", e_err = " << e_err << "\n"
       << "sumerr2 = " << sumerr2 << "\n"
       << "err x 10^3 = " << err * 1e3 << "\n"
       << "delta_x = " << delta_x << "\n";
  */

  return err;

  
}

double GetGammaEEBr3pi(double sigma_V, double M_V, double Gamma_V) {

  double BB =  1e-6 / crx_scale / 12. / pi * sigma_V * M_V * M_V;
  double value = BB * Gamma_V;

  cout << "value x 1e3 = " << value * 1e3 << endl;
  
  return value;
}

double GetGammapi0_err() {// convert to eV

  double Gamma_ratio = mass_pi0 / F_pi;
  double Norm = TMath::Power(alpha * N_c, 2.) / TMath::Power(pi, 3.) / 576;
  double err2 = (9 * TMath::Power(Gamma_ratio, 4) * TMath::Power(mass_pi0_err, 2) + 4 * TMath::Power(Gamma_ratio, 6) * TMath::Power(F_pi_err, 2));

  double err = Norm * TMath::Sqrt(err2) * 1e6;
  
  return err;
}

double GetGammapi0() {

  double width = mass_pi0 * 1e6 * TMath::Power(mass_pi0 / F_pi, 2.) * TMath::Power(alphapi * N_c, 2.) / (576 * pi);

  return width;
}

double GetPValue(double sigma_nb) {
  double p = ROOT::Math::erfc(sigma_nb / TMath::Sqrt(2)) / 2.;

  //cout << "p = " << p << endl;
  return p;
  
}
