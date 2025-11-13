// BW formula in Ref: Physics Reports, 427 (2006) "Precision Electroweak Measurements on the Z Resonance", Eq. 1.41

double vmd_crx3pi_fcn(double *m, double *par) {// VMD 3pi cross section model function

  double M_V = par[0];
  double M_V2 = M_V * M_V;
  double M_phi_GeV = M_phi * .001;

  double Gamma_V = par[1];
  double Sigma_V = par[2];
  double Sigma_bkg_V = par[3];

  //cout << "s_sqrt = " << s_sqrt << endl;

  /*
  double test_mass = 0.5;
  double F3pi_test = F3pi_num_fun(&test_mass, &test_mass);

  cout << "\nM_V_GeV = " << M_V_GeV << " GeV/c², F3pi_V = " << F3pi_V << "\n" 
       << "M_phi_GeV = " << M_phi_GeV << " GeV/c², F3pi_phi = " << F3pi_phi << "\n"
       << "F3pi_test = " << F3pi_test << endl; 
  */
  
  // Get omega, phi and bkg amplitude
  double M_V_GeV = M_V * .001;
  double m_GeV = m[0] * .001;
  double F3pi_V = F3pi_num_fun(&M_V_GeV, &M_V_GeV); 
  double F3pi_m = F3pi_num_fun(&m_GeV, &m_GeV); // m variable
  
  //double Amp_V_par[5] = {M_V, F3pi_V, Gamma_V, Sigma_V, F3pi_m};
  //TComplex Amp_V = Amp_omega_fcn(&m[0], Amp_V_par);

  double Gamma_Vs = Gamma_omegas_fcn(m[0], M_V, F3pi_V, Sigma_V, Gamma_V, F3pi_m);
  //double N1 = m[0] * m[0] * Gamma_V * TMath::Sqrt(Sigma_V * m[0] / F3pi_V);
  double N1 = M_V * M_V * Gamma_V * TMath::Sqrt(Sigma_V * M_V / F3pi_V);
  
  double a = m[0] * m[0] - M_V2;
  double b = m[0] * Gamma_Vs;

  double A = (N1 * a) / (a * a + b * b);
  double B = (N1 * b) / (a * a + b * b);

  //cout << Amp_V.Rho2() << ", " << N1 * N1 / (a * a + b * b) << ", " << A * A + B * B << endl;
  //cout << F3pi_V << endl;

  double F3pi_phi = F3pi_num_fun(&M_phi_GeV, &M_phi_GeV);
  
  //TComplex Amp_phi = Amp_phi_fcn(&m[0], &F3pi_phi);

  double Gamma_phis = Gamma_phis_fcn(m[0], F3pi_phi);
  double N2 =  M_phi * M_phi * Gamma_phi * TMath::Sqrt(Sigma_phi * M_phi / F3pi_phi);
  //double N2 = m[0] * m[0] * Gamma_phi * TMath::Sqrt(Sigma_phi * m[0] / F3pi_phi);
  double c = m[0] * m[0] - M_phi * M_phi;
  double d = m[0] * Gamma_phis;

  double C = (N2 * c) / (c * c + d * d);
  double D = (N2 * d) / (c * c + d * d);

  //cout << Sigma_phi << endl;
  
  //cout << Amp_phi.Rho2() << ", " << N2 * N2 / (c * c + d * d) << ", " << C * C + D * D << endl;

  //TComplex phase = TMath::Cos(phase_angle) + TMath::Sin(phase_angle) * TComplex::I();

  //cout << phase_angle << endl;
  
  double E = TMath::Cos(phase_angle);
  double F = TMath::Sin(phase_angle);

  double G = C * E - D * F;
  double H = D * E + C * F ;
  
  //cout << "|phase|² = " << phase.Rho2() << endl;
 
  //TComplex Amp_bkg = Amp_bkg_fcn(Sigma_bkg_V); //TMath::Power(M_V, 3 / 2) * TMath::Sqrt(Sigma_bkg_V / F3pi_V) + TComplex::I() * 0.;

  double K = TMath::Power(M_V, 3 / 2) * TMath::Sqrt(Sigma_bkg_V / F3pi_V);
  double L = 0.;
  
  //cout << Amp_bkg.Rho2() << endl;
  //cout << "Sigma_bg = " << endl;
  //cout << "|Amp_V|² = " << Amp_V.Rho2() << ", |Amp_phi|² = " << Amp_phi.Rho2() << ", |Amp_bkg|² = " << Amp_bkg.Rho2() << endl;

  //double AmpSum_Rho2 = (A + G + K) * (A + G + K) + (B + H + L) * (B + H + L);
  //double AmpSum_Rho2 = (A + G) * (A + G) + (B + H) * (B + H);
  //double AmpSum_Rho2 = (A) * (A) + (B) * (B);
  double AmpSum_Rho2 = (G) * (G) + (H) * (H);

  
  //cout << "AmpSum_Rho2 = " << AmpSum_Rho2 << ", AmpSum_Rho2 (checked) = " << (Amp_V + phase * Amp_phi + Amp_bkg).Rho2() << endl;

  //double vmd_crx3pi = F3pi_m / TMath::Power(m[0], 3) * (Amp_V + phase * Amp_phi + Amp_bkg).Rho2();  double vmd_crx3pi = F3pi_m / TMath::Power(m[0], 3) * AmpSum_Rho2; 
  //double vmd_crx3pi = F3pi_m / TMath::Power(m[0], 3) * (Amp_phi).Rho2();
  double vmd_crx3pi = F3pi_m / TMath::Power(m[0], 3) * AmpSum_Rho2;
  
  //cout << vmd_crx3pi_phi << ", vmd_crx3pi = " << vmd_crx3pi << ", m = " << m[0] << ", Amp_phi.Im = " << Amp_phi.Im() << ", phase.Im = " << phase.Im() << endl;
  //cout << (Amp_phi * phase).Rho2() << ", " << Amp_phi.Rho2() << ", " << AmpSum_Rho2 << endl;
  
  //cout << "vmd_crx3pi: F3pi_s = " << F3pi_m << endl;
  //cout << vmd_crx3pi << endl;
  //cout << "m_GeV = " << m_GeV  << ", F3pi_m = " << F3pi_m << " [GeV/c2], |Amp_omega|² = " << Amp_V.Rho2() << ", |Amp_phi|² = " << Amp_phi.Rho2() << ", |Amp_bkg|² = " << Amp_bkg.Rho2() << endl;

  /*
  cout << "M_V = " << M_V_GeV << ", s_sqrt = " << s_sqrt << "\n"
       << "F3pi_1 = " << F3pi_1 << ", F3pi_2 = " << F3pi_2 << "\n"
       << "F3pi_2 / F3pi_1 = " << F3pi_2 / F3pi_1 << "\n"
       << "F3pi_V = " << F3pi_V << ", F3pi_m = " << F3pi_m << endl;  
  */
  
  return vmd_crx3pi;
  
}

double bw_crx3pi_fcn(double *m, double *par) {// Breit-Wigner 3pi cross section model function

  double M_V = par[0]; // Vector meson mass
  double M_V2 = M_V * M_V; // Vector meson mass squared
    
  double Gamma_V = par[1]; // Width
  double Gamma_V2 = Gamma_V * Gamma_V;

  double Sigma_V = par[2]; // Peak cross section

  double m2 = m[0] * m[0]; // Mass squared

  //TComplex Amp_bw = M_V * Gamma_V / (M_V2 - m2 - TComplex::I() * m[0] * Gamma_V); // Amplitude

  double Amp_bw_Rho2 = m2 * Gamma_V2 / (TMath::Power(M_V2 - m2, 2) + TMath::Power(m[0] * Gamma_V, 2));
  //double Amp_bw_Rho2 = M_V2 * Gamma_V2 / (TMath::Power(M_V2 - m2, 2) + TMath::Power(M_V * Gamma_V, 2));
    
  //cout << "|Amp_bw|² = " << Amp_bw.Rho2() << ", Amp_bw_Rho2 = " << Amp_bw_Rho2 << endl;
  
  double bw_crx3pi = Sigma_V * Amp_bw_Rho2;
  //double bw_crx3pi = Sigma_V * Amp_bw.Rho2();

  //cout << bw_crx3pi << endl;
  
  return bw_crx3pi;
  
}

double bw_crx3pi_conv_fcn(double *m, double *par) {// Breit-Wigner 3pi cross section model function for convolution BW(s-m2), for integration over mass mean given mass m

  double M_V = par[0]; // Vector meson mass
  double M_V2 = M_V * M_V; // Vector meson mass squared
    
  double Gamma_V = par[1]; // Width
  double Gamma_V2 = Gamma_V * Gamma_V;

  double Sigma_V = par[2]; // Peak cross section

  double m2 = (m[0] - par[3]) * (m[0] - par[3]); // Mass squared

  double Amp_bw_Rho2 = M_V2 * Gamma_V2 / (TMath::Power(M_V2 - m2, 2) + TMath::Power(M_V * Gamma_V, 2));
    
  
  double bw_crx3pi = Sigma_V * Amp_bw_Rho2;
  
  return bw_crx3pi;
  
  
}

double get_frac_fcn(double m, double sqrt_s) {

  double m2 = m * m;
  double s = sqrt_s * sqrt_s;

  double frac = 1 - m2 / s;
  
  return frac;
}

double W0_fcn(double *m, double *par) {// leading order radiator function, nozero theta0

  double x = get_frac_fcn(m[0], par[0]);
  
  double norm_factor = alpha / pi / x;
  double term1 = 2 - 2 * x + x * x;
  double term2 = TMath::Log((1 + TMath::Cos(theta0)) / (1 - TMath::Cos(theta0)));
  double term3 = x * x * TMath::Cos(theta0);
  
  double W0 = norm_factor * (term1 * term2 - term3);
  
  //cout << "sqrt_s = " << par[0] << ", m = " << m[0] << ", x = " << x << ", norm_factor = " << norm_factor << ", W0 = " << W0 << endl;
  
  return W0;

}

double W0_full_fcn(double *m, double *par) {

  double x = get_frac_fcn(m[0], par[0]);

  double alphapi = alpha / pi;
  double sqrt_s = par[0];
  
  double W0_full = (alphapi / x) * (TMath::Log(TMath::Power(sqrt_s / me, 2)) - 1) * (2 - 2 * x + x * x);

  //cout << "W0_full = " << W0_full << endl;
    
  return W0_full;
}

double W0_full_conv_fcn(double *m, double *par) {

  
  double x = get_frac_fcn(m[0] - par[1], par[0]);

  double alphapi = alpha / pi;
  double sqrt_s = par[0];
  
  double W0_full = (alphapi / x) * (TMath::Log(TMath::Power(sqrt_s / me, 2)) - 1) * (2 - 2 * x + x * x);

  //cout << "W0_full = " << W0_full << endl;
    
  return W0_full;
}


double W1_fcn(double *m, double *par) {

  double delta = 0., beta = 0., L = 0.;
  double x = 0., A = 0., a = 0., mm = 0, B = 0., C = 0., d = 0., D = 0.;
  double y1 = TMath::Cos(theta0); //cout << theta0 << endl;
  double W1 = 0.;
  double W0_temp = 0.;
  double alphapi = alpha / pi;
  double sqrtS = par[0];
  
  mm = m[0];
  x = get_frac_fcn(m[0], par[0]);
  //x = 1. - TMath::Power(mm / sqrtS,2.);
  L = TMath::Log(TMath::Power(sqrtS / me, 2.)); 
  beta = 2. * alphapi * (L - 1.);
  a  = beta - 1.;	
  A = beta * TMath::Power(x, a);
  B = 1. + alphapi * (pi * pi / 3. - 1. / 2.) + (3. / 4.) * beta - (beta * beta /24.) * (L / 3. + 2. * pi * pi - 37. / 4.);
  C = beta * (1. - x / 2.);
  d = (1. + 3. * (1. - x) * (1. - x));
  D = (beta * beta /8.) * (4. * (2. - x) * TMath::Log(1. / x) - d * TMath::Log(1. - x) / x - 6. + x);

  W1 = A * B - C + D; 

  //cout << mm << ", " << x << ", " << pi << ", " << beta << ", " << B << endl;
  
  return W1;

}

/*
double bw_crx3pi_ISR_fcn(double *s_sqrt, double *par) {

  double m_max = s_sqrt[0];
  double m_min = par[0];
  
  TF1 *integrand_fcn = new TF1("integrand_fcn", m_min, m_max, 1);
  integrand_fcn -> SetParameter(0, s_sqrt[0]);

  
  double crx3pi = integrand_fcn -> Integral(m_min, m_max);

  cout << "crx3pi ISR distorted = " << crx3pi << ", m_min = " << m_min << ", m_max = " << m_max << endl;
  
  return crx3pi;
  
  
}
*/


double integrand_fcn(double *m, double *par) {

  double sqrt_s = par[0];
  double s = sqrt_s * sqrt_s;
  
  // calculate W0
  double W0 = W0_full_fcn(&m[0], &sqrt_s);

  // calcualte crx3pi
  double bw_crx3pi_par[3] = {M_omega, Gamma_omega, Sigma_omega};
  double bw_crx3pi = bw_crx3pi_fcn(&m[0], bw_crx3pi_par);

  double norm_factor = 2 * m[0] / s; 
  double integrand = norm_factor * W0 * bw_crx3pi; // function of m & s

  //cout << "integrand = " << integrand << endl; 
  //cout << "sqrt_s = " << par[0] << ", m = " << m[0] << ", norm_factor = " << norm_factor << ", W0 = " << W0 << ", bw_crx3pi = " << bw_crx3pi << ", integrand = " << integrand << endl;
  
  return integrand;
  
}
