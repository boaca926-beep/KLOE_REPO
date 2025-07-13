//
double GetCrx3piMax_new(double BB, double Mass_V) {

  double crx3pi = 12 * pi * BB / Mass_V / Mass_V * 1e6 * crx_scale;

  /*
  cout << "\n\tBB = " << BB << "\n"
       << "\tM_V = " << mass_omega << " [MeV]\n"
       << "\tcrx3pi = " << crx3pi << " [nb]\n";
  */
  
  return crx3pi;
  
}

double vmd_crx3pi_fcn(double *m, double *par) {// VMD 3pi cross section model function

  double M_V = par[0];
  double M_V2 = M_V * M_V;
  double M_phi_GeV = M_phi * .001;

  double Gamma_V = par[1];

  double BB = par[2];
  double Sigma_V = 12 * pi * BB / M_V / M_V * 1e6 * crx_scale;

  double Sigma_bg = par[3];
  
  // Get omega, phi amplitude
  double M_V_GeV = M_V * .001;
  double m_GeV = m[0] * .001;

  double F3pi_V = F3pi_num_fun(&M_V_GeV, &M_V_GeV); 
  double F3pi_m = F3pi_num_fun(&m_GeV, &m_GeV); // m variable
  double F3pi_phi = F3pi_num_fun(&M_phi_GeV, &M_phi_GeV);

  //cout << "F3pi_V = " << F3pi_V << ", F3pi_m = " << F3pi_m << ", F3pi_phi = " << F3pi_phi << endl; 
    
  double Gamma_Vs = Gamma_omegas_fcn(m[0], M_V, F3pi_V, Sigma_V, Gamma_V, F3pi_m);

  double N1 = M_V * M_V * Gamma_V * TMath::Sqrt(Sigma_V * M_V / F3pi_V);
  //double N1 = m[0] * m[0] * Gamma_V * TMath::Sqrt(Sigma_V * m[0] / F3pi_V); // equivalent to bw with F3pi = 1, and Gamma_Vs -> Gamma_V in term "b" 
  
  double a = m[0] * m[0] - M_V2;
  double b = m[0] * Gamma_Vs;

  double A = (N1 * a) / (a * a + b * b);
  double B = (N1 * b) / (a * a + b * b);
  
  //cout << Amp_V.Rho2() << ", " << N1 * N1 / (a * a + b * b) << ", " << A * A + B * B << endl;
  //cout << F3pi_V << endl;

  double Amp_V_par[5] = {M_V, F3pi_V, Gamma_V, Sigma_V, F3pi_m};
  TComplex Amp_V = Amp_omega_fcn(&m[0], Amp_V_par);

  /*
  double F3pi_phi = F3pi_num_fun(&M_phi_GeV, &M_phi_GeV);
  TComplex Amp_phi = Amp_phi_fcn(&m[0], &F3pi_phi);

  TComplex phase = TMath::Cos(phase_angle) + TMath::Sin(phase_angle) * TComplex::I();
  */
  
  double Gamma_phis = Gamma_phis_fcn(m[0], F3pi_phi);


  double N2 =  M_phi * M_phi * Gamma_phi * TMath::Sqrt(Sigma_phi * M_phi / F3pi_phi);
  //double N2 = m[0] * m[0] * Gamma_phi * TMath::Sqrt(Sigma_phi * m[0] / F3pi_phi);
  double c = m[0] * m[0] - M_phi * M_phi;
  double d = m[0] * Gamma_phis;

  double C = (N2 * c) / (c * c + d * d);
  double D = (N2 * d) / (c * c + d * d);

  //cout << "Gamma_Vs = " << Gamma_Vs << ", Gamma_phis = " << Gamma_phis << endl;
 
  //cout << Amp_phi.Rho2() << ", " << N2 * N2 / (c * c + d * d) << ", " << C * C + D * D << endl;

  
  //cout << phase_angle << endl;
  
  double E = TMath::Cos(phase_angle);
  double F = TMath::Sin(phase_angle);

  double G = C * E - D * F;
  double H = D * E + C * F ;
  
  //cout << "|phase|² = " << phase.Rho2() << endl;
 
  //TComplex Amp_bkg = Amp_bkg_fcn(Sigma_bg); //TMath::Power(M_V, 3 / 2) * TMath::Sqrt(Sigma_bg / F3pi_V) + TComplex::I() * 0.;

  double K = TMath::Power(M_V, 3 / 2) * TMath::Sqrt(Sigma_bg / F3pi_V);
  double L = 0.;
  
  //cout << Amp_bkg.Rho2() << endl;
  //cout << "Sigma_bg = " << endl;
  //cout << "|Amp_V|² = " << Amp_V.Rho2() << ", |Amp_phi|² = " << Amp_phi.Rho2() << ", |Amp_bkg|² = " << Amp_bkg.Rho2() << endl;

  double AmpSum_Rho2 = (A + G) * (A + G) + (B + H) * (B + H);
  //double AmpSum_Rho2 = (A) * (A) + (B) * (B);
  //double AmpSum_Rho2 = (G) * (G) + (H) * (H);

  
  //cout << "AmpSum_Rho2 = " << AmpSum_Rho2 << ", AmpSum_Rho2 (checked) = " << (Amp_V + phase * Amp_phi + Amp_bkg).Rho2() << endl;
  
  //double vmd_crx3pi = F3pi_m / TMath::Power(m[0], 3) * (Amp_V).Rho2(); // 3pi cross section
  //double vmd_crx3pi = F3pi_m / TMath::Power(m[0], 3) * (Amp_V + phase * Amp_phi).Rho2(); // 3pi cross section
  
  double vmd_crx3pi = F3pi_m / TMath::Power(m[0], 3) * AmpSum_Rho2; 

  // test bw 
  TComplex denomin(m[0] * m[0] - M_V * M_V, M_V * Gamma_V);
  double denomin_sqr = denomin.Rho2();
  
  
  double crx3pi_bw = Sigma_V * (m[0] * m[0]) * (Gamma_V * Gamma_V) / denomin_sqr;


  /*
  cout << "vmd_crx3pi = " << vmd_crx3pi << endl;
  cout << "m = " << m[0] << ", Gamma_V = " << Gamma_V << ", M_V = " << M_V << ", BB = " << BB << ", Sigma_V = " << Sigma_V << ", crx3pi_bw = " << crx3pi_bw << ", denomin.Rho2 = " << denomin.Rho2() << "\n\n";
  */
  
  
  //cout << vmd_crx3pi_phi << ", vmd_crx3pi = " << vmd_crx3pi << ", m = " << m[0] << ", Amp_phi.Im = " << Amp_phi.Im() << ", phase.Im = " << phase.Im() << endl;
  //cout << (Amp_phi * phase).Rho2() << ", " << Amp_phi.Rho2() << ", " << AmpSum_Rho2 << endl;
  
  //cout << "vmd_crx3pi: F3pi_s = " << F3pi_m << endl;
  //cout << "m_GeV = " << m_GeV  << ", F3pi_m = " << F3pi_m << " [GeV/c2], |Amp_omega|² = " << Amp_V.Rho2() << ", |Amp_phi|² = " << Amp_phi.Rho2() << ", |Amp_bkg|² = " << Amp_bkg.Rho2() << endl;

  /*
  cout << "M_V = " << M_V_GeV << ", s_sqrt = " << s_sqrt << "\n"
       << "F3pi_1 = " << F3pi_1 << ", F3pi_2 = " << F3pi_2 << "\n"
       << "F3pi_2 / F3pi_1 = " << F3pi_2 / F3pi_1 << "\n"
       << "F3pi_V = " << F3pi_V << ", F3pi_m = " << F3pi_m << endl;  
  */
  
  return vmd_crx3pi;
  
}


//
double Func_bw(double *x, double *par) {

  double bw = 0.;
  double m_3pi = x[0];
  double m_V = par[0]; // Vector meson mass in GeV/c^2
  double W_V = par[1]; // Vector meson mass width in GeV/c^2
  double BB = par[2]; // B_ee * B_3pi
  double sigma_max = GetCrx3piMax_new(BB, m_V);// 12. * pi * BB / (m_V * m_V); // Normalization factor, 12pi * B_ee * B_3pi / m^2_V

  // s-M^{2}_V + i * sqrt{s} * W_V
  double denomin_sqr_real = TMath::Power(m_3pi * m_3pi - m_V * m_V, 2.) + TMath::Power(m_V * W_V, 2.);
  TComplex denomin(m_3pi * m_3pi - m_V * m_V, m_V * W_V);
  TComplex denomin_sqr = denomin.Rho2();
  
  
  bw = sigma_max * (m_3pi * m_3pi) * (W_V * W_V) / denomin_sqr;
  //bw = 1. / denomin_sqr;

  /*
  cout << "\nm_3pi = " << m_3pi << " [GeV/c^2]\n"
    //<< "sqrtS = " << sqrtS << " [GeV/c^2]\n"
       << "para: m_V = " << m_V << " [GeV/c^2], W_V = " << W_V << " [GeV/c^2], BB = " << BB << "\ n"
       << "denomin = (" << denomin.Re() << ", " << denomin.Im() << "), Rho2 (check) = " << denomin.Rho2() << "(" << denomin_sqr_real << ") [GeV/c^2]\n"
       << "sigma_max = " << sigma_max << "\n" 
       << "bw = " << bw << "\n";
  */

  //cout << "sigma_max = " << sigma_max << ", m_V = " << m_V << ", BB = " << BB << endl;

  //cout << "m3pi = " << m_3pi << ", bw = " << bw << endl;

  return bw;//nb
  
}
//
double GetBB_new(double sigma_max, double Mass_V) {

  double BB = sigma_max * 1e-6 / crx_scale * Mass_V * Mass_V / 12. / pi;

  return BB;

}
