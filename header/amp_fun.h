double Gamma_phis_fcn(double m, double F3pi_V) {// Energy dependent width function for phi meson

  double M_V2 = M_phi * M_phi;
  double m2 = m * m;
  
  // phi -> K+ K-
  double FKpKm_s = F2pi_fcn(m, mK);
  double FKpKm_V = F2pi_fcn(M_phi, mK);
  double ratio_KpKm = Br_KpKm * M_V2 * FKpKm_s / m2 / FKpKm_V;

  //cout << FKpKm_s << ", " << FKpKm_V << ", " << ratio_KpKm << endl; 
  
  // phi -> KS KL
  double FKK_s = F2pi_fcn(m, mK0);
  double FKK_V = F2pi_fcn(M_phi, mK0);
  double ratio_KK = Br_KK * M_V2 * FKK_s / m2 / FKK_V;

  //cout << FKK_s << ", " << FKK_V << ", " << ratio_KK << endl; 

  // phi -> eta gamma
  double Fetagam_s = Fpigam_fcn(m, meta);
  double Fetagam_V = Fpigam_fcn(M_phi, meta);
  double ratio_etagam = Br_etagam * Fetagam_s / Fetagam_V;

  //cout << Fetagam_s << ", " << Fetagam_V << ", " << ratio_etagam << endl; 

  // phi -> pi+ pi- pi0
  double m_GeV = m * .001; //cout << m_GeV << endl;
  double M_phi_GeV = M_phi * .001; //cout << M_phi_GeV << endl;
  double F3pi_s = F3pi_num_fun(&m_GeV, &m_GeV);
  double F3pi_V_tmp = F3pi_num_fun(&M_phi_GeV, &M_phi_GeV);
  //double F3pi_V_tmp = F3pi_num_fun(&M_phi_GeV, &M_phi_GeV);
  
  double ratio_3pi = Br_3pi_phi * m * F3pi_s / M_phi / F3pi_V_tmp;
  
  // Get energy dependent width function
  double Gamma_Vs = Gamma_phi * (ratio_KpKm + ratio_KK + ratio_etagam + ratio_3pi + (1. - Br_3pi_phi - Br_etagam - Br_KK - Br_KpKm));
  //double Gamma_Vs = Gamma_phi * (ratio_KpKm + ratio_KK + ratio_etagam + ratio_3pi);

  //cout << Br_3pi_phi + Br_etagam + Br_KK + Br_KpKm << endl;
  
  //cout << "F3pi_s = " << F3pi_s << ", F3pi_V = " << F3pi_V << endl; 
  //cout << FKpKm_s << ", " << FKpKm_V << ", " << FKK_s << ", " << FKK_V << ", " << ratio_etagam << ", " << ratio_3pi << endl;
    
  return Gamma_Vs;
  
}


double Gamma_omegas_fcn(double m, double M_V, double F3pi_V, double Sigma_V, double Gamma_V, double F3pi_s) {// Energy dependent width function for omega meson

  double M_V2 = M_V * M_V;
  double m2 = m * m;

  //cout << m << ", " << M_V << endl;
  // omega -> pi+ pi-
  double F2pi_s = F2pi_fcn(m, mpi);
  double F2pi_V = F2pi_fcn(M_V, mpi);
  double ratio_2pi = Br_2pi * M_V2 * F2pi_s / m2 / F2pi_V;
  //double ratio_2pi = (1. - Br_pigam - Br_3pi_omega) * M_V2 * F2pi_s / m2 / F2pi_V;

  // omega -> pi0 gamma
  double Fpigam_s = Fpigam_fcn(m, mpi0);
  double Fpigam_V = Fpigam_fcn(M_V, mpi0);
  double ratio_pigam = Br_pigam * Fpigam_s / Fpigam_V;
  //double ratio_pigam = (1. - Br_2pi - Br_3pi_omega) * Fpigam_s / Fpigam_V;

  // omega -> pi+ pi- pi0
  double ratio_3pi = Br_3pi_omega * m * F3pi_s / M_V / F3pi_V;
  //double ratio_3pi = (1. - Br_2pi - Br_pigam) * m * F3pi_s / M_V / F3pi_V;

  //cout << M_V << endl;
  
  //cout << "Gamma_Vs: F3pi_s = " << F3pi_s << ", " << F3pi_V << endl;
  
  //double test_mass = 380.262 * .001;
  //cout << m << ", " << F3pi_s << ", " << M_V << ", " << F3pi_V << ", ratio_3pi = " << ratio_3pi << endl;
  
  double Gamma_Vs = Gamma_V * (ratio_2pi + ratio_pigam + ratio_3pi + (1. - Br_2pi - Br_pigam - Br_3pi_omega));
  //double Gamma_Vs = Gamma_V * (ratio_2pi + ratio_pigam + ratio_3pi);

  //cout << Br_2pi + Br_pigam + Br_3pi_omega << endl;
  //cout << m << ", " << ratio_2pi << ", " << ratio_pigam << ", " << ratio_3pi << ", " << Gamma_Vs << ", ratio_sum = " << ratio_2pi + ratio_pigam + ratio_3pi << endl; 
  
  return Gamma_Vs;
  
}

TComplex Amp_phi_fcn(double *m, double *par) {// Amplitude function; phi meson

  double F3pi_V = par[0];
  
  double M_V2 = M_phi * M_phi; // Vector mass squared

  double Gamma_Vs = Gamma_phis_fcn(m[0], F3pi_V); // Energy dependent width

  TComplex D_V = m[0] * m[0] - M_V2 + TComplex::I() * m[0] * Gamma_Vs; // Complex nominator
  
  TComplex Amp_V = M_V2 * Gamma_phi * TMath::Sqrt(Sigma_phi * M_phi / F3pi_V) / D_V; // Complex amplitude

  //cout << D_V.Rho2() << endl;
  //cout << "m = " << m[0] << ", M_V = " << M_phi << ", Gamma_V = " << Gamma_phi << ", Gamma_Vs = " << Gamma_Vs << ", Sigma_V = " << Sigma_phi << ", F3pi_V = " << F3pi_V << endl;
  
  return Amp_V;

}


TComplex Amp_omega_fcn(double *m, double *par) {// Amplitude function; omega meson

  double M_V = par[0];
  double F3pi_V = par[1];
  double Gamma_V = par[2];
  double Sigma_V = par[3];
  double F3pi_s = par[4];
  
  double M_V2 = M_V * M_V; // Vector mass squared

  double Gamma_Vs = Gamma_omegas_fcn(m[0], M_V, F3pi_V, Sigma_V, Gamma_V, F3pi_s); // Energy dependent width

  TComplex D_V = m[0] * m[0] - M_V2 + TComplex::I() * m[0] * Gamma_Vs; // Complex nominator
  
  TComplex Amp_V = M_V2 * Gamma_V * TMath::Sqrt(Sigma_V * M_V / F3pi_V) / D_V; // Complex amplitude
  //TComplex Amp_V = m[0] * m[0] * Gamma_Vs * TMath::Sqrt(Sigma_V * m[0] / F3pi_V) / D_V; // Complex amplitude

  //cout << Gamma_V << endl;
  //cout << "s_sqrt = " << m[0] << ", M_V = " << M_V << ", Gamma_V = " << Gamma_V << ", Gamma_Vs = " << Gamma_Vs << ", Sigma_V = " << Sigma_V << endl;

  //cout << Amp_V.Rho2() << endl;
  
  return Amp_V;

}

TComplex Amp_bkg_fcn(double Sigma_bkg_V) {// Amplitude function; background

  double M_omega_GeV = M_omega * .001; 
  double F3pi_omega = F3pi_num_fun(&M_omega_GeV, &M_omega_GeV);
  
  TComplex Amp_bkg = TMath::Power(M_omega, 3 / 2) * TMath::Sqrt(Sigma_bkg_V / F3pi_omega) + TComplex::I() * 0.;

  //cout << M_omega_GeV << ", " << Sigma_bkg_V << ", " << F3pi_omega << endl;
  
  return Amp_bkg;

}

/*
double g_rhopipi_fun(double m, double m1) {//rho -> 2pi transition coupling constant, assume to be equal for all processes

  double term = q0_fcn(m, m1);
  
  double g_rhopi = 6 * pi * mrho0 * mrho0 * Gamma_rho0 / TMath::Power(term, 3);

  return g_rhopi;
  
}
*/
