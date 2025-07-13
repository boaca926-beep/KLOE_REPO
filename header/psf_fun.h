double mpi0_GeV = 134.9768 * .001; // neutral pi mass
double mpi_GeV = 139.57039 * .001; // charged pi mass
//mrho0 = mrho0 * .001; // rho mass
//Gamma_rho0 = Gamma_rho0 * .001; // neutral rho width 

double E2_cm_fun(double m12, double m1, double m2) {// E2 is the energy of particle 2 in the m12 c.m. system

  double E2_cm = (m12 * m12 - m1 * m1 + m2 * m2) / (2 * m12);
  
  return E2_cm;
  
}

double E3_cm_fun(double m, double m12, double m3) {// E3 is the energy of particle 3 in the m12 c.m. system

  double E3_cm  = (m * m - m12 * m12 - m3 * m3) / (2 * m12);
  
  return E3_cm;
  
}

double Fpigam_fcn(double m, double M_m) {// M_m is pi0 or eta
  
  double m2 = m * m; // mass squared
  
  double M_m2 = M_m * M_m; // meson mass squared

  double Fpigam = 0.;

  if ((1. - M_m2 / m2) > 0.) {
    Fpigam = TMath::Power(m * (1. - M_m2 / m2), 3);
  }
  else {
    Fpigam = 0.;
  }
  
  return Fpigam;
  
}

double F2pi_fcn(double m, double M_m) {// M_m is charged meson messes: K+/K- or pi+/pi-

  double m2 = m * m; // 3pi mass squared

  double M_m2 = M_m * M_m; // charged meson mass squared

  double F2pi = 0.;

  if ((m2 / 4. - M_m2) > 0.) {
    F2pi = TMath::Power((m2 / 4. - M_m2), 3 / 2);
  }
  else {
    F2pi = 0.;
  }
  
  return F2pi;
  
}

double dalitz_down_fun(double *m, double *par) {

  double m12 = m[0];
  double m12sq = m12 * m12;
  double M = par[0];
  
  double m12sq_min = 4 * mpi_GeV * mpi_GeV; // GeV
  double m12sq_max = (M - mpi0_GeV) * (M - mpi0_GeV); // GeV
 
  double E2_cm = 0., E3_cm = 0.;
  double term1 = 0., term2 = 0., term3 = 0., term4 = 0.;
  double dalitz_down = 0.;

  if (m12sq <= m12sq_max && m12sq >= m12sq_min) {
    E2_cm = E2_cm_fun(m12, mpi_GeV, mpi_GeV);
    E3_cm = E3_cm_fun(M, m12, mpi0_GeV);

    term1 = E2_cm + E3_cm;
    term2 = TMath::Sqrt(E2_cm * E2_cm - mpi_GeV * mpi_GeV);
    term3 = TMath::Sqrt(E3_cm * E3_cm - mpi0_GeV * mpi0_GeV);
    term4 = term2 + term3;

    dalitz_down = term1 * term1 - term4 * term4;

  }
  else {

    dalitz_down = 0.;
    
  }

  /*
  cout << "\nE2_cm = " << E2_cm << ", E3_cm = " << E3_cm << "\n"
       << "m12sq = " << m12sq << ", (m12sq_min, m12sq_max) = (" << m12sq_min << ", " << m12sq_max << ")\n"  
       << "dalitz_down = " << dalitz_down << endl;
  */
  
  return dalitz_down;
  
}

double dalitz_upp_fun(double *m, double *par) {

  double m12 = m[0];
  double m12sq = m12 * m12;
  double M = par[0];
 
  double m12sq_min = 4 * mpi_GeV * mpi_GeV; // GeV
  double m12sq_max = (M - mpi0_GeV) * (M - mpi0_GeV); // GeV
  
  double E2_cm = 0., E3_cm = 0.;
  double term1 = 0., term2 = 0., term3 = 0., term4 = 0.;
  double dalitz_upp = 0.;

  if (m12sq <= m12sq_max && m12sq >= m12sq_min) {
    E2_cm = E2_cm_fun(m12, mpi_GeV, mpi_GeV);
    E3_cm = E3_cm_fun(M, m12, mpi0_GeV);

    term1 = E2_cm + E3_cm;
    term2 = TMath::Sqrt(E2_cm * E2_cm - mpi_GeV * mpi_GeV);
    term3 = TMath::Sqrt(E3_cm * E3_cm - mpi0_GeV * mpi0_GeV);
    term4 = term2 - term3;

    dalitz_upp = term1 * term1 - term4 * term4;

  }
  else {

    dalitz_upp = 0.;
    
  }

  /*
  cout << "\nE2_cm = " << E2_cm << ", E3_cm = " << E3_cm << "\n"
       << "m12sq = " << m12sq << ", (m12sq_min, m12sq_max) = (" << m12sq_min << ", " << m12sq_max << ")\n"  
       << "dalitz_upp = " << dalitz_upp << endl;
  */
  
  return dalitz_upp;
  
}

double getmomsq_fun(double m1, double m2, double m3) {// Momentum squared for pions
  
  double value_temp1 = m1 * m1 - (TMath::Sqrt(m2) + m3) * (TMath::Sqrt(m2) + m3);
  double value_temp2 = m1 * m1 - (TMath::Sqrt(m2) - m3) * (TMath::Sqrt(m2) - m3);
  
  double getmomsq = value_temp1 * value_temp2 / (4 * m1 * m1);	

  /*
  cout << "\nZ = " << getmomsq << ", value_temp1 = " << value_temp1 << "\n"
       << "m1 = " << m1 << ", m2 = " << m2 << ", m3 = " << m3 << endl;
  */
  
  return getmomsq;
}

double Get_k_fcn(double a, double b, double c, double d) {// Get k component for complex number h + ki
  double f = a * c - b * d;
  double l = a * d + b * c;
  
  double k = l / (f * f + l * l);
  
  return k;
  
}

double Get_h_fcn(double a, double b, double c, double d) {// Get h component for complex number h + ki
  double f = a * c - b * d;
  double l = a * d + b * c;
  
  double h = f / (f * f + l * l);
  
  return h;
  
}


double Amp_rho_sq_fcn(double m_0, double m_plus, double m_minus, double m) {// sum of all rho -> 2pi amplitude squared

  //double s_sqrt = m[0];
  
  // Get q0(m_0)
  double q0_m_0 = 0.;
  double q0_mrho0 = 0.;
  
  if ((m_0 * m_0 - 4 * mpi * mpi) > 0.) {
    q0_m_0 = 0.5 * TMath::Sqrt(m_0 * m_0 - 4 * mpi * mpi); 
  }
  else {
    q0_m_0 = 0.;
    //cout << "negative!!!" << endl;
  }

  if ((mrho0 * mrho0 - 4 * mpi * mpi) > 0.) {
    q0_mrho0 = 0.5 * TMath::Sqrt(mrho0 * mrho0 - 4 * mpi * mpi); 
  }
  else {
    q0_mrho0 = 0.;
    //cout << "negative!!!" << endl;
  }

  // rho -> pi pi transition coupling constants
  
  double g_rhopi_sq = 6 * pi * mrho0 * mrho0 * Gamma_rho0 / TMath::Power(q0_mrho0, 3);

  //cout << g_rhopi_sq << endl;
  
  // Get q_plus(m_plus)
  double m_minus_2 = m_minus * m_minus; 

  double q_m_plus_sq = getmomsq_fun(m_plus, mpi0 * mpi0, mpi); 
  double q_m_plus = 0.;

  double q_mrho_plus_sq = getmomsq_fun(mrho0, mpi0 * mpi0, mpi);
  double q_mrho_plus = 0.;

  if (q_mrho_plus_sq > 0.) {
    q_mrho_plus = TMath::Sqrt(q_mrho_plus_sq);
    //cout << q_mrho_plus << endl;
  }
  else {
    q_mrho_plus = 0.;
    //cout << "negative!!!" << endl;
  }
  
  if (q_m_plus_sq > 0.) {
    q_m_plus = TMath::Sqrt(q_m_plus_sq);
    //cout << q_m_plus << endl;
  }
  else {
    q_m_plus = 0.;
    //cout << "negative!!!" << endl;
  }

  double q_m_minus_sq = getmomsq_fun(m_minus, mpi0 * mpi0, mpi); 
  double q_m_minus = 0.;

  double q_mrho_minus_sq = q_mrho_plus_sq;
  double q_mrho_minus = TMath::Sqrt(q_mrho_plus_sq);

  
  if (q_m_minus_sq > 0.) {
    q_m_minus = TMath::Sqrt(q_m_minus_sq);
    //cout << q_m_minus << endl;
  }
  else {
    q_m_minus = 0.;
    //cout << "negative!!!" << endl;
  }

  
  // Get width
  double mass_ratio_rho0 = mrho0 / m_0;  
  double moment_ratio_rho0 = q0_m_0 / q0_mrho0;

  double mass_ratio_rho_plus = mrho0 / m_plus;  
  double moment_ratio_rho_plus = q_m_plus / q_mrho_plus;

  double mass_ratio_rho_minus = mrho0 / m_minus;  
  double moment_ratio_rho_minus = q_m_minus / q_mrho_minus;

  //(cout << q_m_minus << ", " << q_mrho_minus << endl;
  
  
  double Gamma_Vs_rho0 = TMath::Power(mass_ratio_rho0, 2) * Gamma_rho0 * TMath::Power(moment_ratio_rho0, 3); // Energy dependent width
  double Gamma_Vs_rho_plus = TMath::Power(mass_ratio_rho_plus, 2) * Gamma_rho0 * TMath::Power(moment_ratio_rho_plus, 3);
  double Gamma_Vs_rho_minus = TMath::Power(mass_ratio_rho_minus, 2) * Gamma_rho0 * TMath::Power(moment_ratio_rho_minus, 3);

  /*
  cout << "\nmrho0 = " << mrho0 << ", mpi = " << mpi << ", Gamma_rho0 = " << Gamma_rho0 << "\n"
       << "m_0 = " << m_0 << "," << ", m_plus = " << m_plus << ", m_minus = " << m_minus << "\n"
       << "mass_ratio_rho0 = " << mass_ratio_rho0 << ", moment_ratio_rho0 = " << moment_ratio_rho0 << ", Gamma_Vs_rho0 = " << Gamma_Vs_rho0 << "\n"
       << "mass_ratio_rho_plus = " << mass_ratio_rho_plus << ", moment_ratio_rho_plus = " << moment_ratio_rho_plus << ", Gamma_Vs_rho_plus = " << Gamma_Vs_rho_plus << "\n"
       << "mass_ratio_rho_minus = " << mass_ratio_rho_minus << ", moment_ratio_rho_minus = " << moment_ratio_rho_minus << ", Gamma_Vs_rho_minus = " << Gamma_Vs_rho_minus << "\n";
  */
  
  //TComplex D_V = mrho0 * mrho0 - m_0 * m_0 + TComplex::I() * m_0 * Gamma_Vs_rho0; // Complex nominator

  //TComplex Amp_rho0 = 1. / D_V;

  //double Amp_rho02_checked = 1. / (TMath::Power(mrho0 * mrho0 - m_0 * m_0, 2) + TMath::Power(m_0 * Gamma_Vs, 2));// denominator squared 

  //cout << "Amp_rho0.Rho2() = " << Amp_rho0.Rho2() << ", Amp_rho02_checked = " << Amp_rho02_checked << endl;

  // test complex products; 1: rho0, 2: rho+, 3: rho-
  // complex number in form of D = a - bi, Z = c - di
   
  double a1 = mrho0 * mrho0 - m_0 * m_0;
  double b1 = m_0 * Gamma_Vs_rho0;
  double c1 = 1.;
  double d1 = 0.; //s1 * PHI_apprx_fcn(&m_0, &s_sqrt);

  //cout << "d1 = " << d1 << ", m_0 = " << m_0 << ", s_sqrt = " << s_sqrt << endl;
  
  double a2 = mrho0 * mrho0 - m_plus * m_plus;
  double b2 = m_plus * Gamma_Vs_rho_plus;
  double c2 = c1;
  double d2 = 0.; //s1 * PHI_apprx_fcn(&m_plus, &s_sqrt);

  //cout << "d1 = " << d1 << ", m_plus = " << m_plus << ", s_sqrt = " << s_sqrt << endl;
  
  double a3 = mrho0 * mrho0 - m_minus * m_minus;
  double b3 = m_minus * Gamma_Vs_rho_minus;
  double c3 = c1;
  double d3 = 0.; //s1 * PHI_apprx_fcn(&m_minus, &s_sqrt);
  
  //cout << "d1 = " << d1 << ", m_minus = " << m_plus << ", s_sqrt = " << s_sqrt << endl;
  
  //TComplex D1 = a1 - TComplex::I() * b1;
  //TComplex Z1 = c1 - TComplex::I() * d1;
  double h1 = Get_h_fcn(a1, b1, c1, d1);
  double k1 = Get_k_fcn(a1, b1, c1, d1);
    
  //TComplex D2 = a2 - TComplex::I() * b2;
  //TComplex Z2 = c2 - TComplex::I() * d2;
  double h2 = Get_h_fcn(a2, b2, c2, d2);
  double k2 = Get_k_fcn(a2, b2, c2, d2);
  
  //TComplex D3 = a3 - TComplex::I() * b3;
  //TComplex Z3 = c3 - TComplex::I() * d3;
  double h3 = Get_h_fcn(a3, b3, c3, d3);
  double k3 = Get_k_fcn(a3, b3, c3, d3);
  
  //TComplex Amp_rho1 = g_rpp / D1 / Z1;
  //TComplex Amp_rho2 = g_rpp / D2 / Z2;
  //TComplex Amp_rho3 = g_rpp / D3 / Z3;
  
  //TComplex Amp_rho_sum_checked = Amp_rho1 + Amp_rho2 + Amp_rho3;

  //cout << "a3pi = " << a3pi << endl;
  
  double h_sum = h1 + h2 + h3 + a3pi / TMath::Sqrt(gpipi2);
  double k_sum = k1 + k2 + k3;
  
  //double Amp_rho_sq = g_rhopi_sq * (h_sum * h_sum + k_sum * k_sum);
  double Amp_rho_sq = gpipi2 * (h_sum * h_sum + k_sum * k_sum);

  //cout << g_rhopi_sq << endl;
  
  /*
  cout << "h1 = " << h1 << ", k1 = " << k1 << "\n" 
       << "h2 = " << h2 << ", k2 = " << k2 << "\n" 
       << "h3 = " << h3 << ", k3 = " << k3 << "\n"; 
  */
  
  //cout << "Amp_rho_sq (real) = " << Amp_rho_sq << ", Amp_rho_sq (complex) = " << Amp_rho0.Rho2() << endl;

  /*
  cout << "Amp_rho_sum_Rho2 (complex) = " << Amp_rho_sum_checked.Rho2() << "\n"
       << "Amp_rho_sum_Rho2 (real) = " << Amp_rho_sq << "\n"
       << "h1 = " << h1 << ", k1 = " << k1 << "\n" 
       << "h2 = " << h2 << ", k2 = " << k2 << "\n" 
       << "h3 = " << h3 << ", k3 = " << k3 << "\n"; 
  */
  
  return Amp_rho_sq;
  
}

// m: mass. m_0 = |(P_pi_minus + P_pi_plus)²|, m_minus = |(P_pi_minus + P_pi_0)²|, m_plus = |(P_pi_plus + P_pi_0)²|. All masses are in GeV!
double integrand_psf_fun(double *x, double *par){// Integrand for psf for numerical integration;

  double m_0_2 = x[0];
  double m_0 = TMath::Sqrt(x[0]);

  double m_plus_2 = par[0];
  double m_plus = TMath::Sqrt(par[0]);

  double m = par[1];

  double m_minus_2 = m * m + 2. * mpi_GeV * mpi_GeV + mpi0_GeV * mpi0_GeV - m_0_2 - m_plus_2; 
  double m_minus = TMath::Sqrt(m_minus_2);

  //Test PHI funciton
  //double PHI = PHI_apprx_fcn(&m_0, &m);

  //cout << "m_0 = " << m_0 << ", s_sqrt = " << m << ", PHI_apprx = " << PHI_apprx << endl;
  //cout << m << endl;
  
  double mrho0_GeV = mrho0;
  double m_0_min = 2 * mpi_GeV;
  double m_0_max = m - mpi0_GeV; 
  
  double moment_pion_0_2 = getmomsq_fun(m, m_0_2, mpi0_GeV);
  double moment_pion_plus_2 = getmomsq_fun(m, m_minus_2, mpi_GeV);
  double moment_pion_minus_2 = getmomsq_fun(m, m_plus_2, mpi_GeV);
  
  double X = moment_pion_0_2;
  double Y = moment_pion_minus_2;
  double Z = moment_pion_plus_2;
  
  //Test X, Y and Z
  //M = 0.782, (xx_temp_sqrt, yy_temp_sqrt) = (0.513776, 0.417773), X = 0.0364775, Y = 0.0657026, Z = 0.0461931
  
  double moment_pion_cross_sq = Y * Z / 2. + X * (Y + Z) / 2. - (X * X + Y * Y + Z * Z) / 4.; // Squared cross product of pi_plus and pi_minus

  /*
  cout << "\nm = " << m << "\n"
       << "(m_0, m_plus) = (" << m_0 << ", " << m_plus << ")\n"
       << "X = " << X << ", Y = " << Y << ", Z = " << Z << "\n"
       << "moment_pion_cross_sq = " << moment_pion_cross_sq << endl;
  */

  double Amp_rho_sq = Amp_rho_sq_fcn(m_0, m_plus, m_minus, m);

  //cout << Amp_rho_sq << endl;
    
  double integrand_psf = m_0 * m_plus * moment_pion_cross_sq * Amp_rho_sq; 
  //double integrand_psf = m_0 * m_plus * moment_pion_cross_sq; // Absence of amp_rhopi
  
  return integrand_psf;
    
}

double F3pi_num_fun(double *x, double *par) {

  double m = x[0];
  
  //cout << "m value = " << m << ", s_sqrt = " << s_sqrt << endl;
  
  double F3pi_value = 0.;

  const double m12_min = 2 * mpi_GeV; // GeV
  const double m12_max = m - mpi0_GeV; // GeV

  const double m23_max = m - mpi_GeV;
  const double m23_min = mpi_GeV + mpi0_GeV;

  const double m12sq_min = m12_min * m12_min;
  const double m12sq_max = m23_max * m23_max;

  //cout << "mpi_GeV = " << mpi_GeV << endl;

  const int np_m12 = 20;
  const int np_m23 = 20;

  int k = 0;
  double m12 = m12_min;
  double m23 = m23_min;
  double psf = 0.;

  double dm12 = (m12_max -m12_min) / np_m12;
  double dm23 = (m23_max -m23_min) / np_m23;
  
  //auto dt = new TGraph2D(np_m12 * np_m23);
    
  double dalitz_par[1] = {m};
  double dalitz_upp_tmp = 0.;
  double dalitz_down_tmp = 0.;

  const double norm_factor = 1. / 12. / pi / pi / m;

  //cout << "norm_factor = " << norm_factor << endl; 
  
  //dt->SetNpy(41);
  //dt->SetNpx(40);

  double m12sq = 0., m23sq = 0.; 
  double integrand_psf_par[2] = {0., m};
	
  for (Int_t i = 0; i < np_m12; i++) {
    for (Int_t j = 0; j < np_m23; j++) {

      m12sq = m12 * m12;
      m23sq = m23 * m23;
      
      dalitz_upp_tmp = dalitz_upp_fun(&m12, dalitz_par);
      dalitz_down_tmp = dalitz_down_fun(&m12, dalitz_par);

      //cout << "m12 = " << m12 << ", m12sq = " << m12sq << ", m23 = " << m23 << ", m23sq = " << m23sq << ", dalitz_upp_tmp = " << dalitz_upp_tmp << endl;
      if (m23sq <= dalitz_upp_tmp && m23sq >= dalitz_down_tmp && m12sq <= m12sq_max && m12sq >= m12sq_min){
	integrand_psf_par[0] = m23sq;
	psf = integrand_psf_fun(&m12sq, integrand_psf_par);
	//cout << psf << endl;
      }
      else {
	psf = 0.;
      }
      
      //dt -> SetPoint(k, m12sq, m23sq, psf);
      F3pi_value += (psf * dm12 * dm23);
      k++;
      m23 = m23 + dm23;
          
    }
    m12 = m12 + dm12;
    m23 = m23_min;
  }
  
  //gStyle->SetPalette(1);
  //dt->SetMarkerStyle(20);
  //dt->Draw("SURF3");

  return norm_factor * F3pi_value;
  
}




