TFile *f_in = new TFile("result.root");

TH1D *hsig = (TH1D *) f_in -> Get("hsig");
TH1D *hsig_true = (TH1D *) f_in -> Get("hsig_true");
TH1D *hsig_gen = (TH1D *) f_in -> Get("hsig_gen");
TH1D *hefficy = (TH1D *) f_in -> Get("hefficy");

TH2D *hsmearmatr = (TH2D *) f_in -> Get("hsmearmatr");
TH2D *hsmearmatr_err = (TH2D *) f_in -> Get("hsmearmatr_err");
TH2D *hsmearmatr_efficy = (TH2D *) f_in -> Get("hsmearmatr_efficy");
TH2D *hsmearmatr_efficy_err = (TH2D *) f_in -> Get("hsmearmatr_efficy_err");
TH2D *hcorrmatrix = (TH2D *) f_in -> Get("hcorrmatrix");

TH1D* hcorrmatrix_projY = hcorrmatrix -> ProjectionY("hcorrmatrix_projY");
hcorrmatrix_projY -> SetName("hcorrmatrix_projY");

TTree *TINPUT = (TTree *) f_in -> Get("TINPUT");

  
const int list_size = 1000;

double M3PI[list_size], M3PI_ERR[list_size]; 
double EFFICY[list_size]; 
double EFFICY_ERR[list_size]; 
double ISRLUMI[list_size]; 
double CRX3PI_BW[list_size]; 

double N3PI_PRE[list_size];
double N3PI_SMEAR[list_size], N3PI_SMEAR_ERR[list_size];

double get_sigma2_y(double a, double b, double a_err, double b_err) {

  double term1 = TMath::Power(a_err * b_err, 2); // Var(a)V(b)
  double term2 = TMath::Power(a * b_err, 2); // (Ea)^{2}Var(b)
  double term3 = TMath::Power(b * a_err, 2); // (Eb)^{2}Var(a)
  
  double sigma2 = term1 + term2 + term3;
  
  return sigma2;
}

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
