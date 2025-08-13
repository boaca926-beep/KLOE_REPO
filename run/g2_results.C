int g2_results(){

  const double bias = 1165900;

  // sm prediction
  const double amu_sm = 116592033 * 1e-11;
  const double amu_sm_err = 62 * 1e-11;

  const double amu_FNAL_21 = 116592040 * 1e-11;
  const double amu_err_FNAL_21 = TMath::Sqrt(43 * 43 + 15 * 15 + 3 * 3) * 1e-11;

  // Experimental evarage
  const double amu_exp = 116592071.5 * 1e-11;
  const double amu_exp_err = 14.5 * 1e-11;

  double diff_FNAL_21 = (amu_FNAL_21 - amu_sm) * 1e10;
  double diff_err_FNAL_21 = TMath::Sqrt(amu_err_FNAL_21 * amu_err_FNAL_21 + amu_sm_err * amu_sm_err) * 1e10;

  double Delta_amu = (amu_exp - amu_sm);
  double Delta_amu_err = TMath::Sqrt(amu_exp_err * amu_exp_err + amu_sm_err * amu_sm_err);
  
  cout << "bias = " << bias << ", "
       << "amu_FNAL_21 = " << amu_FNAL_21 << ", diff (10^-10) = " << diff_FNAL_21 << "+/-" << diff_err_FNAL_21 << "\n"
       << "Delta amu (X10^10)= " << Delta_amu * 1e10 << "+/-" << Delta_amu_err * 1e10 << endl;
  
  return 0;
  
}
