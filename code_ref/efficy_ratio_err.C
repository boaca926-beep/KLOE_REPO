double get_efficy_ratio_err(double *PARA, double *PARA_ERR, double x) {

  double ratio_err = 0.;
  double ratio_err2 = 0.;

  double a_err = PARA_ERR[0];
  double b_err = PARA_ERR[1];
  double c_err = PARA_ERR[2];
  
  double term1 = TMath::Power(x, 4) * TMath::Power(a_err, 2);
  double term2 = TMath::Power(x, 2) * TMath::Power(b_err, 2);
  double term3 = TMath::Power(c_err, 2);

  ratio_err2 = term1 + term2 + term3;
  ratio_err = TMath::Sqrt(ratio_err2);
  
  cout << "a = " << PARA[0] << "+/-" << PARA_ERR[0] << "\n"
       << "b = " << PARA[1] << "+/-" << PARA_ERR[1] << "\n"
       << "c = " << PARA[2] << "+/-" << PARA_ERR[2] << "\n"
       << "x = " << x << "\n"
       << "ratio_err = " << ratio_err << endl;
 
  return ratio_err;
}

int efficy_ratio_err() {

  double a = 2., a_err = 0.1;
  double b = 3., b_err = 0.05;
  double c = 1., c_err = 0.01;

  double PARA[3] = {2., 3., 1.};
  double PARA_ERR[3] = {0.1, 0.05, 0.01};

  double x = 4.;
  
  double ratio_err = get_efficy_ratio_err(PARA, PARA_ERR, x);
  
       

  return 0;
  
}
