// error propagation for FILFO efficiency
int error_propa() {

  /*
  const double x = 397856.;
  double dx = TMath::Sqrt(x);
  
  const double y = 44040.;
  double dy = TMath::Sqrt(y);
  */

  // test values of x and y
  const double x = 5.0, dx = 0.1;
  const double y = 3.0, dy = 0.2;
  
  double z = 2 * x / (y + 2 * x);
  double dz_dx = 2 * y / TMath::Power(y + 2 * x, 2);
  double dz_dy = -2 * x / TMath::Power(y + 2 * x, 2);

  double z_err2 = TMath::Power(dz_dx * dx, 2) + TMath::Power(dz_dy * dy, 2);
  double z_err = TMath::Sqrt(z_err2);
  
  cout << "x = " << x << "+/-" << dx << "\n"
       << "y = " << y << "+/-" << dy << "\n"
       << "z = " << z << "+/-" << z_err << endl;

  

  return 0;
}
