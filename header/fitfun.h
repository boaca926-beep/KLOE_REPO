double gauss1d(double *x, double *par) {
  double r1 = (x[0] - par[1]) / par[2];
  double fitval = par[0] * TMath::Exp(-0.5 * r1 * r1);
  //cout << r1 << endl;
  return fitval;
}

double log_gauss1d(double *x, double *par) {
  double r1 = (TMath::Log(x[0]) - par[1]) / par[2];
  double fitval = par[0] * TMath::Exp(-0.5 * r1 * r1) / x[0];
  //cout << r1 << endl;
  return fitval;
}

double fun_double(double *x, double *par) {
  double *p1 = &par[0];
  double *p2 = &par[3];
  double result = gauss1d(x, p1) + gauss1d(x, p2);
  //double result = g1(x, p1);
  return result;
}

double fun_double_log(double *x, double *par) {
  double *p1 = &par[0];
  double *p2 = &par[3];
  double result = log_gauss1d(x, p1) + log_gauss1d(x, p2);
  //double result = g1(x, p1);
  return result;
}

double fun_pol(double *x,  double *par) {
  double p0 = par[0];
  double p1 = par[1];
  double p2 = par[2];

  //double f = 1. / (1. - p3 * exp((-x + p1) / p2));
  double f = p0 + p1 * x[0] + p2 * x[0] * x[0];
  //cout << "f = " << f << endl;
  
  return f;

}

Double_t ChiSqr7fitFunc(Double_t *x, Double_t *par) {
	double k=7./2.;
	double f1=ROOT::Math::tgamma(k); 
	double f2=TMath::Power(2,k); 
	double f3=TMath::Power(x[0],k-1.);
	double f4=TMath::Exp(-x[0]/(2.));
	double f5=par[0]*f3*f4;///(f1*f2);
	return f5;
}

Double_t ChiSqr8fitFunc(Double_t *x, Double_t *par) {
	double k=8./2.;
	double f1=ROOT::Math::tgamma(k); 
	double f2=TMath::Power(2,k); 
	double f3=TMath::Power(x[0],k-1.);
	double f4=TMath::Exp(-x[0]/(2.));
	double f5=par[0]*f3*f4;///(f1*f2);
	return f5;
}

Double_t ChiSqr7fitFunc_2par(Double_t *x, Double_t *par) {
	double k=par[0]/2.;
	double f1=ROOT::Math::tgamma(k); 
	double f2=TMath::Power(2,k); 
	double f3=TMath::Power(x[0],k-1.);
	double f4=TMath::Exp(-x[0]/(2.));
	double f5=par[1]*f3*f4;///(f1*f2);
	//double f5=par[1]*f1*f2;
	//cout<<par[0]<<endl;
	//std::cout<<"ndf = "<<par[0]<<", norm = "<<par[1]<<endl;
	return f5;
}
