double SYST_ERRII[2] = {0., 0.}; // list[0]: nega_err, list[1]: plus_err

void get_syst_errII(double err1, double err2){

  double max = getmaxoftwo(err1, err2);
  double min = getminoftwo(err1, err2);
  double plus_err = 0., nega_err = 0.;

    
  if (err1 * err2 < 0.) {

    if (err1 == max) {
      plus_err = err1;
      nega_err = err2;
    }
    else {
      plus_err = err2;
      nega_err = err1;
    }
    //cout << max << endl;
    
  }
  else if (err1 * err2 > 0.) {

    //cout << "min = " << min << ", max = " << max << endl;
    
    if (err1 < 0) {
      plus_err = 0.;
      nega_err = min;
    }
    else {
      plus_err = max;
      nega_err = 0.;
    }
    
  }
  else if (err1 == 0. && err2 != 0.) {

    if (err2 > 0.) {
       plus_err = err2;
       nega_err = 0.;
    }
    else {
      plus_err = 0.;
      nega_err = err2;
    }
    
  }
  else if (err2 == 0. && err1 != 0.) {

    if (err1 > 0.) {
       plus_err = err1;
       nega_err = 0.;
    }
    else {
      plus_err = 0.;
      nega_err = err1;
    }
    
  }
  else if (err1 == 0. && err2 == 0.) {

    plus_err = 0.;
    nega_err = 0.;
    
  }
  else {
    cout << "New case !!!" << endl;
  }
  
  
  SYST_ERRII[0] = nega_err;
  SYST_ERRII[1] = plus_err;
  
  

}
