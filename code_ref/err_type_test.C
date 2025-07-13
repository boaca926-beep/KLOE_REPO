#include "../header/plot_omega_syst.h"
#include "err_type_test.h"

int err_type_test(){

  // different signs
  //double err1 = -3, err2 = 1; // case I
  //double err1 = 20, err2 = -1; // case II

  // same sings
  //double err1 = -1, err2 = -10; // case I
  //double err1 = 2, err2 = 8; // case II

  // at least one 0
  //double err1 = 0., err2 = 19.;
  //double err1 = -10., err2 = 0.;
  double err1 = 0., err2 = 0.;
  
  get_syst_errII(err1, err2);

  cout << "nega_err = " << SYST_ERRII[0] << ", plus_err = " << SYST_ERRII[1] << endl;

  return 0;
  
}
