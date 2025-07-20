#include "../header/sm_para.h"
#include "../header/path.h"
#include "../header/sfw2d.txt"
#include "../header/sfw1d.txt"
#include "../header/psf_fun.h"
#include "../header/amp_fun.h"
#include "../header/model.h"
#include "../header/method.h"
//#include "../header/omega_fit_vmd.h"
//#include "../header/cut_para.h"

int peak_crx3pi() {

  double CRX3PI_LIST[1] = {1635};
  double MASS_LIST[1] = {782.73};
  double BB_LIST[1] = {7.81E-5};
  
  double BB = BB_LIST[0];
  //double BB = GetBB_new(CRX3PI_LIST[0], MASS_LIST[0]);
  double crx3pi_peak = GetCrx3piMax_new(BB_LIST[0], MASS_LIST[0]);

    
  cout << "BB = " << BB_LIST[0] << ", crx3pi_peak = " << crx3pi_peak << endl;
  

  return BB;

}
