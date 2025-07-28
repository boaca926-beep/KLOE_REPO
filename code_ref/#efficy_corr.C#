#include "../header/sm_para.h"
#include "../header/path.h"
#include "../header/sfw2d.txt"
#include "../header/sfw1d.txt" 
#include "../header/method.h"
#include "../header/omega_fit_efficy.h"
#include "../header/cut_para.h"


int efficy_corr () {

  TFile *f_efficy_ratio = new TFile("../../efficy_evtcls/efficy_ratio.root");

  getObj(f_efficy_ratio);
  TGraphErrors* gf_ratio = (TGraphErrors*)f_efficy_ratio -> Get("gf_ratio");
  
  get_efficy(gf_ratio);
  

  return 0;
  
}
