#include "../header/syst_summary.h"
#include "../header/width_summary.h"
//#include "bb_summary.h"


int syst_summary() {

  cout << "Total systematic uncertainty" << endl;
  
  // comment part
  double systerr_plus=0., systerr_plus_sq_sum = 0.;
  double systerr_neg=0., systerr_neg_sq_sum = 0.;

  cout << "para_nm = " << para_nm << endl;
  cout << "syst. err. = (err_plus, err_minus), (err2_sum_plus, err2_sum_neg) \n" << endl;

  for (int i=0; i < list_size; i++) {

    systerr_plus_sq_sum += SYSTERR_plus[i] * SYSTERR_plus[i];
    systerr_neg_sq_sum += SYSTERR_neg[i] * SYSTERR_neg[i];

    //Form(%.2f", variable)
    cout << i + 1 << ". " << CUT_NM[i] << ": (" << Form("%.2f", SYSTERR_plus[i]) << ", " << Form("%.2f", SYSTERR_neg[i]) << ")\n";
    
    //cout << i + 1 << ". " << CUT_NM[i] << ": syst. err = (" << SYSTERR_plus[i] << ", " << SYSTERR_neg[i] << ")\n";
      //<< "\n\t err2_sum = (" << systerr_plus_sq_sum << ", " << systerr_neg_sq_sum << ")\n\n";
  }

  systerr_plus=TMath::Sqrt(systerr_plus_sq_sum);
  systerr_neg=TMath::Sqrt(systerr_neg_sq_sum);
    
  cout << "\nPara. " << para_nm << " = " << value_nomi << "\n"
       << "Stat. Err. = +/-" << error_stat << "\n"
       << "Syst. Err. = (+" << Form("%.2f", systerr_plus) << ", -" << Form("%.2f", systerr_neg) << ")" << endl;
  
  return 0;
}
