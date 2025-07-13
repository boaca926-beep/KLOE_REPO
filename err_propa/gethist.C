#include "path.h"
#include "hist.h"

int gethist(){

  cout << "IM3pi_bin = " << IM3pi_bin << endl;
  
  fillHist();

  TFile * f_output = new TFile(outputHist + "hist.root", "recreate");

  /// save histos
  HIM3pi_fit -> Write("HIM3pi_fit", 1);
  HSFW2D -> Write("HSFW2D", 1);
  HSFW1D -> Write("HSFW1D", 1);
  HSIG -> Write("HSIG",1);  
  HIM3pi_crx -> Write("HIM3pi_crx",1);
    
  f_output -> Close();
  
  return 0;
  
}
