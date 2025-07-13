#include "../header/path.h"
#include "../header/hist.h"

int gethist(){

  // Preparing histos ...
  fillHist();

  /// save histos
  TFile * f_output = new TFile(outputHist + ".root", "recreate");

  HIM3pi_fit -> Write("HIM3pi_fit", 1);
  HSFW2D -> Write("HSFW2D", 1);
  HSFW1D -> Write("HSFW1D", 1);
  HSIG -> Write("HSIG",1);  
  HIM3pi_crx -> Write("HIM3pi_crx",1);

  /*
  HIM3pi_fit -> Write();
  HSFW2D -> Write();
  HSFW1D -> Write();
  HSIG -> Write();  
  HIM3pi_crx -> Write();
  */
  
  f_output -> Close();
  
  return 0;
  
}
