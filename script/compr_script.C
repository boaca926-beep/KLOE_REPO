#include <iostream>
void compr_script() {
  gROOT->ProcessLine(".L ../run/bkg_compr.C");
  gROOT->ProcessLine("bkg_compr()");
}
