#include <iostream>
void compr_script() {
  gROOT->ProcessLine(".L ../run/pulls_compr.C");
  gROOT->ProcessLine("pulls_compr()");
}
