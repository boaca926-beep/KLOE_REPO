#include <iostream>
#include "tree.h"
void run_script() {
  gROOT->ProcessLine(".L MyClass.C");
  gROOT->ProcessLine(".L Analys.C");
  gROOT->ProcessLine("Analys(inputFile, outputFile)");
}
