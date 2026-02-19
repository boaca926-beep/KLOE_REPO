#include <iostream>
void plot_script() {
  gROOT->ProcessLine(".L ../run_plot/plot_omega.C");
  gROOT->ProcessLine("plot_omega()");
}
