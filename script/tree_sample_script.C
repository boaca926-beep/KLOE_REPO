#include <iostream>
void tree_sample_script() {
gROOT->ProcessLine(".L ../run_vertex_bkgrej/tree_sample.C");
gROOT->ProcessLine("tree_sample()");
}
