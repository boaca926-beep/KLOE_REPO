#include <iostream>
void tree_cut_script() {
gROOT->ProcessLine(".L ../run_vertex_bkgrej/tree_cut_bkgrej.C");
gROOT->ProcessLine("tree_cut_bkgrej()");
}
