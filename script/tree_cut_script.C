#include <iostream>
void tree_cut_script() {
gROOT->ProcessLine(".L ../run_ufo/tree_cut_bkgrej.C");
gROOT->ProcessLine("tree_cut_bkgrej()");
}
