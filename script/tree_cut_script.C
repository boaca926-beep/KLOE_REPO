#include <iostream>
void tree_cut_script() {
gROOT->ProcessLine(".L ../run_ufo/tree_cut_evtcls.C");
gROOT->ProcessLine("tree_cut_evtcls()");
}
