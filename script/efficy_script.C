void efficy_script() {
  gROOT->ProcessLine(".L ../run_ufo/efficy_evtcls.C");
  gROOT->ProcessLine("efficy_evtcls()");
}
