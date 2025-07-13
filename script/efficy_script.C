void efficy_script() {
  gROOT->ProcessLine(".L ../run_ufo/efficy_filfo.C");
  gROOT->ProcessLine("efficy_filfo()");
}
