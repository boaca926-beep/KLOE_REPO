void run_script() {
  gROOT->ProcessLine(".L ../run_vertex/MyClass.C");
  gROOT->ProcessLine(".L ../run/Analys_class.C");
  gROOT->ProcessLine("Analys_class(rootFile, sampleFile)");
}
