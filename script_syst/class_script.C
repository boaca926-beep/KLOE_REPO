void class_script() {
  gROOT->ProcessLine(".L ../run_ufo/MyClass.C");
  gROOT->ProcessLine(".L ../run/Analys_class.C");
  gROOT->ProcessLine("Analys_class(rootFile, sampleFile)");
}
