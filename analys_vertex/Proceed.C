void Proceed(const TString input_str = "Analys_class(\"../../path_small/sig_tmp\",\"../../sig_tmp\")"){
//void Proceed(const TString input_str = "Analys_class(\"../../path_small/sig_path1\",\"../../sig_sum\")"){

  cout << "READING INPUT FILES !" << endl;
  //gROOT -> ProcessLine( "gErrorIgnoreLevel = 1001;")
  gROOT -> ProcessLine(".L MyClass.C");
  gROOT -> ProcessLine(".L Analys_class.C");
  gROOT -> ProcessLine(input_str);
}
