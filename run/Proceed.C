//void Proceed(const TString input_str = "Analys_class(\"../../path_small/ksl_path1\",\"../../ksl_pre\")"){
//void Proceed(const TString input_str = "Analys_class(\"../path_chain/SIG_sum\",\"../../sig_sum\")"){
//void Proceed(const TString input_str = "Analys_class(\"../../path_norm/sig_path_sum\",\"../../sig_sum\")"){
//void Proceed(const TString input_str = "Analys_class(\"../../path_small/sig_path\",\"../../sig\")"){
//void Proceed(const TString input_str = "Analys_class(\"../../path_small/sig_tmp\",\"../../sig_tmp1\")"){
//void Proceed(const TString input_str = "Analys_class(\"../path_chain/exp_path\",\"../../exp\")"){
void Proceed(const TString input_str = "Analys_class(\"../path_chain/ufo_path\",\"../../ufo\")"){

  cout << "READING INPUT FILES !" << endl;
  //gROOT -> ProcessLine( "gErrorIgnoreLevel = 1001;")
  //gROOT -> ProcessLine(".L ../run/MyClass.C");
  //gROOT -> ProcessLine(".L ../run/Analys_class.C");

  //gROOT -> ProcessLine(".L ../run_vertex/MyClass.C");
  //gROOT -> ProcessLine(".L ../run/Analys_class.C");

  gROOT -> ProcessLine(".L ../run_ufo/MyClass.C");
  gROOT -> ProcessLine(".L ../run/Analys_class.C");
  
  gROOT -> ProcessLine(input_str);
}
