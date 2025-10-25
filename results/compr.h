const TString file_type = "vmd";
//const TString infile_tmp = "../crx3pi/output_" + file_type + "/crx3pi0.root";
const TString infile_tmp = "/media/bo/Backup/KLOE_OUTPUT/output_" + file_type + "/crx3pi0.root";

const TString model_type = "VMD";

void SetGFAttr(TGraph *gf, const TString x_title = "", const TString y_title = "Lumi_{ISR} [nb^{-1}]") {

  //gf -> SetTitle("systmatics");
  gf -> SetName("gf");
  gf -> SetMarkerStyle(20);
  gf -> SetMarkerSize(0.6);
  //gf -> SetLineColor(0);
  gf -> SetMarkerColor(kBlack);
  gf -> GetXaxis() -> SetRangeUser(750., 820.);
  gf -> GetXaxis() -> SetTitle(x_title);
  gf -> GetXaxis() -> SetLabelSize(0.04);  
  gf -> GetXaxis() -> SetTitleOffset(.8);
  gf -> GetXaxis() -> SetTitleSize(.045);
  gf -> GetXaxis() -> CenterTitle();

  gf -> GetYaxis() -> SetTitle(y_title);
  gf -> GetYaxis() -> CenterTitle();
  gf -> GetYaxis() -> SetLabelSize(0.04);  
  gf -> GetYaxis() -> SetTitleOffset(.7);
  gf -> GetYaxis() -> SetTitleSize(.045);
  //gf -> GetYaxis() -> SetRangeUser(y_min, y_max);
  
}
