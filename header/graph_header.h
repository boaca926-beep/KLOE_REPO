TGraphErrors *get_graph_syst(double xlist[], double ylist[], double xlist_err[], double ylist_err[], int length) {

  TGraphErrors *gf = new TGraphErrors(length, xlist, ylist, xlist_err, ylist_err);

  double upper_band_size = 7.;
  double down_band_size = 5.;

  
  gf -> SetMarkerStyle(20);
  gf -> SetMarkerSize(1.2);
  gf -> SetLineColor(kBlack);
  gf -> SetLineWidth(2);
  gf -> SetMarkerColor(kBlack);
  gf -> GetXaxis() -> SetTitle(cut_title);
  gf -> GetXaxis() -> CenterTitle();
  //gf -> GetXaxis() -> SetRangeUser(-0.5., 0.5.);
  gf -> GetXaxis() -> SetTitleOffset(1.);
  gf -> GetXaxis() -> SetTitleSize(0.06);
  gf -> GetXaxis() -> SetLabelSize(0.05);  

  gf -> GetYaxis() -> SetTitleOffset(1.4);
  gf -> GetYaxis() -> SetRangeUser(Y_NORM[0] - down_band_size * Y_ERR_NORM[0], Y_NORM[0] + upper_band_size * Y_ERR_NORM[0]);
  gf -> GetYaxis() -> SetTitleSize(0.06);
  gf -> GetYaxis() -> SetLabelSize(0.05);  
  gf -> GetYaxis() -> SetTitle(para_title + " " + para_unit);
  gf -> GetYaxis() -> CenterTitle();
  
  //gf -> Draw("AP");
  
  return gf;

}
