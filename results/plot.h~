void SetGFAttr(TGraph * gf_fcn, const TString x_title = "", const TString y_title = "Lumi_{ISR} [nb^{-1}]") {

  //gf_fcn -> SetTitle("systmatics");
  gf_fcn -> SetName("gf");
  gf_fcn -> SetMarkerStyle(20);
  gf_fcn -> SetMarkerSize(0.6);
  //gf_fcn -> SetLineColor(0);
  gf_fcn -> SetMarkerColor(kBlack);
  gf_fcn -> GetXaxis() -> SetRangeUser(750., 820.);
  gf_fcn -> GetXaxis() -> SetTitle(x_title);
  gf_fcn -> GetXaxis() -> SetLabelSize(0.04);  
  gf_fcn -> GetXaxis() -> SetTitleOffset(.8);
  gf_fcn -> GetXaxis() -> SetTitleSize(.045);
  gf_fcn -> GetXaxis() -> CenterTitle();

  gf_fcn -> GetYaxis() -> SetTitle(y_title);
  gf_fcn -> GetYaxis() -> CenterTitle();
  gf_fcn -> GetYaxis() -> SetLabelSize(0.04);  
  gf_fcn -> GetYaxis() -> SetTitleOffset(.7);
  gf_fcn -> GetYaxis() -> SetTitleSize(.045);
  //gf_fcn -> GetYaxis() -> SetRangeUser(y_min, y_max);
  
}

void PteAttr(TPaveText *pt) {

  pt -> SetTextSize(0.04);
  pt -> SetFillColor(0);
  pt -> SetTextAlign(12);

}

void legtextsize(TLegend* l, Double_t size) {
  for(int i = 0 ; i < l -> GetListOfPrimitives() -> GetSize() ; i++) {
    TLegendEntry *header = (TLegendEntry*)l->GetListOfPrimitives()->At(i);
    header->SetTextSize(size);
  }
}

void SetLegend(TLegend* l) {

  l -> SetTextFont(132);
  l -> SetFillStyle(0);
  l -> SetBorderSize(0);
  l -> SetNColumns(2);

}

void SetCVAttr(TH1D *h1d, TString var_symb, TString unit) {

  const double ymax_normed = h1d -> GetMaximum();
  cout << "ymax = " << ymax_normed << endl;
  
  h1d -> GetXaxis() -> SetTitle(var_symb + " " + unit);
  h1d -> GetXaxis() -> CenterTitle();
  h1d -> GetXaxis() -> SetTitleSize(0.07);
  h1d -> GetXaxis() -> SetTitleOffset(.9);
  h1d -> GetXaxis() -> SetLabelOffset(.01);
  h1d -> GetXaxis() -> SetLabelSize(0.05);//0.03
  //h1d -> GetXaxis() -> SetRangeUser(110., 450.);
  //h1d -> GetXaxis() -> SetRangeUser(-500., 50.);
  //h1d -> GetXaxis() -> SetRangeUser(20., 180.);
  
  //h1d -> GetYaxis() -> SetLabelFont(43); // Absolute front size in pixel (precision 3)
  h1d -> GetYaxis() -> SetLabelSize(0.05);
  //h1d -> GetYaxis() -> SetRangeUser(1e-8, ymax_normed * 60.); //*1e2, *1.3, *30 *60
  h1d -> GetYaxis() -> SetRangeUser(1e-5, ymax_normed * 30.); //*1e2, *1.3, *30 *60, *1.5
  
  //h1d -> GetYaxis() -> SetTitle(TString::Format("Entries/%0.2f", binwidth) + " " + unit);
  //h1d -> GetYaxis() -> SetTitle("Events");
  h1d -> GetYaxis() -> CenterTitle();
  h1d -> GetYaxis() -> SetTitleSize(0.07);
  h1d -> GetYaxis() -> SetTitleOffset(1.1);
  
  //hist_data -> SetStats(0);
  
}
