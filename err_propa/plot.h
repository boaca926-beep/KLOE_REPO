void format_h(TH1D* h, Int_t linecolor, Int_t width) {
  h->SetLineColor(linecolor);
  //cout << "histo format" << endl;
  h->SetLineWidth(width);
}

void formatfill_h(TH1D* h, Int_t fillcolor, Int_t fillstyle) {
  h -> SetFillStyle(fillstyle);
  h -> SetFillColor(fillcolor);
  h -> SetLineColor(0);
}

double getbinwidth(TH1D* h) {
  Int_t binsize=0;
  double width=0.;
  double xmax=0., xmin=0.;
  xmax = h->GetXaxis()->GetXmax(); //cout<<xmax<<endl;
  xmin = h->GetXaxis()->GetXmin(); //cout<<xmin<<endl;
  binsize=h->GetNbinsX(); //cout<<"binsize = " << binsize << endl;
  width=(xmax-xmin)/binsize; //cout<<width<<endl;

  return width;
}

void SetLegend(TLegend* l) {

  l -> SetTextFont(132);
  l -> SetFillStyle(0);
  l -> SetBorderSize(0);
  l -> SetNColumns(2);

}

void legtextsize(TLegend* l, Double_t size) {
  for(int i = 0 ; i < l -> GetListOfPrimitives() -> GetSize() ; i++) {
    TLegendEntry *header = (TLegendEntry*)l->GetListOfPrimitives()->At(i);
    header->SetTextSize(size);
  }
}

//
void PteAttr(TPaveText *pt) {

  pt -> SetTextSize(0.04);
  pt -> SetFillColor(0);
  pt -> SetTextAlign(12);

}

