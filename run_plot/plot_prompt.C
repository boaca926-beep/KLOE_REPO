#include "../header/plot.h"

int plot_prompt() {

  gErrorIgnoreLevel = kError;
  TGaxis::SetMaxDigits(3);
  gStyle->SetOptTitle(0);
  gStyle->SetStatBorderSize(0);
  gStyle->SetOptStat(0);

  // get histo
  TFile *intree = new TFile("/home/bo/Desktop/input_norm_TDATA/input/sig.root");

  TH1D *prompt_distr = (TH1D *)intree -> Get("prompt_distr");
  
  formatfill_h(prompt_distr, 4, 3001);

  /*
  TIter next_tree(intree -> GetListOfKeys());
  
  TString objnm_tree, classnm_tree;
  
  int i = 0;

  TKey *key;
  
  while ( (key = (TKey *) next_tree() ) ) {// start tree while lop
    
    i ++;
    
    objnm_tree   =  key -> GetName();
    classnm_tree = key -> GetClassName();
    key -> GetSeekKey();
    
    cout << " tree" << i << ": classnm = " << classnm_tree << ", objnm = " << objnm_tree << endl;
    
  }
  */

  // normalize histo
  const double bins = prompt_distr -> GetNbinsX(); 
  double entries  = prompt_distr -> Integral(1, bins);
  prompt_distr -> Scale(1. / entries);

  cout << "bins = " << bins << ", entries = " << entries << endl;

  TH1D * hprompt = new TH1D("hprompt", "", prompt_distr->GetNbinsX(), 0, prompt_distr->GetNbinsX());

  int k = 0;
  const double bar_width = 0.3 * (prompt_distr -> GetBinCenter(2) - prompt_distr -> GetBinCenter(1));
  double XC[10]; // X centers
  double YH[10]; // Y height
      
  cout << "bar_width = " << bar_width << endl;
  
  for (int i = 1; i <= prompt_distr->GetNbinsX(); i++) {
    XC[i-1] = i;
    YH[i-1] = prompt_distr -> GetBinContent(i);
	
    //hprompt -> Fill(i-1, prompt_distr -> GetBinContent(i));
    cout << i << ", bin content = " << XC[i-1] << ", bin center = " << YH[i-1] << endl;

    // Create a bar (TBox: x1, y1, x2, y2)
    TBox *bar = new TBox(
			 XC[i] - bar_width/2, // Left edge
			 (YH[i] > 0) ? 0 : YH[i], // Bottom (handles negative values)
			 XC[i] + bar_width/2, // Right edge
			 (YH[i] > 0) ? YH[i] : 0  // Top
			 );
    
    bar->SetFillColor(kBlue);
    bar->SetLineWidth(1);
    bar->Draw();
  }

  TGraph *gf_prompt = new TGraph(10, XC, YH);

  gf_prompt -> SetMarkerStyle(21);
  gf_prompt -> SetMarkerSize(1.0);
  gf_prompt -> SetLineColor(kBlack);
  gf_prompt -> SetLineWidth(2);
  gf_prompt -> SetMarkerColor(kBlack);
  gf_prompt -> GetXaxis() -> CenterTitle();
  //gf_prompt -> GetXaxis() -> SetRangeUser(-0.5., 0.5.);
  gf_prompt -> GetXaxis() -> SetTitleOffset(1.);
  gf_prompt -> GetXaxis() -> SetTitleSize(0.06);
  gf_prompt -> GetXaxis() -> SetLabelSize(0.05);  
  gf_prompt -> GetXaxis() -> SetNdivisions(20);
  gf_prompt -> GetXaxis() -> SetTitle("Number of prompt photons n_{#gamma}");
  
  gf_prompt -> GetYaxis() -> SetTitleOffset(.8);
  gf_prompt -> GetYaxis() -> SetTitleSize(0.08);
  gf_prompt -> GetYaxis() -> SetLabelSize(0.05);  
  gf_prompt -> GetYaxis() -> CenterTitle();
  gf_prompt -> GetYaxis() -> SetTitle("Probabilities");
  
  //formatfill_h(gf_prompt, 4, 3001);

  // plots
  TCanvas *cv = new TCanvas("prompt_distr", "", 0, 0, 1000, 800);
  cv -> SetBottomMargin(0.15);//0.007
  cv -> SetLeftMargin(0.15);

  /*
  gf_prompt -> SetMarkerStyle(21);
  gf_prompt -> SetMarkerSize(0.7);
  gf_prompt -> GetYaxis() -> SetTitle("Probabilities");
  gf_prompt -> GetYaxis() -> SetTitleOffset(1.);
  gf_prompt -> GetYaxis() -> SetTitleSize(0.07);
  gf_prompt -> GetYaxis() -> SetLabelSize(0.05);
  gf_prompt -> GetYaxis() -> CenterTitle();
  
  //gf_prompt -> GetYaxis() -> SetRangeUser(0.1, ymax);
  gf_prompt -> GetXaxis() -> SetTitle("Number of prompt photons n_{#gamma}");
  gf_prompt -> GetXaxis() -> SetRangeUser(0, 10);
  gf_prompt -> GetXaxis() -> SetTitleOffset(0.9);
  gf_prompt -> GetXaxis() -> SetTitleSize(0.07);
  gf_prompt -> GetXaxis() -> SetLabelSize(0.05);
  gf_prompt -> GetXaxis() -> CenterTitle();
  */
  
  gf_prompt -> Draw("AP");
  gf_prompt -> SetMaximum(10.);
  
  gPad -> SetLogy();

  cv -> Modified();
  cv -> Update();
    
  cv -> SaveAs("prompt.pdf");
  

  return 0;
  
}


