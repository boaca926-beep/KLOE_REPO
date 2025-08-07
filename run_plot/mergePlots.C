#include "../header/plot_omega_syst.h"

int mergePlots() {

  gErrorIgnoreLevel = kError;
  TGaxis::SetMaxDigits(3);
  //gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);  

  const int para_nb = 3;
  
  TCanvas *cv = new TCanvas("cv", "", 0, 0, 1000, 1500);
  cv -> Divide(1, para_nb, 0, 0);

  const TString para_str[para_nb] = {"omega_mass", "omega_width", "BB"};
  const TString para_title[para_nb] = {"M_{#omega} [MeV/c^{2}]", "#Gamma_{#omega} [MeV]", "B_{ee}B_{3#pi}"};

  int indx = 0;
  TString finput_tmp = "";
  
  for (Int_t i = 0; i < 1; i ++) {// column index

    for (Int_t j = 0; j < para_nb; j ++) {// row index 

      finput_tmp = outputPlot + para_str[indx] + "_" + cut_label + ".root";
      TFile* file_input = new TFile(finput_tmp);
      
      cout << para_str[indx] << " from " << finput_tmp << endl;

      TIter next_tree(file_input -> GetListOfKeys());
      TString objnm_tree, classnm_tree;
      
      //int k = 0;
      TKey *key;
      
      while ( (key = (TKey *) next_tree() ) ) {// start cutlist_tree while loop
	
	objnm_tree   = key -> GetName();
	classnm_tree = key -> GetClassName();
	key -> GetSeekKey();
	
	cout << "classnm = " << classnm_tree << ", objnm = " << objnm_tree << endl;
	
      }

      indx ++;
      
      // plot
      cv -> cd(indx);

      gPad -> SetLeftMargin(0.15);
      gPad -> SetRightMargin(0.017);
      gPad -> SetBottomMargin(0.07);
      gPad -> SetTopMargin(0.02);
	

      TGraphErrors *gf_Z = (TGraphErrors *) file_input -> Get("gf_Z");
      TGraphErrors *gf_Z_norm = (TGraphErrors *) file_input -> Get("gf_Z_norm");
      TGraphErrors *gf_mg1 = (TGraphErrors *) file_input -> Get("gf_mg1");

      TPaveText *pt = new TPaveText(0.7, 0.9, 0.7, 0.9, "NDC");
  
      pt -> SetTextSize(0.06);
      pt -> SetFillColor(0);
      pt -> SetTextAlign(12);
      pt -> SetTextColor(kBlack);
      pt -> AddText(para_title[indx - 1]);

      if (indx == 3) {
	cout << "here" << endl;
	//gf_Z -> GetXaxis() -> SetLabelOffset(.5);
	gPad -> SetBottomMargin(0.25);
      }
      else {
	gf_Z -> GetXaxis() -> SetLabelOffset(.5);
      }
      
      gf_Z -> GetYaxis() -> SetTitleSize(0.1);
      gf_Z -> GetYaxis() -> SetTitleOffset(.7);
      gf_Z -> GetYaxis() -> SetLabelSize(0.1);

      gf_Z -> GetXaxis() -> SetTitleSize(0.1);
      gf_Z -> GetXaxis() -> SetTitleOffset(1.0);
      gf_Z -> GetXaxis() -> SetLabelSize(0.07);
      
      gf_Z -> SetMarkerColorAlpha(kWhite, 0.6);
      gf_Z -> Draw("AP");
      gf_Z_norm -> Draw("P");
      //gf_mg1 -> Draw("C3");
      //pt -> Draw("Same");

      const int list_size = 3; //gf_Z -> GetN();
      const double bar_width = 0.1 * (gf_Z -> GetX()[1] - gf_Z -> GetX()[0]);
      double XC[list_size]; // X centers
      double YH[list_size]; // Y height
      cout << "list_size = " << list_size << ", bar_width = " << bar_width << endl;
      
      for (int i = 0; i < gf_Z -> GetN(); i++) {
	XC[i] = gf_Z -> GetX()[i];
	YH[i] = gf_Z -> GetY()[i];
	
	cout << XC[i] << ", " << YH[i] << endl;
	//hgf_Z -> SetBinContent(i+1, gf_Z -> GetY()[i]);
	
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
  
      gf_mg1 -> Draw("C3");
      pt -> Draw("Same");
      gPad -> SetGrid(); // Add grid lines
      
      file_input -> Close();
    }

  }

  cv -> SaveAs(outputPlot + cut_label + "_merged.pdf");

  return 0;
  
}
