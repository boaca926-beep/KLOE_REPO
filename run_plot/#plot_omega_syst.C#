//#include "../header/path.h"
#include "../header/sm_para.h"
#include "../header/plot_omega_syst.h"
#include "../header/method.h"
#include "../header/graph.h"

int plot_omega_syst() {

  gErrorIgnoreLevel = kError;
  TGaxis::SetMaxDigits(4);
  //gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);  

  // loop over input files
  string line;
  int f_indx = 0;

  //double w_y = 0.;
  //double w_sum = 0.;
  //double wy_sum = 0.;
  
  ifstream filelist(syst_path + "path_output.txt");
  cout << syst_path << endl;
  //ifstream filelist("../../result_syst/norm_lumi_nb/path_output.txt");
  //ifstream filelist("../../presel_syst_egammamin/path_output.txt");

  //TFile* f_tmp = new TFile("../../presel_syst_egammamin/result0/omega_fit/omegafit.root");
  cout << syst_path + "path_output.txt" << endl;
  
  if (filelist.is_open()) {
    while (!filelist.eof()) {
      if (getline(filelist, line, '\n')) {
	if (line[0] != '!') {
	  TFile* f_in = new TFile(line.data());
	  TString fname(line.data());
	  //getObj(f_in);
	  //cout << line.data() << " !!!" << endl;
	
	  
	  double cut_value[1];
	  double x = 0., x_err = 0.;
	  
	  double OMEGA_PARA[3]; // sfw2d MC scaling factors
	  double OMEGA_PARA_ERR[3];
	  
	  TTree* tree_info = (TTree*)f_in -> Get("TCRX3PI");
	  tree_info -> SetBranchAddress("Br_OMEGA_PARA", OMEGA_PARA);
	  tree_info -> SetBranchAddress("Br_OMEGA_PARA_ERR", OMEGA_PARA_ERR);
	  tree_info -> SetBranchAddress("Br_cut_value", cut_value);

	  // Branches
	  int entries = tree_info -> GetEntries();
	  for (Int_t irow = 0; irow < entries; irow++) {// tree-loop
	    
	    tree_info -> GetEntry(irow);
	    //entries = tree_info -> GetLeaf("Br_entries") -> GetValue(0);
	    
	  }// end tree-loop

	  //isr3pi_sfw = OMEGA_PARA[1];
	  //isr3pi_sfw_err = OMEGA_PARA[1];

	  //omegapi_sfw = OMEGA_PARA[2];
	  //omegapi_sfw_err = OMEGA_PARA[2];

	  //etagam_sfw = OMEGA_PARA[3];
	  //etagam_sfw_err = OMEGA_PARA[3];

	  //ksl_sfw = OMEGA_PARA[4];
	  //ksl_sfw_err = OMEGA_PARA[4];

	  //mcrest_sfw = OMEGA_PARA[5];
	  //mcrest_sfw_err = OMEGA_PARA[5];
	  
	  //<<"\tisr3pi_sfw = "<< isr3pi_sfw << "+/-" << isr3pi_sfw_err << "\n"
	  //<<"\tomegapi_sfw = "<< omegapi_sfw << "+/-" << omegapi_sfw_err << "\n"
	  //<<"\tetagam_sfw = "<< etagam_sfw << "+/-" << etagam_sfw_err << "\n"
	  //<<"\tksl_sfw = "<< ksl_sfw << "+/-" << ksl_sfw_err << "\n"
	  //<<"\tmcrest_sfw = "<< mcrest_sfw << "+/-" << mcrest_sfw_err << "\n";
	       
	  // Fill list
	  XLIST[f_indx] = cut_value[0];
	  XLIST_ERR[f_indx] = 0;

	  if (para_indx == 2) {
	    YLIST[f_indx] = OMEGA_PARA[para_indx] * 1e5;
	    YLIST_ERR[f_indx] = OMEGA_PARA_ERR[para_indx] * 1e5;
	  }
	  else {
	    YLIST[f_indx] = OMEGA_PARA[para_indx];
	    YLIST_ERR[f_indx] = OMEGA_PARA_ERR[para_indx];
	  }
	  
	  //cout << f_indx << ", " << line.data() << ", " << cut_label << " = " << cut_value[0] << ", " << para_label << " = " << YLIST[f_indx] << "+/-" << YLIST_ERR[f_indx] << ", para_indx = " << para_indx << "\n";

	  // normial value
	  
	  f_indx ++;
	  
	}
      }
    }
  }

  // Determine NORM
  norm_indx = f_indx / 2;
  step1_indx = norm_indx - 1;
  step2_indx = norm_indx + 1;
  
  getnorm(Y_NORM, Y_ERR_NORM, YLIST, YLIST_ERR, norm_indx, f_indx);
  getnorm(X_NORM, X_ERR_NORM, XLIST, YLIST_ERR, norm_indx, f_indx);
  
  cout << para_label << "\n"
       << "f_indx = " << f_indx << ", norm_indx = " << norm_indx << "\n"
       << "nominal [X0, Y0] = [" << X_NORM[0] << ", " << Y_NORM[0] << "+/-" << YLIST_ERR[0] << "]\n" ; 

  double plus_err = 0., nega_err = 0.;
  
  // Type I syst. error. Weighted average +/- weighted standard deviation
  double waverage = get_waverage(YLIST, YLIST_ERR, f_indx);
  double wy_dv = get_wy_sdv(YLIST, YLIST_ERR, f_indx, waverage);
  plus_err = wy_dv;
  nega_err = wy_dv * -1;
  cout << "Type I syst. err: "<< "w_average = " << waverage << "+/-" << wy_dv << "\n"
       << "- err. = " << nega_err << ", " << "+ err = " << plus_err << "\n";
  
  //int length = sizeof(YLIST) / sizeof(YLIST[0]);
  //double syst_err = syst_typeI(YLIST, length);

  // Type II syst. error.
  fill_uncorr(UNCORR_ERR, YLIST, YLIST_ERR, f_indx);
  cout << "norm - 1: " << step1_diff << ", norm + 1: " << step2_diff << endl;

  //double syst_err = getmaxoftwo(TMath::Abs(step1_diff), TMath::Abs(step2_diff));

  //cout << "TypeII Syst. Err." << endl;
  
  //if (syst_err == TMath::Abs(step1_diff)) {
  //cout << "-1 sigma determines the SYST." << endl;
  //cout << "Syst. Err. = " << syst_err << ", z = " << step1_Z << endl;
  //}
  //else if (syst_err == TMath::Abs(step2_diff)) {
  //cout << "+1 sigma determines the SYST." << endl;
  //cout << "Syst. Err. = " << syst_err << ", z = " << step2_Z << endl;
  //}

  // Create graphs
  // gf_syst: variation as function of cut variable, uncorr. errors
  gf_syst = get_graph_syst(XLIST, YLIST, XLIST_ERR, UNCORR_ERR, f_indx);
  double upper_band_size = 3.;
  double down_band_size = 3.;
  gf_syst -> GetYaxis() -> SetRangeUser(Y_NORM[0] - down_band_size * Y_ERR_NORM[0], Y_NORM[0] + upper_band_size * Y_ERR_NORM[0]);
  gf_syst -> GetXaxis() -> SetTitle(cut_title);
  gf_syst -> GetYaxis() -> SetTitle(para_title + " " + para_unit);
  gf_syst -> SetTitle("Variation");
  gf_syst -> SetName("gf_syst");
 
  // gf_band: norm
  gf_band = get_graph_band(XLIST, BAND, X_ERR_NORM, BAND_ERR, f_indx);
  gf_norm = get_graph_norm(X_NORM, Y_NORM, X_ERR_NORM, Y_ERR_NORM, 1); //new TGraphErrors(1, );
  
  gf_syst -> RemovePoint(norm_indx);
  TMultiGraph *gf_mg = new TMultiGraph();
  gf_mg -> SetTitle("Exclusion graphs");
  gf_mg -> SetName("gf_mg");  
  gf_mg -> Add(gf_band);

  // 2sigma band
  gf_sigmaband_plus = get_graph_band(XLIST, SIGMABAND, X_ERR_NORM, SIGMABAND_ERR, f_indx);
  gf_sigmaband_nega = get_graph_band(XLIST, SIGMABAND, X_ERR_NORM, SIGMABAND_ERR, f_indx);
  
  TMultiGraph *gf_mg1 = new TMultiGraph();
  gf_mg1 -> SetTitle("");
  gf_mg1 -> SetName("gf_mg1");  
  
  // gf_Z
  gf_Z = get_graph_syst(XLIST, ZLIST, XLIST_ERR, ZLIST_ERR, f_indx);
  //gf_Z = new TGraphErrors(f_indx, XLIST, ZLIST, XLIST_ERR, ZLIST_ERR);
  gf_Z -> SetFillColor(38);
  gf_Z -> GetXaxis() -> SetTitle(cut_title);
  gf_Z -> GetYaxis() -> SetTitle("Significance Z [#sigma_{#delta}]");
  //gf_Z -> GetYaxis() -> SetTitle(para_title + " " + para_unit);
  gf_Z -> SetTitle("");
  gf_Z -> SetName("gf_Z");

  gf_Z_norm = get_graph_syst(X_NORM, Z_NORM, 0, Z_ERR_NORM, 1);
  gf_Z_norm -> SetName("gf_Z_norm");

  // plot
  TCanvas *cv_syst = new TCanvas("cv_syst", "Syst. " + cut_label, 0, 0, 1500, 800);
  cv_syst -> SetLeftMargin(0.16);
  cv_syst -> SetBottomMargin(0.15);

  gf_syst -> Draw("AP");
  gf_norm -> Draw("PZ");
  gf_mg -> Draw("C3");

  TCanvas * cv_Z = new TCanvas("cv_Z", "Significance Z " + cut_label, 0, 0, 1500, 800);
  cv_Z -> SetLeftMargin(0.16);
  cv_Z -> SetBottomMargin(0.15);

  TPaveText *pt = new TPaveText(0.4, 0.82, 0.5, 0.86, "NDC");
  
  pt -> SetTextSize(0.06);
  pt -> SetFillColor(0);
  pt -> SetTextAlign(12);
  pt -> SetTextColor(kBlack);

  pt -> AddText(para_title);

  gf_Z -> SetMarkerColorAlpha(kWhite, 0);
  gf_Z -> Draw("AP");
  gf_Z_norm -> Draw("P");
  
  //TH1F *hgf_Z = new TH1F("hgf_Z", "Converted Bar Chart", gf_Z -> GetN(), 0, gf_Z -> GetN());

  const int list_size = gf_Z -> GetN();
  const double bar_width = 0.1 * (gf_Z -> GetX()[1] - gf_Z -> GetX()[0]);
  double XC[list_size]; // X centers
  double YH[list_size]; // Y height
  cout << "list_size = " << list_size << ", bar_width = " << bar_width << endl;
  
  for (int i = 0; i < gf_Z -> GetN(); i++) {
    XC[i] = gf_Z -> GetX()[i];
    YH[i] = gf_Z -> GetY()[i];
    
    //cout << XC[i] << ", " << YH[i] << endl;
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
  
  gf_Z -> GetYaxis() -> SetRangeUser(-band_limit, band_limit);
  //cout << band_limit << endl;
  //gf_Z -> SetBarWidth(0.5); // Adjust width (0.5 = 50% of bin size)

  //gf_Z -> Draw("AP AXIS SAME"); // Draw axes without points
  //gf_Z -> Draw("AB");
  //hgf_Z -> Draw("BAR");
  gPad -> SetGrid(); // Add grid lines
  
  //gf_sigmaband_plus -> SetLineWidth(1500);
  //gf_sigmaband_plus -> SetFillStyle(3005);
  //gf_sigmaband_nega -> SetLineWidth(-1500);
  //gf_sigmaband_nega -> SetFillStyle(3005);
  gf_mg1 -> Add(gf_sigmaband_plus);
  gf_mg1 -> Add(gf_sigmaband_nega);
  gf_mg1 -> Draw("C3");
  pt -> Draw("Same");

  // Final syst. err

  TFile *f_output = new TFile(outputPlot + para_label + "_" + cut_label + ".root", "recreate");


  TTree* TRESULT = new TTree("TRESULT", "recreate");
  TRESULT ->SetAutoSave(0);
  
  TRESULT -> Branch("Br_err_type", &err_type, "Br_err_type/I");
  TRESULT -> Branch("Br_SYST_ERR", &SYST_ERR, "Br_SYST_ERR[2]/D");

  if (err_type == 1) {
    SYST_ERR[0] = nega_err;
    SYST_ERR[1] = plus_err;
    cout << "Type " << err_type << ", syst. err. (" << SYST_ERR[0] << ", " << SYST_ERR[1] << ")" << endl;
  }
  else if (err_type == 2) {
    get_syst_errII(step1_diff, step2_diff);

    cout << "Type " << err_type << ", syst. err. (" << SYST_ERR[0] << ", " << SYST_ERR[1] << ")" << endl;
  }
  else{
    cout << "new err type? !!!" << endl;
  };

  TRESULT -> Fill();
  
  
  // save
  cout << outputPlot << endl;
  cv_syst -> SaveAs(outputPlot + para_label + "_" + cut_label + ".pdf");
  cv_Z -> SaveAs(outputPlot + para_label + "_" + cut_label + "_signifi.pdf");

  gf_Z -> Write();
  gf_mg1 -> Write();
  gf_Z_norm -> Write();

  TRESULT -> Write();
  f_output -> Close();
  
  return 0;
  
}
