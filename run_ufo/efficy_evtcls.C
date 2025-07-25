#include "../header/efficy.h"
#include "../header/sm_para.h"
#include "../header/method.h"
#include "../header/graph.h"

double binomial_err(double nb_true, double nb_gen) {
  double error = 0.;
  double ratio = 0.; 

  if (nb_gen != 0.) {
    ratio = nb_true / nb_gen;
    error = TMath::Sqrt(ratio * (1. - ratio) / nb_gen);
  }
   
  //cout << "true = " << nb_true << ", gen = " << nb_gen << ", ratio = " << ratio << ", error = " << error << endl;

  return error;
}

TRandom *generator = new TRandom();
rnd=new TRandom3();

int efficy_evtcls() {

  cout << "Calculate efficiency ..." << endl; 

  cout << "input path: " << inputFile << ", tree type = " << treeType << endl;

  TFile *f_input = new TFile(inputFile + ".root");
  //getObj(f_input);

  TTree *ALLCHAIN_CUT = (TTree*)f_input -> Get(treeType);

  int evtcls_indx = -1, trigger_indx = -1, filfo_indx = -1, sel_indx = -1;
  int filfo28_indx = -1;
  
  double evnt_trigger = 0.; // number of events after the trigger
  double evnt_filfo = 0., N0 = 0.; // number of events after the filfo
  double evnt_bkg = 0.; // number of backgroud events
  double evnt_sel = 0.; // number of events after selection
  double evnt_evtcls = 0.; // number of events after event classification
  
  double m3pi = 0., m3pi_true = 0.;
  double efficy_cut = 0.;
  
  TFile *f_output = new TFile(outputCut + "efficy.root", "update");
  
  TTree* TRESULT = new TTree(treeType, "recreate");
  TRESULT -> SetAutoSave(0);

  TRESULT -> Branch("Br_evnt_sel", &evnt_sel, "Br_evnt_sel/D");
  TRESULT -> Branch("Br_evnt_evtcls", &evnt_evtcls, "Br_evnt_evtcls/D");
  TRESULT -> Branch("Br_efficy_cut", &efficy_cut, "Br_efficy_cut/D");

  // initialize histos
  for (int i = 0; i < list_size; i ++) {

    H1DLIST[i] = new TH1D("h1d_IM3pi_" + CUT_TYPE[i] + "_" + treeType, "", IM3pi_bin, IM3pi_min, IM3pi_max);
    H1DLIST[i] -> Sumw2();

    //cout << CUT_TYPE[i] << endl; 
  
  }

  //ALLCHAIN_CUT -> Print();
    
  for (Int_t irow = 0; irow < ALLCHAIN_CUT -> GetEntries(); irow ++) {// loop trees
	  
    ALLCHAIN_CUT -> GetEntry(irow);

    if (treeType == "TISR3PI_SIG") {//
      //m3pi = ALLCHAIN_CUT -> GetLeaf("Br_IM3pi_7C") -> GetValue(0);
      m3pi_true = ALLCHAIN_CUT -> GetLeaf("Br_IM3pi_true") -> GetValue(0);
      m3pi = DetectorEvent(TMath::Abs(m3pi_true));
    }
    else {
      m3pi = ALLCHAIN_CUT -> GetLeaf("Br_IM3pi_7C") -> GetValue(0);
    }
    
    trigger_indx = ALLCHAIN_CUT -> GetLeaf("Br_trigger_indx") -> GetValue(0);
    filfo_indx = ALLCHAIN_CUT -> GetLeaf("Br_filfo_indx") -> GetValue(0);
    filfo28_indx = ALLCHAIN_CUT -> GetLeaf("Br_filfo28_indx") -> GetValue(0);
    evtcls_indx = ALLCHAIN_CUT -> GetLeaf("Br_evtcls_indx") -> GetValue(0);
    sel_indx = ALLCHAIN_CUT -> GetLeaf("Br_sel_indx") -> GetValue(0);

    //cout << sel_indx << endl;

    //if (m3pi > 760 && m3pi < 800) {
      evnt_sel ++;
    
      //if (trigger_indx == 0) continue; // trigger
      evnt_trigger ++;

      if (filfo28_indx == 0) {// filfo28
	N0 ++;
	continue; 
      }
      evnt_filfo ++;
      
      //cout << trigger_indx << " " << filfo_indx << endl;
      //cout << filfo28_indx << endl;
      
      if (sel_indx == 0) continue; // background rejection
      evnt_bkg ++;
    
      //if (evtcls_indx == 0) continue; // evnt classification
      //evnt_evtcls ++;
    
      H1DLIST[0] -> Fill(m3pi);
    
      if (evtcls_indx == 0) continue; // evnt classification 
      evnt_evtcls ++;
      
      H1DLIST[1] -> Fill(m3pi);
    
      //if (evtcls_indx == 0) cout << evtcls_indx << endl;
      //if (filfo_indx == 0) cout << filfo_indx << endl;
      //if (trigger_indx == 0) cout << trigger_indx << endl;
    
      //cout << m3pi << endl;

      //}
    
  }

  //h1d_IM3pi -> Draw();

  efficy_cut = evnt_evtcls / evnt_bkg;
  //efficy_cut = evnt_bkg / evnt_evtcls;

  cout << "evnt_sel = " << evnt_sel << "\n"
       << "evnt_trigger = " << evnt_trigger << "\n"
       << "evnt_filfo (N28) = " << evnt_filfo << ", N28 = " << evnt_filfo << ", N0+N28 = " << N0 + evnt_filfo << ", efficy_filfo = " << N0 / evnt_filfo / 9 <<  "\n"
       << "evnt_bkg = " << evnt_bkg << "\n"
       << "evnt_evtcls = " << evnt_evtcls << "\n"
       << "efficy_cut = " << efficy_cut << endl;

  TRESULT -> Fill();

  HIM3PI -> Add(H1DLIST[0]);
  HIM3PI -> Add(H1DLIST[1]);

  // efficiency
  double efficy_tmp = 0., efficy_err_tmp = 0.;
  double nb_sel = 0., nb_sel_err = 0.;
  double nb_evtcls = 0., nb_evtcls_err = 0.;

  const int binsize = H1DLIST[0] -> GetNbinsX();
  double xmin = H1DLIST[0] -> GetXaxis() -> GetXmin();
  double xmax = H1DLIST[0] -> GetXaxis() -> GetXmax();

  double M3PI[binsize], M3PI_ERR[binsize];
  double EFFICY[binsize], EFFICY_ERR[binsize];
  double NB_SEL[binsize], NB_SEL_ERR[binsize];
  double NB_EVTCLS[binsize], NB_EVTCLS_ERR[binsize];
  
  for (int i = 1; i <= binsize; i ++) {

    nb_sel = H1DLIST[0] -> GetBinContent(i);
    nb_sel_err = H1DLIST[0] -> GetBinError(i);
    
    nb_evtcls = H1DLIST[1] -> GetBinContent(i);
    nb_evtcls_err = H1DLIST[1] -> GetBinError(i);

    //cout << "nb_sel = " << nb_sel << "+/-" << nb_sel_err << endl;
    //cout << "nb_evtcls = " << nb_evtcls << "+/-" << nb_evtcls_err << endl;

    NB_SEL[i - 1] = nb_sel;
    NB_SEL_ERR[i - 1] = nb_sel_err;

    NB_EVTCLS[i - 1] = nb_evtcls;
    NB_EVTCLS_ERR[i - 1] = nb_evtcls_err;
    
  
    M3PI[i - 1] = H1DLIST[0] -> GetBinCenter(i);
    M3PI_ERR[i - 1] = 0.;

    if (nb_sel == 0. || nb_evtcls == 0.) {

      EFFICY[i - 1] = 0.; 
      EFFICY_ERR[i - 1] = 0.;
      
    }
    else {

      EFFICY[i - 1] = nb_evtcls / nb_sel; 
      EFFICY_ERR[i - 1] = binomial_err(nb_evtcls, nb_sel);
    
    }

    //cout << "mass bin " << i << ", mass = " << M3PI[i - 1] << ", nb_evtcls = " << nb_evtcls << ", nb_sel = " << nb_sel << ", efficy_tmp = " << EFFICY[i - 1] << "+/-" << EFFICY_ERR[i - 1] << endl; 
    
  }

  TGraphErrors *gf_efficy = get_graph_syst(M3PI, EFFICY, M3PI_ERR, EFFICY_ERR, binsize);
  gf_efficy -> SetName("gf_efficy_" + treeType);

  TGraphErrors *gf_nb_sel = get_graph_syst(M3PI, NB_SEL, M3PI_ERR, NB_SEL_ERR, binsize);
  gf_nb_sel -> SetName("gf_nb_sel_" + treeType);

  TGraphErrors *gf_nb_evtcls = get_graph_syst(M3PI, NB_EVTCLS, M3PI_ERR, NB_EVTCLS_ERR, binsize);
  gf_nb_evtcls -> SetName("gf_nb_evtcls_" + treeType);

  /*
  // plot
  TCanvas *cv_efficy = new TCanvas("cv_efficy", "cv_efficy", 1000, 800);
  cv_efficy -> SetLeftMargin(0.1);

  gf_efficy -> GetYaxis() -> SetNdivisions(512);
  gf_efficy -> GetYaxis() -> SetRangeUser(0., 1.);
  gf_efficy -> GetYaxis() -> SetTitle("#varepsilon");
  gf_efficy -> GetYaxis() -> SetTitleSize(0.04);
  gf_efficy -> GetYaxis() -> SetTitleOffset(1.2);
  gf_efficy -> GetYaxis() -> SetLabelSize(0.035);
  gf_efficy -> GetYaxis() -> CenterTitle();
  
  gf_efficy -> GetXaxis() -> SetTitle("M_{3#pi} [MeV/c^{2}]");
  gf_efficy -> GetXaxis() -> SetTitleOffset(1.1);
  gf_efficy -> GetXaxis() -> SetTitleSize(0.04);
  gf_efficy -> GetXaxis() -> SetLabelSize(0.04);
  gf_efficy -> GetXaxis() -> CenterTitle();
  //gf_efficy -> GetXaxis() -> SetRangeUser(760., 800.);
  
  */

  //gf_efficy -> Draw("AP");
  gf_nb_sel -> Draw("AP");
  gf_nb_evtcls -> Draw("P");
  
  /// save
  TRESULT -> Write();

  HIM3PI -> Write("HIM3PI_" + treeType, 1);
  gf_efficy -> Write();
  gf_nb_sel -> Write();
  gf_nb_evtcls -> Write();
  
  f_output -> Close();

  return 0;
  
}
