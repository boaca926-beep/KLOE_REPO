#include "../header/bkg_compr.h"
#include "../header/plot.h"
//#include "../hist.h"
//#include "../sf_2g/sfw1d_ref.h"

/*
TRandom *rnd=0;

Double_t DetectorEvent(Double_t mTrue) {
  // smear by double-gaussian
  if(rnd->Rndm()>frac) {
    return rnd->Gaus(mTrue+smallBias,smallSigma);
  } else {
    return rnd->Gaus(mTrue+wideBias,wideSigma);
  }
}
*/

int bkg_compr(){

  //gROOT->SetBatch(kTRUE);  
  gErrorIgnoreLevel = kError;
  TGaxis::SetMaxDigits(3);
  gStyle->SetOptStat(1110);
  gStyle->SetOptTitle(0);

  // switch on histogram errors
  TH1::SetDefaultSumw2();

  // random generator
  //rnd=new TRandom3();

  // Inspect the input tree
  //cout << infile_nm << endl;

  // get generated signal
  //const TString input_folder = "/home/bo/Desktop/input_chain_TDATA";
  //const TString input_folder = "/home/bo/Desktop/input_etagam_TDATA";
  const TString input_folder = "/home/bo/Desktop/input_vertex_TDATA";
  //const TString input_folder = "/home/bo/Desktop/input_vertex_TDATA_bkgrej_more";
  cout << input_folder << endl;
  
  TFile* intree_gen = new TFile(input_folder + "/gen/tree_gen.root");
  //TFile* intree_gen = new TFile("/home/bo/Desktop/input_chain_TUFO/gen/tree_gen.root");

  TIter next_tree1(intree_gen -> GetListOfKeys());

  TString objnm_tree, classnm_tree;

  int i = 0;
  TKey *key;
  
  while ( (key = (TKey *) next_tree1() ) ) {
    
    i ++;
    
    objnm_tree   =  key -> GetName();
    classnm_tree = key -> GetClassName();
    //key -> GetSeekKey();
    
    cout << "tree" << i << ": classnm = " << classnm_tree << ", objnm = " << objnm_tree << endl;
    
  }

  /*
  double Eisr_gen = 0., angle_isr_gen = 0.;

  TH1D* hEisr_gen = new TH1D("hEisr_gen", "", 200, 0., 500.);
  TH1D* hangle_isr_gen = new TH1D("hangle_isr_gen", "", 400, -1., 1.);

  for (Int_t irow = 0; irow < ALLCHAIN_GEN -> GetEntries(); irow++) {
      
      ALLCHAIN_GEN -> GetEntry(irow);

      if (irow > 1e7) break;
      Eisr_gen = ALLCHAIN_GEN -> GetLeaf("Br_E_pho_isr") -> GetValue(0);
      angle_isr_gen = ALLCHAIN_GEN -> GetLeaf("Br_Angle_pho_isr") -> GetValue(0);

      
      hEisr_gen -> Fill(Eisr_gen);
      hangle_isr_gen -> Fill(angle_isr_gen);
      
      //cout << Eisr_gen << endl;
      //cout << angle_isr_gen << endl;
      
      cout << irow << endl;
  
  }
  */

  TTree* ALLCHAIN_GEN = (TTree*)intree_gen -> Get("TISR3PI_SIG_GEN");
	  
  const double evnb_sig_gen = ALLCHAIN_GEN -> GetEntries(); //number of generated signal events

  // get MC recon.
  //TFile* intree = new TFile("./tree_final.root");
  //TFile* intree = new TFile("/home/bo/Desktop/analysis/crx3pi/output_angle_cut/tree_cut0.root");
  //TFile* intree = new TFile("/media/bo/8E97-E8DD/KLOE_OUTPUT/output_norm/tree_cut0.root");
     
  //TFile* intree = new TFile("/home/bo/Desktop/input_vertex_TDATA/cut/tree_pre.root");
  //TFile* intree = new TFile("/home/bo/Desktop/input_chain_TUFO/cut/tree_pre.root");
  TFile* intree = new TFile(input_folder + "/cut/tree_pre.root");
    
  //TFile* intree = new TFile("/home/bo/Desktop/analysis/crx3pi/output_pre/tree_cut0.root");

  TIter next_tree(intree -> GetListOfKeys());

  //TString objnm_tree, classnm_tree;

  i = 0;
  
  while ( (key = (TKey *) next_tree() ) ) {
    
    i ++;
    
    objnm_tree   =  key -> GetName();
    classnm_tree = key -> GetClassName();
    //key -> GetSeekKey();
    
    cout << "tree" << i << ": classnm = " << classnm_tree << ", objnm = " << objnm_tree << endl;
    
  }

  // check branches
  TTree* TDATA = (TTree*)intree -> Get("TDATA");
  TTree* TEEG = (TTree*)intree -> Get("TEEG");
  TTree* TOMEGAPI = (TTree*)intree -> Get("TOMEGAPI");
  TTree* TKSL = (TTree*)intree -> Get("TKSL");
  TTree* TKPM = (TTree*)intree -> Get("TKPM");
  TTree* TRHOPI = (TTree*)intree -> Get("TRHOPI");
  TTree* TETAGAM = (TTree*)intree -> Get("TETAGAM");
  TTree* TBKGREST = (TTree*)intree -> Get("TBKGREST");
  TTree* TISR3PI_SIG = (TTree*)intree -> Get("TISR3PI_SIG");

  //TDATA -> GetListOfLeaves() -> Print();

  
  //bkg
  const int TLSize = 9;
  char name[TLSize], title[TLSize];

  TTree *TrList[TLSize] = {TDATA, TEEG, TOMEGAPI, TKSL, TKPM, TRHOPI, TETAGAM, TBKGREST, TISR3PI_SIG};
  const TString TrNm[TLSize] = {"data", "eeg", "omegapi" , "ksl", "kpm", "rhopi" , "etagam" , "bkgrest", "isr3pi"};

  int color_list[TLSize] = {1, 6, 7, 28, 46, 42, 3, 37, 4};
  
  // Create arrays of histograms.
  TObjArray Hlist(100);
  TObjArray H2dlist(100);

  TH1D* h;
  TH2D* h2d;
  TH2D* h2d_pchi2;

  // fixed variables
  double var_value = 0., IM3pi = 0., IM3pi_true = 0., IM3pi_det = 0.;
  double evnb_data = 0., evnb_eeg = 0., evnb_omegapi = 0., evnb_ksl = 0., evnb_kpm = 0., evnb_rhopi = 0., evnb_etagam = 0., evnb_bkgrest = 0., evnb_isr3pi = 0.;
  double betapi0 = 0., IM2pi = 0.;
  double chi2 = 0., pvalue = 0.;
  
  //if (var_nm.Contains("IM3pi_7C")) cout << var_nm << endl;
  
  for (int i = 0; i < TLSize; i ++) {// start MC type loop

    //sprintf(name,"h%d",i);
    //sprintf(title,"histo nr:%d",i);
    h = new TH1D("hist_" + TrNm[i], "", binsize, var_min, var_max);
    h2d = new TH2D("h2d_discrp_" + TrNm[i], "", 200, 0.25, 0.65, 200, 0.3, 1.);
    h2d_pchi2 = new TH2D("h2d_pchi2_" + TrNm[i], "", 200, 0., 1., 200, 0., 100.);

    cout << TrNm[i] << ", Tree: " << TrList[i] -> GetName() << endl;

    for (Int_t irow = 0; irow < TrList[i] -> GetEntries(); irow++) {
      
      TrList[i] -> GetEntry(irow);
      
      var_value = TrList[i] -> GetLeaf("Br_" + var_nm) -> GetValue(0);
      
      if (TrNm[i] == "data") evnb_data ++;
      else if (TrNm[i] == "eeg") evnb_eeg ++;
      else if (TrNm[i] == "omegapi") evnb_omegapi ++;
      else if (TrNm[i] == "ksl") evnb_ksl ++;
      else if (TrNm[i] == "kpm") evnb_kpm ++;
      else if (TrNm[i] == "rhopi") evnb_rhopi ++;
      else if (TrNm[i] == "etagam") evnb_etagam ++;
      else if (TrNm[i] == "bkgrest") evnb_bkgrest ++;
      else if (TrNm[i] == "isr3pi") evnb_isr3pi ++;
      
      IM3pi_true = TrList[i] -> GetLeaf("Br_IM3pi_true") -> GetValue(0);
      IM3pi = TrList[i] -> GetLeaf("Br_IM3pi_7C") -> GetValue(0);
      //IM3pi_det = DetectorEvent(TMath::Abs(IM3pi_true));
      
      betapi0 = TrList[i] -> GetLeaf("Br_betapi0") -> GetValue(0);
      IM2pi = TrList[i] -> GetLeaf("Br_ppIM") -> GetValue(0) * 1e-3;
      
      chi2 = TrList[i] -> GetLeaf("Br_lagvalue_min_7C") -> GetValue(0);
      //pvalue = TrList[i] -> GetLeaf("Br_pvalue") -> GetValue(0);

      h2d -> Fill(IM2pi, betapi0);
      h2d_pchi2 -> Fill(pvalue, chi2);
      //cout << "betapi0 = " << betapi0 << ", IM2pi = " << IM2pi << " [GeV/c^2]" << endl;
      //if (IM3pi < IM3pi_max && IM3pi > IM3pi_min) {
	
      //cout << var_nm << " = " <<  var_value << endl;
      
      //if (var_nm.Contains("IM3pi_7C") && TrNm[i].Contains("isr3pi")) {
	
	//var_value = IM3pi_det;
	//cout << "IM3pi_true = " << IM3pi_true << ", IM3pi_det = " << IM3pi_det << "\n";
	  
      //}

      //cout << var_value << endl;

      //if (IM3pi < IM3pi_max && IM3pi > IM3pi_min) {

	h -> Fill(var_value);
	
	//}
	
    }
    
    format_h(h, color_list[i], 2);
    Hlist.Add(h);

    //cout << h -> GetName() << endl;
    
    H2dlist.Add(h2d);
    H2dlist.Add(h2d_pchi2);
    
  }

  
  //double evnb_data = 0., evnb_eeg = 0., evnb_omegapi = 0., evnb_ksl = 0., evnb_kpm = 0., evnb_rhopi = 0., evnb_etagam = 0., evnb_bkgrest = 0., evnb_isr3pi = 0.;

  //double evnb_isr3pi = evnb_isr3pi;
  evnb_eeg = evnb_eeg * 2.;
  const double evnb_mcrest = evnb_kpm + evnb_rhopi + evnb_bkgrest;
  const double evnb_mcsum = evnb_eeg + evnb_omegapi + evnb_ksl + evnb_etagam + evnb_mcrest + evnb_isr3pi;
  const double evnb_bkgsum = evnb_eeg + evnb_omegapi + evnb_ksl + evnb_etagam + evnb_mcrest;
  const double sb_ratio = evnb_isr3pi / TMath::Sqrt(evnb_mcsum);
  const double s_frac = evnb_isr3pi / evnb_mcsum * 100.;
  const double b_frac = evnb_bkgsum / evnb_mcsum * 100.;
    
  ofstream myfile;
  myfile.open ("./output/compr.txt");
  
  myfile << "Number of events\n"
	 << "data = " << evnb_data << "\n"
	 << "1. eeg = " << evnb_eeg << "\n"
	 << "2. omegapi = " << evnb_omegapi << "\n"
	 << "3. ksl = " << evnb_ksl << "\n"
	 << "4. etagam = " << evnb_etagam << "\n"
	 << "5. mcrest = " << evnb_mcrest << "\n"
	 << "\tbkgrest = " << evnb_bkgrest << "\n"
	 << "\tkpm = " << evnb_kpm << "\n"
	 << "\trhopi = " << evnb_rhopi << "\n"
	 << "6. isr3pi (gen) = " << evnb_isr3pi << " (" << evnb_sig_gen << ")\n"
	 << "mcsum = " << evnb_mcsum << "\n"
	 << "s_frac [%] = " << s_frac << "\n"
	 << "b_frac [%] = " << b_frac << "\n"  
	 << "sb_ratio = " << sb_ratio << "\n";
  
  // MC rest, merge bkgrest, kpm and rhopi
  TH1D *hist_bkgrest = (TH1D *) Hlist.FindObject("hist_bkgrest");  
  TH1D *hist_kpm = (TH1D *) Hlist.FindObject("hist_kpm");  
  TH1D *hist_rhopi = (TH1D *) Hlist.FindObject("hist_rhopi");  
  TH1D *hist_eeg = (TH1D *) Hlist.FindObject("hist_eeg");  
  TH1D *hist_isr3pi = (TH1D *) Hlist.FindObject("hist_isr3pi");  
  TH1D *hist_omegapi = (TH1D *) Hlist.FindObject("hist_omegapi");  
  TH1D *hist_etagam = (TH1D *) Hlist.FindObject("hist_etagam");  
  TH1D *hist_ksl = (TH1D *) Hlist.FindObject("hist_ksl");  
 
  TH1D *hist_mcrest = (TH1D*) hist_bkgrest -> Clone();
  hist_mcrest -> Add(hist_kpm, 1.);
  hist_mcrest -> Add(hist_rhopi, 1.);
  hist_mcrest -> SetName("hist_mcrest");

  Hlist.Add(hist_mcrest);

  // add 2d
  //TH2D* h2d_discrp_mcrest = (TH2D*) h2d_discrp_bkgrest -> Clone();
  //h2d_discrp_mcrest -> Add(h2d_discrp_kpm, 1.);
  //h2d_discrp_mcrest -> Add(h2d_discrp_rhopi, 1.);
  //h2d_discrp_mcrest -> SetName("h2d_discrp_mcrest");

  TH2D *h2d_discrp_rhopi = (TH2D *) H2dlist.FindObject("h2d_discrp_rhopi");  
  TH2D *h2d_discrp_eeg = (TH2D *) H2dlist.FindObject("h2d_discrp_eeg");  
  TH2D *h2d_discrp_omegapi = (TH2D *) H2dlist.FindObject("h2d_discrp_omegapi");  
  TH2D *h2d_discrp_ksl = (TH2D *) H2dlist.FindObject("h2d_discrp_ksl");  
  TH2D *h2d_discrp_isr3pi = (TH2D *) H2dlist.FindObject("h2d_discrp_isr3pi");  
  TH2D *h2d_discrp_etagam = (TH2D *) H2dlist.FindObject("h2d_discrp_etagam");  

  TH2D *h2d_discrp_mcrest = (TH2D*) h2d_discrp_rhopi -> Clone();
  
  H2dlist.Add(h2d_discrp_mcrest);
  // eeg
  TH1D * hist_eeg_sc = (TH1D*) hist_eeg -> Clone();
  hist_eeg_sc -> Scale(2.);
  hist_eeg_sc -> SetName("hist_eeg_sc");

  TH2D * h2d_discrp_eeg_sc = (TH2D*) h2d_discrp_eeg -> Clone();
  h2d_discrp_eeg_sc -> Scale(2.);
  h2d_discrp_eeg_sc -> SetName("h2d_discrp_eeg_sc");

  // isr3pi
  TH1D * hist_isr3pi_sc = (TH1D*) hist_isr3pi -> Clone();
  hist_isr3pi_sc -> Scale(1); //sfw1d_isr3pi
  hist_isr3pi_sc -> SetName("hist_isr3pi_sc");

  TH2D * h2d_discrp_isr3pi_sc = (TH2D*) h2d_discrp_isr3pi -> Clone();
  //h2d_discrp_isr3pi_sc -> Scale(sfw1d_isr3pi);
  h2d_discrp_isr3pi_sc -> SetName("h2d_discrp_isr3pi_sc");

  //cout << "sfw1d_isr3pi = " << sfw1d_isr3pi << endl;
  
  // omegapi
  TH1D * hist_omegapi_sc = (TH1D*) hist_omegapi -> Clone();
  hist_omegapi_sc -> Scale(1.);
  hist_omegapi_sc -> SetName("hist_omegapi_sc");

  TH2D * h2d_discrp_omegapi_sc = (TH2D*) h2d_discrp_omegapi -> Clone();
  h2d_discrp_omegapi_sc -> Scale(1.);
  h2d_discrp_omegapi_sc -> SetName("h2d_discrp_omegapi_sc");

  // etagam
  TH1D * hist_etagam_sc = (TH1D*) hist_etagam -> Clone();
  hist_etagam_sc -> Scale(1.);
  hist_etagam_sc -> SetName("hist_etagam_sc");

  TH2D * h2d_discrp_etagam_sc = (TH2D*) h2d_discrp_etagam -> Clone();
  h2d_discrp_etagam_sc -> Scale(1.);
  h2d_discrp_etagam_sc -> SetName("h2d_discrp_etagam_sc");

  // ksl
  TH1D * hist_ksl_sc = (TH1D*) hist_ksl -> Clone();
  hist_ksl_sc -> Scale(1.);
  hist_ksl_sc -> SetName("hist_ksl_sc");

  TH2D * h2d_discrp_ksl_sc = (TH2D*) h2d_discrp_ksl -> Clone();
  h2d_discrp_ksl_sc -> Scale(1.);
  h2d_discrp_ksl_sc -> SetName("h2d_discrp_ksl_sc");

  // mcrest
  TH1D * hist_mcrest_sc = (TH1D*) hist_mcrest -> Clone();
  hist_mcrest_sc -> Scale(1.);
  hist_mcrest_sc -> SetName("hist_mcrest_sc");

  TH2D * h2d_discrp_mcrest_sc = (TH2D*) h2d_discrp_mcrest -> Clone();
  h2d_discrp_mcrest_sc -> Scale(1.);
  h2d_discrp_mcrest_sc -> SetName("h2d_discrp_mcrest_sc");

  // bkgsum
  TH1D* hist_bkgsum_sc = (TH1D*) hist_eeg_sc -> Clone();
  hist_bkgsum_sc -> Add(hist_omegapi_sc, 1.);
  hist_bkgsum_sc -> Add(hist_ksl_sc, 1.);
  hist_bkgsum_sc -> Add(hist_etagam_sc, 1.);
  hist_bkgsum_sc -> Add(hist_mcrest_sc, 1.);
  hist_bkgsum_sc -> SetName("hist_bkgsum_sc");

  format_h(hist_bkgsum_sc, 6, 2);

  TH2D* h2d_discrp_bkgsum_sc = (TH2D*) h2d_discrp_eeg_sc -> Clone();
  h2d_discrp_bkgsum_sc -> Add(h2d_discrp_omegapi_sc, 1.);
  h2d_discrp_bkgsum_sc -> Add(h2d_discrp_ksl_sc, 1.);
  h2d_discrp_bkgsum_sc -> Add(h2d_discrp_etagam_sc, 1.);
  h2d_discrp_bkgsum_sc -> Add(h2d_discrp_mcrest_sc, 1.);
  h2d_discrp_bkgsum_sc -> SetName("h2d_discrp_bkgsum_sc");

  // nomalization and add to the histo list
  
  Hlist.Add(hist_eeg_sc);
  Hlist.Add(hist_isr3pi_sc);
  Hlist.Add(hist_omegapi_sc);
  Hlist.Add(hist_etagam_sc);
  Hlist.Add(hist_ksl_sc);
  Hlist.Add(hist_mcrest_sc);
  Hlist.Add(hist_bkgsum_sc);

  H2dlist.Add(h2d_discrp_eeg_sc);
  H2dlist.Add(h2d_discrp_isr3pi_sc);
  H2dlist.Add(h2d_discrp_omegapi_sc);
  H2dlist.Add(h2d_discrp_etagam_sc);
  H2dlist.Add(h2d_discrp_ksl_sc);
  H2dlist.Add(h2d_discrp_mcrest_sc);
  H2dlist.Add(h2d_discrp_bkgsum_sc);
  
  
  // save
  TFile *f_out = new TFile(output_folder + "/hist_" + var_nm + ".root", "recreate");

  Hlist.Write();
  H2dlist.Write();
  //hEisr_gen -> Write();
  //hangle_isr_gen -> Write();
    
  f_out -> Close();
  
  return 0;
  
}


