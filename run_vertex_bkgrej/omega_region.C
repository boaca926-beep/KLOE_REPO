#include "../header/sm_para.h"
#include "../header/plot.h"
#include "../header/method.h"

TRandom *generator = new TRandom();
rnd=new TRandom3();

TObjArray *Hlist = new TObjArray(100);
TObjArray *Hlist_sc = new TObjArray(100);

void checkArray(TObjArray *array){

  // Create a TIter object for the TObjArray
  TIter next(array);

  TObject* object = 0;
  int obj_indx = 0;
  while ((object = next()))
    {

      cout << "[" << obj_indx << "]: " << object -> GetName() << endl;
      obj_indx ++;
      
    }
  
}

void inspect_input(TFile *f){// File inspection

  TIter next_tree1(f -> GetListOfKeys());

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
  
}

void fill_hist(TString TrNm[], TTree *TrList[], int TLSize, int color_list[]) {// Fill histos

  // Define variables
  double var_value = 0., IM3pi = 0., IM3pi_true = 0., IM3pi_det = 0.;
  double evnb_data = 0., evnb_eeg = 0., evnb_omegapi = 0., evnb_ksl = 0., evnb_kpm = 0., evnb_rhopi = 0., evnb_etagam = 0., evnb_bkgrest = 0., evnb_isr3pi = 0.;
  double betapi0 = 0., IM2pi = 0.;
  double chi2 = 0., pvalue = 0.;

  const double mass_min = 770; //380, 770
  const double mass_max = 800; //1020, 800
  const double mass_sigma_nb = 1; //1, 0.25
  const double IM3pi_sigma = 2.65; //2.615, 2.65
  //const int binsize = TMath::Nint((mass_max - mass_min) / mass_sigma_nb / IM3pi_sigma);
  const int binsize = 60;
  const TString var_nm = "IM3pi_7C";

  TH1D* h, *hsig;

  for (int i = 0; i < TLSize; i ++) {// start MC type loop

    h = new TH1D("hist_" + TrNm[i], "", binsize, mass_min, mass_max);

    //cout << TrNm[i] << ", Tree: " << TrList[i] -> GetName() << endl;

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

      h -> Fill(var_value);

      //if (irow > 1000) break;

    }

    format_h(h, color_list[i], 2);
    Hlist -> Add(h);

    // Smearing signal
    double m3pi_true = 0., m3pi_corr = 0.;
    
    if (std::string(TrList[i] -> GetName()) == "TISR3PI_SIG") {
      
      hsig = new TH1D("hsig", "", binsize, mass_min, mass_max);
      
      cout << "Starting smearing the signal ..." << endl;
      
      for (Int_t irow = 0; irow < TrList[i] -> GetEntries(); irow++) {// loop chain
	TrList[i] -> GetEntry(irow);
	
	m3pi_true = TrList[i] -> GetLeaf("Br_IM3pi_true") -> GetValue(0);
	
	m3pi_corr = DetectorEvent(TMath::Abs(m3pi_true));
	
	//cout << m3pi_true << endl;
	
	hsig -> Fill(m3pi_corr);
      }
      
      format_h(hsig, color_list[i], 2);
      Hlist -> Add(hsig);

    }
      
  }// end MC type loop

  
  // Summary
  evnb_eeg = evnb_eeg * 2.;
  const double evnb_mcrest = evnb_kpm + evnb_rhopi + evnb_bkgrest;
  const double evnb_mcsum = evnb_eeg + evnb_omegapi + evnb_ksl + evnb_etagam + evnb_mcrest + evnb_isr3pi;
  const double evnb_bkgsum = evnb_eeg + evnb_omegapi + evnb_ksl + evnb_etagam + evnb_mcrest;
  const double sb_ratio = evnb_isr3pi / TMath::Sqrt(evnb_mcsum);
  const double s_frac = evnb_isr3pi / evnb_mcsum * 100.;
  const double b_frac = evnb_bkgsum / evnb_mcsum * 100.;

  cout << "Summary of events\n"
       << "data = " << evnb_data << "\n"
       << "1. eeg = " << evnb_eeg << "\n"
       << "2. omegapi = " << evnb_omegapi << "\n"
       << "3. ksl = " << evnb_ksl << "\n"
       << "4. etagam = " << evnb_etagam << "\n"
       << "5. mcrest = " << evnb_mcrest << "\n"
       << "\tbkgrest = " << evnb_bkgrest << "\n"
       << "\tkpm = " << evnb_kpm << "\n"
       << "\trhopi = " << evnb_rhopi << "\n"
    //<< "6. isr3pi (gen) = " << evnb_isr3pi << " (" << evnb_sig_gen << ")\n"
       << "mcsum = " << evnb_mcsum << "\n"
       << "s_frac = " << s_frac << "%\n"
       << "b_frac = " << b_frac << "%\n"  
       << "sb_ratio = " << sb_ratio << "\n";

  // MC rest, merge bkgrest, kpm and rhopi
  TH1D *hist_bkgrest = (TH1D *) Hlist ->  FindObject("hist_bkgrest");  
  TH1D *hist_kpm = (TH1D *) Hlist ->  FindObject("hist_kpm");  
  TH1D *hist_rhopi = (TH1D *) Hlist ->  FindObject("hist_rhopi");  
  //TH1D *hist_eeg = (TH1D *) Hlist ->  FindObject("hist_eeg");  
  //TH1D *hist_isr3pi = (TH1D *) Hlist ->  FindObject("hist_isr3pi");  
  //TH1D *hist_omegapi = (TH1D *) Hlist ->  FindObject("hist_omegapi");  
  //TH1D *hist_etagam = (TH1D *) Hlist ->  FindObject("hist_etagam");  
  //TH1D *hist_ksl = (TH1D *) Hlist ->  FindObject("hist_ksl");  
 
  TH1D *hist_mcrest = (TH1D*) hist_bkgrest -> Clone();
  hist_mcrest -> Add(hist_kpm, 1.);
  hist_mcrest -> Add(hist_rhopi, 1.);
  hist_mcrest -> SetName("hist_mcrest");

  Hlist ->  Add(hist_mcrest);

}

double eeg_sfw = 0.; 
double isr3pi_sfw = 0.; 
double omegapi_sfw = 0.; 
double etagam_sfw = 0.; 
double ksl_sfw = 0.; 
double mcrest_sfw = 0.; 

double sig_sfw = 0.; 

  
int omega_region() {
  //// To understand more of MC and data in the omega region M3pi = [770, 800] MeV/c^{2}

  //gROOT->SetBatch(kTRUE);  
  gErrorIgnoreLevel = kError;
  TGaxis::SetMaxDigits(3);
  gStyle->SetOptStat(1110);
  gStyle->SetOptTitle(0);
  TH1::SetDefaultSumw2(); // switch on histogram errors

  const TString input_folder = "/home/bo/Desktop/analysis_root_v6/input_norm_TDATA";
  //const TString input_folder = "/home/bo/Desktop/input_norm_TDATA";
  
  //gROOT->SetBatch(kTRUE);  
  gErrorIgnoreLevel = kError;
  TGaxis::SetMaxDigits(3);
  gStyle->SetOptStat(1110);
  gStyle->SetOptTitle(0);
  gStyle->SetStatBorderSize(0);
  gStyle->SetFitFormat("6.2g");

  //// Input files
  TFile* file_gen = new TFile(input_folder + "/gen/tree_gen.root");
  TTree* ALLCHAIN_GEN = (TTree*)file_gen -> Get("TISR3PI_SIG_GEN");

  TFile* file_cut = new TFile(input_folder + "/cut/tree_pre.root");
  TFile *file_sfw2d = new TFile(input_folder + "/sfw2d/sfw2d.root");  
  TFile *file_sfw1d = new TFile(input_folder + "/sfw1d/sfw1d.root");  

  //inspect_input(file_gen);
  //inspect_input(file_cut);

  // Scaling factors

  // sfw1d
  TTree *TSFW1D = (TTree*)file_sfw1d -> Get("TRESULT");

  for (Int_t irow = 0; irow < TSFW1D -> GetEntries(); irow++) {// loop chain

    TSFW1D -> GetEntry(irow);

    sig_sfw = TSFW1D -> GetLeaf("Br_sig_sfw") -> GetValue(0);
    
  }
  
  // sfw2d 
  TTree *TSFW2D = (TTree*)file_sfw2d -> Get("TRESULT");

  for (Int_t irow = 0; irow < TSFW2D -> GetEntries(); irow++) {// loop chain

    TSFW2D -> GetEntry(irow);
    
    eeg_sfw = TSFW2D -> GetLeaf("Br_eeg_sfw") -> GetValue(0);
    isr3pi_sfw = TSFW2D -> GetLeaf("Br_isr3pi_sfw") -> GetValue(0);
    omegapi_sfw = TSFW2D -> GetLeaf("Br_omegapi_sfw") -> GetValue(0);
    etagam_sfw = TSFW2D -> GetLeaf("Br_etagam_sfw") -> GetValue(0);
    ksl_sfw = TSFW2D -> GetLeaf("Br_ksl_sfw") -> GetValue(0);
    mcrest_sfw = TSFW2D -> GetLeaf("Br_mcrest_sfw") -> GetValue(0);

  }

  cout << "Saling factors \n"
       << "eeg_sfw = " << eeg_sfw << "\n"
       << "isr3pi_sfw = " << isr3pi_sfw << "\n"
       << "omegapi_sfw = " << omegapi_sfw << "\n" 
       << "etagam_sfw = " << etagam_sfw << "\n" 
       << "ksl_sfw = " << ksl_sfw << "\n"
       << "mcrest_sfw = " << mcrest_sfw << "\n"
       << "sig_sfw = " << sig_sfw << endl;
  
  /// Get branches
  TTree* TDATA = (TTree*)file_cut -> Get("TDATA");
  TTree* TEEG = (TTree*)file_cut -> Get("TEEG");
  TTree* TOMEGAPI = (TTree*)file_cut -> Get("TOMEGAPI");
  TTree* TKSL = (TTree*)file_cut -> Get("TKSL");
  TTree* TKPM = (TTree*)file_cut -> Get("TKPM");
  TTree* TRHOPI = (TTree*)file_cut -> Get("TRHOPI");
  TTree* TETAGAM = (TTree*)file_cut -> Get("TETAGAM");
  TTree* TBKGREST = (TTree*)file_cut -> Get("TBKGREST");
  TTree* TISR3PI_SIG = (TTree*)file_cut -> Get("TISR3PI_SIG");

  //// Preparing histograms
  const int TLSize = 9;
  TString TrNm[TLSize] = {"data", "eeg", "omegapi" , "ksl", "kpm", "rhopi" , "etagam" , "bkgrest", "isr3pi"};
  TTree *TrList[TLSize] = {TDATA, TEEG, TOMEGAPI, TKSL, TKPM, TRHOPI, TETAGAM, TBKGREST, TISR3PI_SIG};
  int color_list[TLSize] = {1, 6, 7, 28, 46, 42, 3, 37, 4};
  
  fill_hist(TrNm, TrList, TLSize, color_list);
  
  //// Check histo array

  checkArray(Hlist); // Check array

  /// Histo normalization

  // Background
  TH1D *hist_eeg_sc = (TH1D *) Hlist -> FindObject("hist_eeg") -> Clone();
  hist_eeg_sc -> Scale(eeg_sfw * 2);
  hist_eeg_sc -> SetName("hist_eeg_sc");
  //cout << "nb_eeg = " << hist_eeg -> Integral(1, hist_eeg -> GetNbinsX()) << endl;
  Hlist_sc -> Add(hist_eeg_sc);
  
  TH1D *hist_omegapi_sc = (TH1D *) Hlist -> FindObject("hist_omegapi") -> Clone();
  hist_omegapi_sc -> Scale(omegapi_sfw);
  hist_omegapi_sc -> SetName("hist_omegapi_sc");
  Hlist_sc -> Add(hist_omegapi_sc);
  
  TH1D *hist_ksl_sc = (TH1D *) Hlist -> FindObject("hist_ksl") -> Clone();
  hist_ksl_sc -> Scale(ksl_sfw);
  hist_ksl_sc -> SetName("hist_ksl_sc");
  Hlist_sc -> Add(hist_ksl_sc);
  
  TH1D *hist_etagam_sc = (TH1D *) Hlist -> FindObject("hist_etagam") -> Clone();
  hist_etagam_sc -> Scale(etagam_sfw);
  hist_etagam_sc -> SetName("hist_etagam_sc");
  Hlist_sc -> Add(hist_etagam_sc);
    
  TH1D *hist_mcrest_sc = (TH1D *) Hlist -> FindObject("hist_mcrest") -> Clone();
  hist_mcrest_sc -> Scale(mcrest_sfw);
  hist_mcrest_sc -> SetName("hist_mcrest_sc");
  Hlist_sc -> Add(hist_mcrest_sc);

  // Signal
  TH1D *hist_isr3pi_sc = (TH1D *) Hlist -> FindObject("hist_isr3pi") -> Clone();
  hist_isr3pi_sc -> Scale(isr3pi_sfw);
  hist_isr3pi_sc -> SetName("hist_isr3pi_sc");
  Hlist_sc -> Add(hist_isr3pi_sc);

  TH1D *hsig_sc = (TH1D *) Hlist -> FindObject("hsig") -> Clone();
  hsig_sc -> Scale(sig_sfw);
  hsig_sc -> SetName("hsig_sc");
  Hlist_sc -> Add(hsig_sc);

  // MC sum
  TH1D* hist_mcsum_sc = (TH1D*) hist_mcrest_sc -> Clone();
  hist_mcsum_sc -> Add(hist_etagam_sc, 1.);
  hist_mcsum_sc -> Add(hist_ksl_sc, 1.);
  hist_mcsum_sc -> Add(hist_omegapi_sc, 1.);
  hist_mcsum_sc -> Add(hist_eeg_sc, 1.);
  //hist_mcsum_sc -> Add(hsig_sc, 1.);
  hist_mcsum_sc -> Add(hist_isr3pi_sc, 1.);
  Hlist_sc -> Add(hist_mcsum_sc);
  
  // Data
  TH1D *hist_data = (TH1D *) Hlist -> FindObject("hist_data");
  
  
  hist_mcsum_sc -> SetName("hist_mcsum_sc");
  format_h(hist_mcsum_sc, 2, 2);

  // test plot
  hist_data -> Draw();
  hist_mcsum_sc -> Draw("SameHist");
    
  // save
  TFile *f_out = new TFile("../../hist.root", "recreate");

  hist_data -> Write();
  
  //Hlist -> Write();
  Hlist_sc -> Write();
    
  f_out -> Close();
  
  return 0;
  
}
