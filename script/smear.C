// using pi0 invariant mass to smear photon 4-momentum resolution
// input sample:

#include "../header/plot.h"

const TString infile_nm="/home/bo/Desktop/analysis_root_v6/crx3pi_norm/output/tree_cut0.root";

const double IM3pi_min = 770;
const double IM3pi_max = 800;
const double IM3pi_sigma = 2.615;
const double mass_sigma_nb = 0.25; //0.25
const int IM3pi_bin = 60; //TMath::Nint((IM3pi_max - IM3pi_min) / mass_sigma_nb / IM3pi_sigma);

int smear() {

  //gROOT->SetBatch(kTRUE);  
  gErrorIgnoreLevel = kError;
  TGaxis::SetMaxDigits(3);
  //gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  // switch on histogram errors
  TH1::SetDefaultSumw2();

  /// Inspect the input tree
  TFile* intree = new TFile(infile_nm);
  
  TIter next_tree(intree -> GetListOfKeys());

  TString objnm_tree, classnm_tree;

  int i = 0;
  TKey *key;
  
  while ( (key = (TKey *) next_tree() ) ) {
    
    i ++;
    
    objnm_tree   =  key -> GetName();
    classnm_tree = key -> GetClassName();
    key -> GetSeekKey();
    
    //cout << "tree" << i << ": classnm = " << classnm_tree << ", objnm = " << objnm_tree << endl;
    
  }

  /// Fill histos
  
  cout << "IM3pi range (" << IM3pi_min << ", " << IM3pi_max << ") MeV/c^2 \n"
       << "bins = " << IM3pi_bin << "\n";

  TTree *TDATA = (TTree*)intree -> Get("TDATA");
  TTree *TEEG = (TTree*)intree -> Get("TEEG");
  TTree *TOMEGAPI = (TTree*)intree -> Get("TOMEGAPI");
  TTree *TKSL = (TTree*)intree -> Get("TKSL");
  TTree *TKPM = (TTree*)intree -> Get("TKPM");
  TTree *TRHOPI = (TTree*)intree -> Get("TRHOPI");
  TTree *TETAGAM = (TTree*)intree -> Get("TETAGAM");
  TTree *TBKGREST = (TTree*)intree -> Get("TBKGREST");
  TTree* TISR3PI_SIG = (TTree*)intree -> Get("TISR3PI_SIG");

  // signal
  TH1D * hIM3pi_signal_true = new TH1D("hIM3pi_signal_true", "", IM3pi_bin, IM3pi_min, IM3pi_max);

  double IM3pi_true = 0., IM3pi_rec = 0.;
  
  for (Int_t irow = 0; irow < TISR3PI_SIG -> GetEntries(); irow++) {// loop chain
    
    TISR3PI_SIG -> GetEntry(irow);

    IM3pi_true = TISR3PI_SIG -> GetLeaf("Br_IM3pi_true") -> GetValue(0);
   
    // filling histos

    hIM3pi_signal_true -> Fill(IM3pi_true);

    if (irow > 1e4) break;
   
  }

  // attributes
  format_h(hIM3pi_signal_true, 4, 2);
  
  hIM3pi_signal_true -> Draw();
  
  // backgrounds
  TObjArray *Hlist = new TObjArray();
  
  const int TLSize = 8;
  
  TTree * TrList[TLSize] = {TDATA, TEEG, TOMEGAPI, TKSL, TKPM, TRHOPI, TETAGAM, TBKGREST};
  const TString TrNm[TLSize] = {"data", "eeg", "omegapi" , "ksl", "kpm", "rhopi" , "etagam" , "bkgrest"};

  int color_list[TLSize] = {1, 6, 7, 28, 46, 42, 3, 37};

  TH1D* h;

  IM3pi_rec = 0.;
  
  for (int i = 0; i < TLSize; i ++) {// start MC type loop

    h = new TH1D("hIM3pi_" + TrNm[i], "", IM3pi_bin, IM3pi_min, IM3pi_max);
    //h -> Sumw2();
    
    for (Int_t irow = 0; irow < TrList[i] -> GetEntries(); irow++) {
      
      TrList[i] -> GetEntry(irow);

      IM3pi_rec = TrList[i] -> GetLeaf("Br_IM3pi_7C") -> GetValue(0);

      h -> Fill(IM3pi_rec);

      if (irow > 1e3) break;
      
    }

    
    format_h(h, color_list[i], 2);

    Hlist->Add(h);
    
  }

  //TH1D *hIM3pi_bkgrest = (TH1D *) Hlist -> FindObject("hIM3pi_bkgrest");
  //TH1D *hIM3pi_kpm = (TH1D *) Hlist -> FindObject("hIM3pi_kpm");
  //TH1D *hIM3pi_rhopi = (TH1D *) Hlist -> FindObject("hIM3pi_rhopi");
  //TH1D *hIM3pi_eeg = (TH1D *) Hlist -> FindObject("hIM3pi_eeg");
  //TH1D *hIM3pi_omegapi = (TH1D *) Hlist -> FindObject("hIM3pi_omegapi");
  //TH1D *hIM3pi_etagam = (TH1D *) Hlist -> FindObject("hIM3pi_etagam");
  //TH1D *hIM3pi_ksl = (TH1D *) Hlist -> FindObject("hIM3pi_ksl");
  TH1D *hIM3pi_data = (TH1D *) Hlist -> FindObject("hIM3pi_data");

  hIM3pi_data -> Draw();
  
  return 0;
  
}
