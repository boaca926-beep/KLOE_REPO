#include "../hist.h"
#include "TFractionFitter.h"

int sfw1d_test() {

  //TFractionFitter
  //m3pi region [650, 900] MeV/c^{2}

  //Get trees
  TFile *intree_cut = new TFile("../crx3pi/output8373_finebin/tree_cut0.root");

  TIter next_tree(intree_cut -> GetListOfKeys());

  TString objnm_tree, classnm_tree;

  int i = 0;
  
  while ( (key = (TKey *) next_tree() ) ) {// start tree while lop
    
    i ++;
    
    objnm_tree   =  key -> GetName();
    classnm_tree = key -> GetClassName();
    key -> GetSeekKey();
    
    cout << "tree" << i << ": classnm = " << classnm_tree << ", objnm = " << objnm_tree << endl;
    
  }

  //Fill histos
  const int TL_size = 9;
  
  TTree * TreeList[TL_size] = {TDATA, TEEG, TOMEGAPI, TKSL, TKPM, TRHOPI, TETAGAM, TBKGREST, TISR3PI_SIG};

  int color_list[TL_size] = {1, 6, 7, 28, 46, 42, 3, 37, 4};

  double IM3pi = 0.;
  
  TH1D * H1D_IM3PI[TL_size];
  TH1D * H1D_IM3PI_TRUE[TL_size];

  TString tree_type = "";

  const double IM3pi_min = 650., IM3pi_max = 900.;
  const int IM3pi_bin = 200;
  
  for (int i = 0; i < TL_size; i ++) {// start MC type loop

    tree_type = TreeList[i] -> GetName();

    H1D_IM3PI[i] = new TH1D("h1d_IM3pi_" + tree_type, "", IM3pi_bin, IM3pi_min, IM3pi_max);
    H1D_IM3PI[i] -> Sumw2();

    for (Int_t irow = 0; irow < TreeList[i] -> GetEntries(); irow++) {// loop chain

      TreeList[i] -> GetEntry(irow);

      IM3pi = TreeList[i] -> GetLeaf("Br_IM3pi_7C") -> GetValue(0);

      H1D_IM3PI[i] -> Fill(IM3pi);

      

    }

    format_h(H1D_IM3PI[i], color_list[i], 2);
    
  }

  for (int i = 0; i < TL_size; i ++) {

    cout << H1D_IM3PI[i] -> GetName() << endl;

  }

  h1d_IM3pi_TEEG -> Scale(2.); // scale eeg with respect to lumi. scaling factor

  // merge mc bkgsum

  // Merge bkg sum
  TH1D * h1d_IM3pi_MCREST = (TH1D*) h1d_IM3pi_TBKGREST -> Clone();
  h1d_IM3pi_MCREST -> Add(h1d_IM3pi_TKPM, 1.);
  h1d_IM3pi_MCREST -> Add(h1d_IM3pi_TRHOPI, 1.);
  h1d_IM3pi_MCREST -> SetName("h1d_IM3pi_MCREST");

  // Merge MC sum
  TH1D * h1d_IM3pi_MCSUM = (TH1D*) h1d_IM3pi_TEEG -> Clone();
  h1d_IM3pi_MCSUM -> Add(h1d_IM3pi_TOMEGAPI, 1.);
  h1d_IM3pi_MCSUM -> Add(h1d_IM3pi_TKSL, 1.);
  h1d_IM3pi_MCSUM -> Add(h1d_IM3pi_TETAGAM, 1.);
  h1d_IM3pi_MCSUM -> Add(h1d_IM3pi_TISR3PI_SIG, 1.);
  h1d_IM3pi_MCSUM -> Add(h1d_IM3pi_MCREST, 1.);
  h1d_IM3pi_MCSUM -> SetName("h1d_IM3pi_MCSUM");

  format_h(h1d_IM3pi_MCSUM, 1, 2);
  
  //h1d_IM3pi_TDATA
  
  //h1d_IM3pi_TEEG
  //h1d_IM3pi_TOMEGAPI
  //h1d_IM3pi_TETAGAM
  //h1d_IM3pi_TKSL
  
  //h1d_IM3pi_TKPM
  //h1d_IM3pi_TRHOPI
  //h1d_IM3pi_TBKGREST

  //h1d_IM3pi_TISR3PI_SIG

  TObjArray *mctot = new TObjArray(6);
  // MC histograms are put in this array
  mctot -> Add(h1d_IM3pi_TISR3PI_SIG);
  mctot -> Add(h1d_IM3pi_TEEG);
  mctot -> Add(h1d_IM3pi_TOMEGAPI);
  mctot -> Add(h1d_IM3pi_TETAGAM);
  mctot -> Add(h1d_IM3pi_TKSL);
  mctot -> Add(h1d_IM3pi_MCREST);

  TFractionFitter* myfit = new TFractionFitter(h1d_IM3pi_TDATA, h1d_IM3pi_MCSUM);

  // initialize
  //myfit -> SetRangeX(1,15);
  //Int_t status = myfit -> Fit();
  // use only the first 15 bins in the fit
  // perform the fit
  
  // Plot histos

  h1d_IM3pi_TISR3PI_SIG -> Draw();
  h1d_IM3pi_TDATA -> Draw("Same");
  h1d_IM3pi_TEEG -> Draw("Same");
  h1d_IM3pi_TOMEGAPI -> Draw("Same");
  h1d_IM3pi_TETAGAM -> Draw("Same");
  h1d_IM3pi_TKSL -> Draw("Same");
  h1d_IM3pi_MCREST -> Draw("Same");
    
  //h1d_IM3pi_MCSUM -> Draw("Same");
  
  return 0;
  
}
