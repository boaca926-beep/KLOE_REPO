#include "../header/fitfun.h"
#include "../header/plot.h"

gROOT->ForceStyle();  

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

TFile* f_cut = new TFile("../../input_norm_TDATA/cut/tree_pre.root");
TH1D* HPullList[100];
TObjArray *Hlist = new TObjArray(100);

void fillHist() {

  // data and MC background
  TIter next_tree(f_cut -> GetListOfKeys());

  TString objnm_tree, classnm_tree;

  const int hist_size = 4;
  const int bins = 200;
  const double x_min = -5.;
  const double x_max = 5.;
  
  TKey *key;

  int hist_indx = 0;
  
  while ( (key = (TKey *) next_tree() ) ) {// start tree while loop

    objnm_tree   =  key -> GetName();
    classnm_tree = key -> GetClassName();
    key -> GetSeekKey();

    //cout << "classnm = " << classnm_tree << ", objnm = " << objnm_tree << endl;

    TTree *tree_tmp = (TTree*)f_cut -> Get(objnm_tree);
    //cout << tree_tmp -> GetName() << endl;

    //TString name = TString::Format("hpull_%d", hist_indx);
    TString name = TString::Format("hist_%s", tree_tmp -> GetName());
    HPullList[hist_indx] = new TH1D(name, "", bins, x_min, x_max);
    cout << "Filling histo " << hist_indx + 1 << ", name: " << HPullList[hist_indx] -> GetName() << endl;

    hist_indx += 1;
    
  }
  
}

int Getpulls_compr() {

  gErrorIgnoreLevel = kError;
  TGaxis::SetMaxDigits(4);
  //gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  //gStyle->SetStatBorderSize(0);
  gStyle->SetFitFormat("6.3g");

  cout << "Get pull comparison from the kin. fit output ...\n";

  //inspect_input(f_cut);

  fillHist();
  
  return 0;

}
