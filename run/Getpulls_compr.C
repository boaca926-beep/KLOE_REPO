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


int Getpulls_compr() {

  gErrorIgnoreLevel = kError;
  TGaxis::SetMaxDigits(4);
  //gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  //gStyle->SetStatBorderSize(0);
  gStyle->SetFitFormat("6.3g");

  cout << "Get pull comparison from the kin. fit output ...\n";

  //inspect_input(f_cut);

  // Create a list of input root files from folder named after branch names:
  // pull_E1, pull_x1

  const int list_size = 2;
  TString BR_LIST[list_size] = {"pull_E1", "pull_x1"};

  for (int i = 0; i < list_size; i++) {// Loop over branches

    TString f_nm = "../../pulls_compr/" + BR_LIST[i] + ".root";  
    cout << i << ": " << f_nm << endl;

    TFile* f_input = new TFile(f_nm);
  
    TIter next_tree1(f_input -> GetListOfKeys());

    TString objnm_tree, classnm_tree;

    int j = 0;
    TKey *key;
  
    while ( (key = (TKey *) next_tree1() ) ) {// Loop over trees
    
      j ++;
    
      objnm_tree   =  key -> GetName();
      classnm_tree = key -> GetClassName();
      //key -> GetSeekKey();
      
      cout << "tree" << j << ": classnm = " << classnm_tree << ", objnm = " << objnm_tree << endl;
    
    } // end loop trees
    
  }// end loop files
  return 0;

}
