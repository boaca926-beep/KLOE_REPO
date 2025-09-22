#include "../header/sm_para.h"
#include "../header/path.h"
#include "../header/cut_para.h"
#include "../header/method.h"
#include <TStopwatch.h>

TStopwatch timer;
timer.Start();
  
int ntuples(){

  //const TString f_path = folder_in + "/" + data_type + ".root";
  cout << "input path: " << sampleFile << endl;

  TFile *f_input = new TFile(sampleFile + ".root");

  //checkFile(f_input);

  TTree *ALLCHAIN_CUT = (TTree*)f_input -> Get("ALLCHAIN_CUT");
  
  // define variables
  double lagvalue_min_7C = 0.;
  double deltaE = 0.;
  double Ephi_miss = 0.;
  double angle_pi0gam12 = 0.;  
  double betapi0 = 0.;
  double ppIM = 0.;
  double IM3pi_7C = 0., IM3pi_true = 0.;
  double IM_pi0_7C = 0.;
  double Eisr = 0., Epi0_pho1 = 0., Epi0_pho2 = 0.;
  double pull_E1 = 0.;
  double pull_x1 = 0.;
  double pull_y1 = 0.;
  double pull_z1 = 0.;
  double pull_t1 = 0.;
  double m02 = 0., mplus2 = 0.;
  double trkmass = 0.;
  double Epho_sum_recoil = 0.;
  
  int phid = 0, sig_type = 0;
  //int trigger_indx = 0;
  //int filfo_indx = 0;
  //int evtcls_indx = 0;
  //int fstate_indx = 0;
  int bkg_indx = 0, recon_indx = 0;

  double evnt_tot = 0; // total number of events
  double evnt_trigger = 0; // number of events after trigger
  double evnt_filfo = 0; // number of events after fiflo
  double evnt_evtcls = 0; // number of events after event classification
  double evnt_fstate = 0; // number of events after 2 tracks and 3 prompt photons

  double Eprompt_max = 0.;
  
  TFile *f_output = new TFile(outputCut + "tree_pre.root", "update");

  // ksl stream
  const int list_size = 11;
  const TString TNM[list_size] = {"TDATA", "TOMEGAPI", "TKPM", "TKSL", "T3PIGAM", "TRHOPI", "TETAGAM", "TBKGREST", "TUFO", "TEEG", "TISR3PI_SIG"};
  
  TTree *TTList[list_size];
  
  TCollection* tree_list = new TList;

  for (int i = 0; i < list_size; i ++) {// loop over the Tree list

    TTList[i] = new TTree(TNM[i], "recreate");
    TTList[i] -> SetAutoSave(0);
    
    tree_list -> Add(TTList[i]);
    
  } // end the Tree list
  
  // define branches
  TObject* treeout=0;
  TIter treeliter(tree_list);
  
  while((treeout=treeliter.Next()) != 0) {
    
    //treeout->Print();
    TTree* tree_tmp=dynamic_cast<TTree*>(treeout);

    tree_tmp -> Branch("Br_bkg_indx", &bkg_indx, "Br_bkg_indx/I");
    tree_tmp -> Branch("Br_recon_indx", &recon_indx, "Br_recon_indx/I");
    	  
  }
  

  for (Int_t irow = 0; irow < ALLCHAIN_CUT -> GetEntries(); irow ++) {// loop trees
	  
    ALLCHAIN_CUT -> GetEntry(irow);

    bkg_indx = ALLCHAIN_CUT -> GetLeaf("Br_bkg_indx") -> GetValue(0);
    recon_indx = ALLCHAIN_CUT -> GetLeaf("Br_recon_indx") -> GetValue(0);

    evnt_tot ++;

    //cout << deltaE << endl;

    TTList[0]-> Fill(); // data
    TTList[8]-> Fill(); // ufo
    TTList[9]-> Fill(); // eeg
    TTList[10]-> Fill(); // sig
    

    if (phid == 0) {// omegapi
      TTList[1] -> Fill();
    }
    else if (phid == 1) {// kpm
      TTList[2] -> Fill();
    }
    else if (phid == 2) {// ksl
      TTList[3] -> Fill();
    }
    else if (phid == 3) {
      if (sig_type == 1) {// 3pi
	TTList[4] -> Fill();
      }
      else {// rho pi
	TTList[5] -> Fill();
      }
    }
    else if (phid == 5) {
      if (sig_type == 1) {// etagam 3pi
	TTList[6] -> Fill();
      }
      else {// etagam rest
	TTList[7] -> Fill();
      }
    }
    else {// bkg rest
      TTList[7] -> Fill();
    }
  }

  // save
  if (data_type == "exp") {
    TTList[0]-> Write();
    cout << "TDATA saved" << endl;
  }
  else if (data_type == "ufo") {
    TTList[8]-> Write();
    cout << "TUFO saved" << endl;

    // Fill and save cut parameters
    //TTree* TPARA_CUT = new TTree("TPARA_CUT", "recreate");
    //TPARA_CUT -> Branch("Br_egammamin", &egammamin, "Br_egammamin/D");
    //TPARA_CUT -> Fill();
    //TPARA_CUT -> Write();
  }
  else if (data_type == "eeg") {
    TTList[9]-> Write("TEEG");
    cout << "TEEG saved" << endl;
  }
  else if (data_type == "sig") {
    TTList[10]-> Write("TISR3PI_SIG");
    cout << "TISR3PI_SIG saved" << endl;
  }
  else if (data_type == "ksl") {
    TTList[1] -> Write();
    TTList[2] -> Write();
    TTList[3] -> Write();
    TTList[4] -> Write();
    TTList[5] -> Write();
    TTList[6] -> Write();
    TTList[7] -> Write();
    cout << "All trees in KSL saved" << endl;
  }

    
  
  // Summary
  cout << "=========================================\n"
       << f_output -> GetName() << endl;
  
  cout << "evnt_tot = " << evnt_tot << "\n"
       << "2 tacks and 3 prompt photons: " << evnt_fstate << "\n";

  cout << "Cut parameters: \n"
       << "chi2_cut = " << chi2_cut << "\n"
       << "angle_cut = " << angle_cut << "\n"
       << "deltaE_cut = " << deltaE_cut << "\n"
       << "beta_cut = " << beta_cut << "\n";
       
  /// histos

  /// save
  f_output -> Close();

  double realTime = timer.RealTime();  // Wall clock time in seconds
  double expect_time = realTime / 60; // Expected time used in mins
  //cout << "Real time = " << realTime << endl;
  cout << "Expected time: " << expect_time << " mins" << endl;

  timer.Stop();
  timer.Print();
  
  return 0;
  
}

  
  
