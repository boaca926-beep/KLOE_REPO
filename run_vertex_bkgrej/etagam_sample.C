TLorentzVector Get4vector(double E, double px, double py, double pz) {

  //given a cluster index returns the 4-mom of a photon
  TVector3 gamma(px, py, pz);
  
  Double_t scale1;
  scale1=E/gamma.Mag();
  TLorentzVector gamma4mom(scale1*gamma, E);
  //cout << gamma4mom.M() << endl;
  return gamma4mom;

}

int etagam_sample() {
  // from e+ e- -> eta gamma get 3 photon final state 4-momentum

  TString sampleFile = "/home/bo/Desktop/analysis_root_v6/ksl.root";

  cout << "Input root file: " << sampleFile << endl;
  
  TFile *f_input = new TFile(sampleFile);

  TTree *PHOTONS = (TTree*)f_input -> Get("PHOTONS");

  TFile *f_output = new TFile("../../etagam.root", "recreate");

  const int list_size = 1;
  const TString TNM[list_size] = {"TETAGAM"};

  TTree *TTList[list_size];
  
  TCollection* tree_list = new TList;

  for (int i = 0; i < list_size; i ++) {// loop over the Tree list

    TTList[i] = new TTree(TNM[i], "recreate");
    TTList[i] -> SetAutoSave(0);
    
    tree_list -> Add(TTList[i]);
    
  } // end the Tree list

  // define variables
  double pho_E1 = 0., pho_px1 = 0., pho_py1 = 0., pho_pz1 = 0.;
  double pho_E2 = 0., pho_px2 = 0., pho_py2 = 0., pho_pz2 = 0.;
  double pho_E3 = 0., pho_px3 = 0., pho_py3 = 0., pho_pz3 = 0.;
  double mpi0 = 0., mpi0_good = 0., mpi0_bad = 0.;
  
  int phid = -1, sig_type = -1;
  
  // define branches
  TObject* treeout=0;
  TIter treeliter(tree_list);
  
  while((treeout=treeliter.Next()) != 0) {
    
    //treeout->Print();
    TTree* tree_tmp=dynamic_cast<TTree*>(treeout);

    tree_tmp -> Branch("Br_E1", &pho_E1, "Br_E1/D");
    tree_tmp -> Branch("Br_px1", &pho_px1, "Br_px1/D");
    tree_tmp -> Branch("Br_py1", &pho_py1, "Br_py1/D");
    tree_tmp -> Branch("Br_pz1", &pho_pz1, "Br_pz1/D");

    tree_tmp -> Branch("Br_E2", &pho_E2, "Br_E2/D");
    tree_tmp -> Branch("Br_px2", &pho_px2, "Br_px2/D");
    tree_tmp -> Branch("Br_py2", &pho_py2, "Br_py2/D");
    tree_tmp -> Branch("Br_pz2", &pho_pz2, "Br_pz2/D");

    tree_tmp -> Branch("Br_mpi0", &mpi0, "Br_mpi0/D");
    tree_tmp -> Branch("Br_mpi0_good", &mpi0_good, "Br_mpi0_good/D");
    tree_tmp -> Branch("Br_mpi0_bad", &mpi0_bad, "Br_mpi0_bad/D");

    tree_tmp -> Branch("Br_phid", &phid, "Br_phid/I");
    tree_tmp -> Branch("Br_sig_type", &sig_type, "Br_sig_type/I");
    	  
  }

  TLorentzVector pi0gam1, pi0gam2, isrgam, trkplus, trkmin;
  
  for (Int_t irow = 0; irow < PHOTONS -> GetEntries(); irow ++) {// loop trees
	  
    PHOTONS -> GetEntry(irow);

    phid = PHOTONS -> GetLeaf("Br_phid") -> GetValue(0);
    sig_type = PHOTONS -> GetLeaf("Br_sig_type") -> GetValue(0);

    pho_E1 = PHOTONS -> GetLeaf("Br_E1") -> GetValue(0);
    pho_px1 = PHOTONS -> GetLeaf("Br_px1") -> GetValue(0);
    pho_py1 = PHOTONS -> GetLeaf("Br_py1") -> GetValue(0);
    pho_pz1 = PHOTONS -> GetLeaf("Br_pz1") -> GetValue(0);

    pho_E2 = PHOTONS -> GetLeaf("Br_E2") -> GetValue(0);
    pho_px2 = PHOTONS -> GetLeaf("Br_px2") -> GetValue(0);
    pho_py2 = PHOTONS -> GetLeaf("Br_py2") -> GetValue(0);
    pho_pz2 = PHOTONS -> GetLeaf("Br_pz2") -> GetValue(0);

    pi0gam1 = Get4vector(pho_E1, pho_px1, pho_py1, pho_pz1);
    pi0gam2 = Get4vector(pho_E2, pho_px2, pho_py2, pho_pz2);

    //cout << pi0gam1.E() << ", " << endl;
    
    // pi0 invariant mass, identified pi0 photons
    mpi0 = (pi0gam1 + pi0gam2).M();

    //cout << "pho_E1 = " << pho_E1 << endl;
      
    // fill trees for eta gamma
    if (!TMath::IsNaN(mpi0)) {

      if (phid == 5 && sig_type == 1 ) {// good events
	mpi0_good = mpi0;
	//cout << mpi0_good << endl;
      }
      else {// bad
	mpi0_bad = mpi0
      }
    
    }

    TTList[0]-> Fill();

  }

  // save
  TTList[0]-> Write();
      
  return 0;
  
}
