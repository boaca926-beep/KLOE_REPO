#define MyClass_cxx
#include "../header/MyClass.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void MyClass::Main()
{
  //cout << "PROCESSING ..." << endl;
  int sig_type = 999;
  int bkg_indx = 999;
  int recon_indx = 999;
  int evtcls_indx = 999;
  int trigger_indx = 999;
  int filfo_indx = 999;
  int test_indx = 0;
  int fstate_indx = 999;
  
  double evnt_tot = 0, evnt_pre_filtred = 0, evnt_pre_pass = 0;
  double evnt_unstr = 0, evnt_cls = 0;
  double evnt_sig = 0, evnt_bkg = 0;
  double evnt_pre = 0;
  double evnt_str2 = 0, evnt_no_str2 = 0, evnt_class = 0;
  double evnt_trigger = 0, evnt_no_trigger = 0;
  double evnt_filfo = 0, evnt_no_filfo = 0;
  double evnt_fstate = 0, evnt_no_fstate = 0;
  double evnt_cutted = 0;
  
  double evnt_trk = 0, evnt_final = 0.;
  double evnt_phi5 = 0;
  double evnt_recon_type0 = 0., evnt_recon_type1 = 0., evnt_recon_type2 = 0.;
  // reconstruction
  int trknb_prompt = 0, trkindx1 = 0, trkindx2 = 0;
  int promptnb = 0;
  // mass
  double IM_3pi = 0., IM_3pi_nofit = 0., IM3pi_7C = 0., IM3pi_true = 0., IM3pi = 0., IM3pi_beta = 0.;
  double IM_pi0 = 0., IM_pi0_nofit = 0., IM_pi0_7C = 0.;
  double trkmass_true = 0., trkmass = 0.;
  double ppIM_true = 0., ppIM = 0., ppIM_beta = 0.;
  double IMisrpho_miss = 0.;
  double IM3pi_pi12 = 0., IM3pi_pi13 = 0., IM3pi_pi23 = 0.; 
  double MASSLIST[100];
  double ANGLELIST[100];
  double PULLIST[100];
  double RESOLIST[100];
  double MOM4PHO1[4];
  double MOM4PHO2[4];
  double MOM4PHO3[4];
  double MOM4TRKPLUS[4];
  double MOM4TRKMINS[4];
  
  // energy
  double ENERGYLIST[100];
  double PI0PHORESD[100];
  double E_pho_isr = 0., Emax_pho = 0.;
  double EPI0GAM[2], EPI0NTMC[2], E_pi0gam1 = 0., E_pi0gam2 = 0.;
  double deltaE_true = 0.;
  double Emax_clust = 0., Esum_clust = 0.;
  double Esum = 0., E_radiv1 = 0., E_radiv2 = 0.;
  double deltaE = 0.;
  // angle
  double Angle_pho_isr = 0.;
  double Angle_clust = 0.;
  double angle_pi0gam12_true = 0.;
  double angle_ppl_3piboost = 0.;
  double minangle = 23.;
  //double egammamin = 5.;
  //double egammamin = 10.;
  //double egammamin = 15.; // minumum photon energy cut, standard 15 MeV 
  //double egammamin = 20.;
  //double egammamin = 25.;
  
  // time
  double Tof_clust = 0.;
  double Sigma2_T_clust = 0., Sigma_T_clust = 0.;
  // kinfit
  const double xyres = 1.2, zres = 1.4, tres = 0.057, tcst = 0.140, eres = 0.057;
  // statistics
  double lagvalue_min_7C = 0., pvalue = 0.;
  double chi2mgg_min_7C = 0.;
  // others
  double tclusdiff = 0.;
  double beta_clu = 0.;
  double betapi0 = 0., betapi0_true = 0.;
  //
  double test_value = 0.;
  
  // define tree
  TTree ALLCHAIN_GEN("ALLCHAIN_GEN", "recreate"); ALLCHAIN_GEN.SetAutoSave(0);

  ALLCHAIN_GEN.Branch("Br_sig_type", &sig_type, "Br_sig_type/I");
  ALLCHAIN_GEN.Branch("Br_sel_type", &bit_select, "Br_sel_type/I");
  ALLCHAIN_GEN.Branch("Br_IM_3pi", &IM_3pi, "Br_IM_3pi/D");
  ALLCHAIN_GEN.Branch("Br_E_pho_isr", &E_pho_isr, "Br_E_pho_isr/D");
  ALLCHAIN_GEN.Branch("Br_Angle_pho_isr", &Angle_pho_isr, "Br_Angle_pho_isr/D");
  ALLCHAIN_GEN.Branch("Br_bpx", &bpx, "Br_bpx/F");
  ALLCHAIN_GEN.Branch("Br_Esum", &Esum, "Br_Esum/D");

  TTree ALLCHAIN_STR2("ALLCHAIN_STR2", "recreate"); ALLCHAIN_STR2.SetAutoSave(0);
  
  ALLCHAIN_STR2.Branch("Br_sig_type", &sig_type, "Br_sig_type/I");
  ALLCHAIN_STR2.Branch("Br_sel_type", &bit_select, "Br_sel_type/I");
  ALLCHAIN_STR2.Branch("Br_IM_3pi", &IM_3pi, "Br_IM_3pi/D");
  
  TTree TrSample ("TrSample", "recreate"); TrSample.SetAutoSave(0);
  TrSample.Branch("Br_IM_3pi", &IM_3pi, "Br_IM_3pi/D");
  TrSample.Branch("Br_E_pho_isr", &E_pho_isr, "Br_E_pho_isr/D");
  TrSample.Branch("Br_Angle_pho_isr", &Angle_pho_isr, "Br_Angle_pho_isr/D");
  TrSample.Branch("Br_bpx", &bpx, "Br_bpx/F");
  TrSample.Branch("Br_Esum", &Esum, "Br_Esum/D");

  //
  TTree ALLCHAIN_CUT ("ALLCHAIN_CUT", "recreate"); ALLCHAIN_CUT.SetAutoSave(0);
  ALLCHAIN_CUT.Branch("Br_evtcls_indx", &evtcls_indx, "Br_evtcls_indx/I");
  ALLCHAIN_CUT.Branch("Br_sig_type", &sig_type, "Br_sig_type/I");
  ALLCHAIN_CUT.Branch("Br_sel_type", &bit_select, "Br_sel_type/I");
  ALLCHAIN_CUT.Branch("Br_phid", &phid, "Br_phid/I");
  ALLCHAIN_CUT.Branch("Br_kineid", &kineid, "Br_kineid/I");
  ALLCHAIN_CUT.Branch("Br_bkg_indx", &bkg_indx, "Br_bkg_indx/I");
  ALLCHAIN_CUT.Branch("Br_recon_indx", &recon_indx, "Br_recon_indx/I");
  ALLCHAIN_CUT.Branch("Br_trigger_indx", &trigger_indx, "Br_trigger_indx/I");
  ALLCHAIN_CUT.Branch("Br_evtcls_indx", &evtcls_indx, "Br_evtcls_indx/I");
  ALLCHAIN_CUT.Branch("Br_filfo_indx", &filfo_indx, "Br_filfo_indx/I");
  ALLCHAIN_CUT.Branch("Br_fstate_indx", &fstate_indx, "Br_fstate_indx/I");
  
  
  //
  
  ALLCHAIN_CUT.Branch("Br_IM_3pi", &IM_3pi, "Br_IM_3pi/D");
  ALLCHAIN_CUT.Branch("Br_IM3pi_7C", &IM3pi_7C, "Br_IM3pi_7C/D");
  ALLCHAIN_CUT.Branch("Br_IM_3pi_nofit", &IM_3pi_nofit, "Br_IM_3pi_nofit/D");  
  ALLCHAIN_CUT.Branch("Br_IM3pi_pi12", &IM3pi_pi12, "Br_IM3pi_pi12/D"); 
  ALLCHAIN_CUT.Branch("Br_IM3pi_pi13", &IM3pi_pi13, "Br_IM3pi_pi13/D"); 
  ALLCHAIN_CUT.Branch("Br_IM3pi_pi23", &IM3pi_pi23, "Br_IM3pi_pi23/D"); 
  ALLCHAIN_CUT.Branch("Br_IM3pi_true", &IM3pi_true, "Br_IM3pi_true/D"); 
  ALLCHAIN_CUT.Branch("Br_IM_pi0", &IM_pi0, "Br_IM_pi0/D");
  ALLCHAIN_CUT.Branch("Br_IM_pi0_nofit", &IM_pi0_nofit, "Br_IM_pi0_nofit/D");  
  ALLCHAIN_CUT.Branch("Br_IM_pi0_7C", &IM_pi0_7C, "Br_IM_pi0_7C/D");
  ALLCHAIN_CUT.Branch("Br_trkmass_true", &trkmass_true, "Br_trkmass_true/D");
  ALLCHAIN_CUT.Branch("Br_trkmass", &trkmass, "Br_trkmass/D");
  ALLCHAIN_CUT.Branch("Br_ppIM_true", &ppIM_true, "Br_ppIM_true/D");
  ALLCHAIN_CUT.Branch("Br_ppIM", &ppIM, "Br_ppIM/D");
  ALLCHAIN_CUT.Branch("Br_ppIM_beta", &ppIM_beta, "Br_ppIM_beta/D");
  ALLCHAIN_CUT.Branch("Br_E_pho_isr", &E_pho_isr, "Br_E_pho_isr/D");
  ALLCHAIN_CUT.Branch("Br_deltaE_true", &deltaE_true, "Br_deltaE_true/D");
  ALLCHAIN_CUT.Branch("Br_Angle_pho_isr", &Angle_pho_isr, "Br_Angle_pho_isr/D"); 
  ALLCHAIN_CUT.Branch("Br_angle_pi0gam12_true", &angle_pi0gam12_true, "Br_angle_pi0gam12_true/D");
  ALLCHAIN_CUT.Branch("Br_Emax_clust", &Emax_clust, "Br_Emax_clust/D");
  ALLCHAIN_CUT.Branch("Br_Esum_clust", &Esum_clust, "Br_Esum_clust/D");

  //
  ALLCHAIN_CUT.Branch("Br_egammamin", &egammamin, "Br_egammamin/D");  
  ALLCHAIN_CUT.Branch("Br_Zv_trk1", &Zv_trk1, "Br_Zv_trk1/D");
  ALLCHAIN_CUT.Branch("Br_Zv_trk2", &Zv_trk2, "Br_Zv_trk2/D");
  ALLCHAIN_CUT.Branch("Br_Zvmin", &Zvmin, "Br_Zvmin/D");
 
  ALLCHAIN_CUT.Branch("Br_Rhov_trk1", &Rhov_trk1, "Br_Rhov_trk1/D");
  ALLCHAIN_CUT.Branch("Br_Rhov_trk2", &Rhov_trk2, "Br_Rhov_trk2/D");
  ALLCHAIN_CUT.Branch("Br_Rhovmin", &Rhovmin, "Br_Rhovmin/D");
  
  //zv[trkindx1]

  //
  ALLCHAIN_CUT.Branch("Br_lagvalue_min_7C", &lagvalue_min_7C, "Br_lagvalue_min_7C/D");
  ALLCHAIN_CUT.Branch("Br_pvalue", &pvalue, "Br_pvalue/D");
  ALLCHAIN_CUT.Branch("Br_betapi0", &betapi0, "Br_betapi0/D");
  ALLCHAIN_CUT.Branch("Br_betapi0_true", &betapi0_true, "Br_betapi0_true/D");
  ALLCHAIN_CUT.Branch("Br_test", &test_value, "Br_test/D");
  ALLCHAIN_CUT.Branch("Br_chi2mgg_min_7C", &chi2mgg_min_7C, "Br_chi2mgg_min_7C/D");
  ALLCHAIN_CUT.Branch("Br_beta_clu", &beta_clu, "Br_beta_clu/D");
  ALLCHAIN_CUT.Branch("Br_tclusdiff", &tclusdiff, "Br_tclusdiff/D");
  
  //
  ALLCHAIN_CUT.Branch("Br_ENERGYLIST", &ENERGYLIST, "Br_ENERGYLIST[100]/D");
  //ENERGYLIST[0]: ISR photon energy, 7C kin. fit
  
  ALLCHAIN_CUT.Branch("Br_PI0PHORESD", &PI0PHORESD, "Br_PI0PHORESD[100]/D");
  
  ALLCHAIN_CUT.Branch("Br_MASSLIST", &MASSLIST, "Br_MASSLIST[100]/D");
  //MASSLIST[0]: pi0 mass no fit
  //MASSLIST[1]: 3pi mass no fit
  //MASSLIST[2]: pi0 mass 7C kin. fit
  //MASSLIST[3]: 3pi mass 7C kin. fit
  //MASSLIST[4]: track mass
  //MASSLIST[5]: pi+ pi- invariant mass

  ALLCHAIN_CUT.Branch("Br_ANGLELIST", &ANGLELIST, "Br_ANGLELIST[100]/D");
  //ANGLELIST[0]: pi0 gamma12
  
  ALLCHAIN_CUT.Branch("Br_RESOLIST", &RESOLIST, "Br_RESOLIST[100]/D");
  //RESOLIST[0]: pi0 mass 7C kin. fit
  //RESOLIST[1]: 3pi mass 7C kin. fit
  //RESOLIST[2]: track mass 7C kin. fit
  //RESOLIST[3]: pi+ pi- invariant mass
  //RESOLIST[4]: betapi0

  ALLCHAIN_CUT.Branch("Br_PULLIST", &PULLIST, "Br_PULLIST[100]/D");

  //ALLCHAIN_CUT.Branch("Br_MOM4TRKPLUS", &MOM4TRKPLUS, "MOM4TRKPLUS[4]/D");

  //
  TTree ALLCHAIN_TEST ("ALLCHAIN_TEST", "recreate"); ALLCHAIN_TEST.SetAutoSave(0);
  ALLCHAIN_TEST.Branch("Br_bkg_indx", &bkg_indx, "Br_bkg_indx/I");
  ALLCHAIN_TEST.Branch("Br_IM3pi_7C", &IM3pi_7C, "Br_IM3pi_7C/D");
  ALLCHAIN_TEST.Branch("Br_IMisrpho_miss", &IMisrpho_miss, "Br_IMisrpho_miss/D");
  ALLCHAIN_TEST.Branch("Br_angle_ppl_3piboost", &angle_ppl_3piboost, "Br_angle_ppl_3piboost/D");
  
  
  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntriesFast();
  
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    // if (Cut(ientry) < 0) continue;

    // check MC or Data
    //cout << "kineid = " << kineid << endl;
    //if (kineid != -999) continue;
    
    evnt_tot ++;

    if (bit_select == 1) {// pre-filtred
      evnt_pre_filtred ++;
    }
    else {
      evnt_pre_pass ++;
    }

    /// beam
    Beam.SetPxPyPzE(bpx, bpy, bpz, bene);
    //TVector3 Boost3vector = -Beam.BoostVector();

    //cout << Beam.M() << endl;

    //cout << "beam (bpx, bpy, bpz, bene) = " << "(" << bpx << ", " << bpy << ", " << bpz << ", " << bene << ") \n";

    /// define isr 3pi signals
    double nb_pho_radiv = 0, nb_pi = 0, nb_pi0pho = 0;
    int pi0gam1_ntmc = 0, pi0gam2_ntmc = 0; 
    TVector3 MC_vect;
    TLorentzVector pi0MC_TLvect, piplusMC_TLvect, piminusMC_TLvect, threepi_TLvect, isrpho_TLvect, finalstate_TLvect, pho_radiv1_TLvect, pho_radiv2_TLvect, pi0gam1_TLvect, pi0gam2_TLvect;
   
      
    /// paticle infomation
    //cout << "\n\n" << endl;
    for (int i = 0; i < ntmc; i ++) {// loop over emc
      MC_vect.SetXYZ(pxmc[i], pymc[i], pzmc[i]);

      //if (phid == 5) {
      //evnt_phi5 ++;
      //cout << i << ": kineid = " << kineid << ", phid = " << phid << ", pidmc = " << pidmc[i] << ", virmom = " << virmom[i] << ", mother = " << mother[indv[i] - 1] << ", pxmc = " << pxmc[i] << ", pymc = " << pymc[i] << ", pzmc = " << pzmc[i] << "\n";
      //}
      
      // radiative photons (isr str3)
      if (pidmc[i] == 1 && virmom[i] == 1 && mother[indv[i] - 1] == 50 && nb_pho_radiv == 0) {// first radiative photon
	nb_pho_radiv ++;
	pho_radiv1_TLvect = GetLorentzVector(MC_vect, 0.);
      }
      else if (pidmc[i] == 1 && virmom[i] == 1 && mother[indv[i] - 1] == 50 && nb_pho_radiv == 1) {// second radiative photon
	nb_pho_radiv ++;
	pho_radiv2_TLvect = GetLorentzVector(MC_vect, 0.);
      }
      
      // radiative photons (isr phok5)
      else if (pidmc[i] == 1 && virmom[i] == 0 && mother[indv[i] - 1] == 50 && nb_pho_radiv == 0) {// first radiative photon
	nb_pho_radiv ++;
	pho_radiv1_TLvect = GetLorentzVector(MC_vect, 0.);
      }
      else if (pidmc[i] == 1 && virmom[i] == 0 && mother[indv[i] - 1] == 50 && nb_pho_radiv == 1) {// second radiative photon
	nb_pho_radiv ++;
	pho_radiv2_TLvect = GetLorentzVector(MC_vect, 0.);
      }
      
      // pions (isr str3)
      else if (pidmc[i] == 7 && virmom[i] == 54 && mother[indv[i] - 1] == 50) {// pi0
	nb_pi ++;
	pi0MC_TLvect = GetLorentzVector(MC_vect, massneupion);
	//cout << "vect.X = " << MC_vect.X() << endl;
	//cout << "TLvect.X = " << pi0MC_TLvect.X() << endl;
      }
      else if (pidmc[i] == 8 && virmom[i] == 54 && mother[indv[i] - 1] == 50) {// piplus
	nb_pi ++;
	piplusMC_TLvect = GetLorentzVector(MC_vect, masschpion);
      }
      else if (pidmc[i] == 9 && virmom[i] == 54 && mother[indv[i] - 1] == 50) {// piminus
	nb_pi ++;
	//cout << "vect.X = " << MC_vect.X() << endl;
	piminusMC_TLvect = GetLorentzVector(MC_vect, masschpion);
      }

      // pions (isr phok5)
      else if (pidmc[i] == 7 && virmom[i] == 0 && mother[indv[i] - 1] == 50) {// pi0
	nb_pi ++;
	pi0MC_TLvect = GetLorentzVector(MC_vect, massneupion);
	//cout << "vect.X = " << MC_vect.X() << endl;
	//cout << "TLvect.X = " << pi0MC_TLvect.X() << endl;
      }
      else if (pidmc[i] == 8 && virmom[i] == 0 && mother[indv[i] - 1] == 50) {// piplus
	nb_pi ++;
	piplusMC_TLvect = GetLorentzVector(MC_vect, masschpion);
      }
      else if (pidmc[i] == 9 && virmom[i] == 0 && mother[indv[i] - 1] == 50) {// piminus
	nb_pi ++;
	//cout << "vect.X = " << MC_vect.X() << endl;
	piminusMC_TLvect = GetLorentzVector(MC_vect, masschpion);
      }

      // pions (eta->3pi)
      else if (pidmc[i] == 7 && virmom[i] == 17 && mother[indv[i] - 1] == 50) {// pi0
	nb_pi ++;
	pi0MC_TLvect = GetLorentzVector(MC_vect, massneupion);
	//cout << "vect.X = " << MC_vect.X() << endl;
	//cout << "TLvect.X = " << pi0MC_TLvect.X() << endl;
      }
      else if (pidmc[i] == 8 && virmom[i] == 17 && mother[indv[i] - 1] == 50) {// piplus
	nb_pi ++;
	piplusMC_TLvect = GetLorentzVector(MC_vect, masschpion);
      }
      else if (pidmc[i] == 9 && virmom[i] == 17 && mother[indv[i] - 1] == 50) {// piminus
	nb_pi ++;
	//cout << "vect.X = " << MC_vect.X() << endl;
	piminusMC_TLvect = GetLorentzVector(MC_vect, masschpion);
      }
      
      // pi0 photons
      else if (pidmc[i] == 1 && virmom[i] == 0 && mother[indv[i] - 1] == 7 && nb_pi0pho == 0) {
	nb_pi0pho ++;
	pi0gam1_TLvect = GetLorentzVector(MC_vect, 0.);
	pi0gam1_ntmc = i;
	//cout << pxmc[i] << endl;
	
      }
      else if (pidmc[i] == 1 && virmom[i] == 0 && mother[indv[i] - 1] == 7 && nb_pi0pho == 1) {
	nb_pi0pho ++;
	pi0gam2_TLvect = GetLorentzVector(MC_vect, 0.);
	pi0gam2_ntmc = i;

	//cout << pxmc[i] << endl;
	
      }

    }// end loop emc

    /// MC info
    finalstate_TLvect = piminusMC_TLvect + piplusMC_TLvect + pi0MC_TLvect + pho_radiv1_TLvect + pho_radiv2_TLvect;
    threepi_TLvect = piminusMC_TLvect + piplusMC_TLvect + pi0gam1_TLvect + pi0gam2_TLvect;
    //threepi_TLvect = piminusMC_TLvect + piplusMC_TLvect + pi0MC_TLvect;
    isrpho_TLvect = Beam - (pi0MC_TLvect + piplusMC_TLvect + piminusMC_TLvect);

    IM_3pi = threepi_TLvect.M();//
    //IM_3pi = (Beam - (pho_radiv1_TLvect + pho_radiv2_TLvect)).M(); //threepi_TLvect.M();
    IM_pi0 = (pi0gam1_TLvect + pi0gam2_TLvect).M();//(Beam - (pho_radiv1_TLvect + pho_radiv2_TLvect)).M(); //threepi_TLvect.M();
    ppIM_true = (piplusMC_TLvect + piminusMC_TLvect).M();
    
    E_radiv1 = pho_radiv1_TLvect.E();
    E_radiv2 = pho_radiv2_TLvect.E();
    E_pho_isr = isrpho_TLvect.E();
    Esum = finalstate_TLvect.E();
    EPI0GAM[0] = pi0gam1_TLvect.E();
    EPI0GAM[1] = pi0gam2_TLvect.E();
    EPI0NTMC[0] = pi0gam1_ntmc; 
    EPI0NTMC[1] = pi0gam2_ntmc;
    E_pi0gam1 = pi0gam1_TLvect.E();
    E_pi0gam2 = pi0gam2_TLvect.E();
    betapi0_true = (pi0MC_TLvect.Vect()).Mag() / pi0MC_TLvect.E();
    deltaE_true = DeltaE(piplusMC_TLvect, piminusMC_TLvect);
    trkmass_true = Trkmass(piplusMC_TLvect, piminusMC_TLvect); 
    Angle_pho_isr = TMath::Cos(isrpho_TLvect.Theta());
    angle_pi0gam12_true = pi0gam1_TLvect.Angle(pi0gam2_TLvect.Vect())*TMath::RadToDeg();
    //cout << "deltaE_true = " << deltaE_true << endl;
    //cout << "IM_3pi = " << IM_3pi << ", E_pho_isr = " << E_pho_isr << ", Angle_pho_isr = " << Angle_pho_isr << "\n";
    //cout << "nb_pho_radiv = " << nb_pho_radiv << ", nb_pi = " << nb_pi << ", nb_pi0pho = " << nb_pi0pho << "\n";
    //cout << "pi0gam1 E = " << E_pi0gam1 << ", pi0gam2 E = " << E_pi0gam2 << endl;
      
    /// define ISR signal

    if (nb_pho_radiv > 0 && nb_pi == 3 && nb_pi0pho == 2) {
      sig_type = 1;
      evnt_sig ++;
      
      //cout << phid << endl;
      //cout << nb_pho_radiv << endl;
      //cout << pho_radiv1_TLvect.E() << "," << pho_radiv2_TLvect.E() << ", checked " << E_radiv1 << ", " << E_radiv2 << "\n";
      
      //cout << "jentry = " << jentry << ", finalstate E = " << finalstate_TLvect.E() << endl;
      //if (phid == 5) cout << "nb_pho_radiv = " << nb_pho_radiv << ", nb_pi = " << nb_pi << ", nb_pi0pho = " << nb_pi0pho << endl;
      //TrSample.Fill();
    }
    else {
      sig_type = 0;
      evnt_bkg ++;
    }
    
    // check etagam->3pigamma
    //if (phid == 5) {
    //}

    //cout << sig_type << endl;
    ALLCHAIN_GEN.Fill();

    /// cluster info
    TVector3 Pos_clust(0., 0., 0.);
    int promptnb = 0;
    double Emax_clust_tmp = 0., rt_clust_tmp = 0.;
    double Esum_clust_tmp = 0.;
    double E2_tmp = 0.;
    double pxmc_tmp = 0.;

    // test condition > 0
    int Pnum1_size = sizeof(Pnum1) / sizeof(Pnum1[0]);

    for (int i = 0; i < nclu; i ++) {// looping over clusters

      Pos_clust.SetXYZ(xcl[i], ycl[i], zcl[i]); // cluster position vector
      Angle_clust = Pos_clust.Theta() * TMath::RadToDeg(); // cluster polar angle
      Tof_clust = tcl[i] - Pos_clust.Mag() / speedc; // cluster tof
      Sigma2_T_clust = TMath::Power(tres * TMath::Sqrt(1000. / enecl[i]), 2) + TMath::Power(tcst, 2);
      Sigma_T_clust = TMath::Sqrt(Sigma2_T_clust);
      Esum_clust_tmp += enecl[i];
      
      // maximum cluster energy
      if (enecl[i] > Emax_clust_tmp) {
	Emax_clust_tmp = enecl[i];
      }

      // select prompt clusters
      if (charged_new[i] == 0 && intime[i] && enecl[i] > egammamin && Angle_clust > minangle && Angle_clust < (180. - minangle) && TMath::Abs(Tof_clust) < 3 * Sigma_T_clust && mmclu[i] != 5) {

	beta_clu = Pos_clust.Mag() / tcl[i] / speedc;
	tclusdiff = Tof_clust;

	// fill prompt cluster info: E, X, Y, Z and T
	E_list[promptnb] = enecl[i]; 
	X_list[promptnb] = xcl[i]; 
	Y_list[promptnb] = ycl[i];
	Z_list[promptnb] = zcl[i];
	T_list[promptnb] = tcl[i];

	pho_indx[promptnb] = Pnum1[i] - 1;
	clu_indx[promptnb] = i;
	pid_indx[promptnb] = Pid1[i];

	
	if (kineid != -999) {
	  //cout << "MC" << endl;
	  E_true_list[promptnb] = TMath::Sqrt(pxmc[Pnum1[i]-1] * pxmc[Pnum1[i]-1] + pymc[Pnum1[i]-1] * pymc[Pnum1[i]-1] + pzmc[Pnum1[i]-1] * pzmc[Pnum1[i]-1]);
	}
	//cout << E_true_list[promptnb] << endl;

	// sigma square
	rt_clust_tmp = TMath::Sqrt(TMath::Power(xcl[i], 2) + TMath::Power(ycl[i], 2));
	
	if (rt_clust_tmp > 200.) {// hits on barrel
	  Sigma2_X_list[promptnb] = xyres * xyres;
	  Sigma2_Y_list[promptnb] = xyres * xyres;
	  Sigma2_Z_list[promptnb] = TMath::Power(zres * TMath::Sqrt(1000. / enecl[i]), 2);
	}
	else {// hits on end-cup
	  Sigma2_X_list[promptnb] = xyres * xyres;
	  Sigma2_Y_list[promptnb] = TMath::Power(zres * TMath::Sqrt(1000. / enecl[i]), 2);
	  Sigma2_Z_list[promptnb] = xyres * xyres;
	}
	Sigma2_E_list[promptnb] = TMath::Power(eres * TMath::Sqrt(1000. * enecl[i]), 2);
	Sigma2_T_list[promptnb] = TMath::Power(tres * TMath::Sqrt(1000. / enecl[i]), 2) + TMath::Power(tcst, 2);

	promptnb ++;
      }// end prompt clusters section
      
      
    }// end looping over clusters
    Emax_clust = Emax_clust_tmp;
    Esum_clust = Esum_clust_tmp;
    //cout << "promptnb = " << promptnb << "\n";

    /// trck info
    TVectorD trkvect(3);
    trkvect = Getpiontrnb();
    trknb_prompt = trkvect(0), trkindx1 = trkvect(1), trkindx2 = trkvect(2);
    Bool_t ifbroken = IfBroken(trkindx1, trkindx2);

    Zv_trk1 = TMath::Abs(zv[trkindx1]);
    Zv_trk2 = TMath::Abs(zv[trkindx2]);
    //cout << zv[trkindx1] << endl;

    Rhov_trk1 = TMath::Sqrt(xv[trkindx1] * xv[trkindx1] + yv[trkindx1] * yv[trkindx1]);
    Rhov_trk2 = TMath::Sqrt(xv[trkindx2] * xv[trkindx2] + yv[trkindx2] * yv[trkindx2]);

    // CUT_COND 0, vertex cut
    if (Zv_trk1 > Zvmin || Zv_trk2 > Zvmin || Rhov_trk1 > Rhovmin || Rhov_trk2 > Rhovmin) continue;
    
    // CUT_COND 1, bit cut
    if(bit_select == 1) continue; 
    evnt_pre ++;

    //cout << "!!!!!!!!!!!!!!!!!!!!!!!!!" << promptnb << endl;
    hprompt_distr -> Fill(promptnb);
    
    // Event classification

    if (necls == 0) {// UFO
      evnt_unstr ++;
    }
    else {// streamed
      evnt_cls ++;
      for (int NrEC=0; NrEC<necls; NrEC++) {     

	//if (eclstream[NrEC] == 2) {
	hstr_distr -> Fill(eclstream[NrEC]);
	//}
	//cout<<eclstream[NrEC]<<endl;
	//if (eclstream[NrEC]==4) {

	//}
      }
    }


    // CUT_COND 2, ksl stream
    /*
    if (!IfStreamed(pstrnb)) continue;
    evtcls_indx = 1;
    evnt_str2 ++;
    */
    
    if (!IfStreamed(pstrnb)) {
      evtcls_indx = 1;
      evnt_str2 ++;
    }
    else {
      evtcls_indx = 999;
      evnt_no_str2 ++;
    }
    
    //if (evtcls_indx == 1) continue;
    //cout << evtcls_indx << endl;
    
    //cout << pstrnb << endl;

    ALLCHAIN_STR2.Fill();

    // CUT_COND 3, trigger
    //if (IfTriggered()) continue; // SYST. CHECK _EMC: pass EMC trigger,  _DC(): pass DC trigger, _EMCandDC(): pass EMC AND DC triggers     
    if (IfTriggered()) {
      trigger_indx = 1;
      evnt_trigger ++;
    }
    else {
      trigger_indx = 999;
      evnt_no_trigger ++;
    }
    
    // CUT4_COND4, FILFO
    //if (IfFilfoed()) continue; // SYST. CHECK, bit 20, standard FILFO cut
    if (IfFilfoed()) {
      filfo_indx = 1;
      evnt_filfo ++;
      //cout << filfo_indx << endl;
    }
    else {
      filfo_indx = 999;
      evnt_no_filfo ++;
    }
      
    // two tracks and 3 prompt photons
    //if (trknb_prompt != 2 || ifbroken ) continue; evnt_trk ++; 
    //if (promptnb != 3 ) continue; evnt_class ++;

    //if ((trknb_prompt != 2 || ifbroken ) || (promptnb != 3)) continue;
    if ((trknb_prompt != 2 || ifbroken ) || (promptnb != 3)) {
      fstate_indx = 1;
      evnt_fstate ++;
    }
    else {
      fstate_indx = 999;
      evnt_no_fstate ++;
    }

    if (fstate_indx == 999 && filfo_indx == 999 && trigger_indx == 999 && evtcls_indx == 999) evnt_cutted ++;
    
    /// 7C kinematical fit
    const int Row = 15, Col = 15, row = 7;
    const int nfloop_7C = 5;
	  
    TVectorD inputvect_7C(Row), sigma2vect_7C(Row), etakinfit_min_7C(Row);
    TVectorD sigma2vectorkinfit_min_7C(Row), pullkinfit(Row);

    for (int nr1 = 0; nr1 < promptnb - 2; nr1 ++) {// select 3 photon loop
      for (int nr2 = nr1 + 1; nr2 < promptnb - 1; nr2 ++) {
	for (int nr3 = nr2 + 1; nr3 < promptnb; nr3 ++) {
	  
	  int step = 0;

	  double lagravaluekinfittildeloop = 1e14, lagravaluekinfitloop_temp = 1e14, lagravaluekinfitloop = 1e14;
    
	  TVectorD inputvect(Row), sigma2vector(Row), sigma2vectorkinfitloop(Row);
	  TVectorD etatildekinfitloop(Row), etakinfitloop_temp(Row), etakinfitloop(Row), diffvectorkinfitloop(Row);
	  TVectorD gkinfitloop(row), lambdavectorloop(row), pullsvectorloop(Row);
  
	  TMatrixD Getakinfitloop(row, Col), Getakinfittransloop(Col, row), Vmatrixloop(Row, Col), Skinfitloop(row, row), InSkinfitloop(row, row), Vkinfitloop(Row, Col);

	  inputvect = FillInputvector_7C(Row, nr1, nr2, nr3, trkindx1, trkindx2);
	  //cout << "!!!! inspect inputvector before kin. fit. " << endl;
	  //cout << "nr1 = " << nr1 << ", nr2 = " << nr2 << ", nr3 = " << nr3 << endl;
	  //inputvect.Print();
	  sigma2vector = FillSigma2vector_7C(Row, nr1, nr2, nr3, trkindx1, trkindx2); //sigma2vector.Print();

	  // initialization
	  etakinfitloop_temp = inputvect, etakinfitloop = inputvect;
	  Vmatrixloop = CovMatrix(sigma2vector, Row, Col);
  
	  while (step < nfloop_7C) {
	    etatildekinfitloop = etakinfitloop_temp; //etatildekinfitloop.Print();
	    etakinfitloop_temp = etakinfitloop;

	    diffvectorkinfitloop = inputvect - etatildekinfitloop; //diffvectorkinfitloop.Print();
	    gkinfitloop = Gfunc_7C(etatildekinfitloop, row, Beam, trkindx1, trkindx2); //gkinfitloop.Print();
	    
	    lagravaluekinfittildeloop=lagravaluekinfitloop_temp;
	    lagravaluekinfitloop_temp=lagravaluekinfitloop;
	    
	    Getakinfitloop = Getafunc_7C(etatildekinfitloop, row, Col); //Getakinfitloop.Print();
	    Getakinfittransloop = Mtrans(Getakinfitloop);
	    Skinfitloop = Getakinfitloop * (Vmatrixloop * Getakinfittransloop); //Skinfitloop.Print();
	    InSkinfitloop = MInvert(Skinfitloop);
	    lambdavectorloop = Lambdavector(row, InSkinfitloop, Getakinfitloop, gkinfitloop, diffvectorkinfitloop); //lambdavectorloop.Print();

	    lagravaluekinfitloop = lambdavectorloop * (Getakinfitloop * diffvectorkinfitloop + gkinfitloop); //cout << lagravaluekinfitloop << endl;
	    
	    lagravaluekinfitloop_temp = lagravaluekinfitloop;
	    lagravaluekinfitloop = lagravaluekinfittildeloop;
	    
	    etakinfitloop = etakinfitfunc(Vmatrixloop, Getakinfittransloop, lambdavectorloop, inputvect, Row);
	    etakinfitloop_temp = etakinfitloop; //etakinfitloop_temp.Print();
	    etakinfitloop = etatildekinfitloop;
	    
	    // statistics
	    Vkinfitloop = Vmatrixloop - Vmatrixloop * (Getakinfittransloop * InSkinfitloop * Getakinfitloop) * Vmatrixloop; //Vkinfitloop.Print();
	    sigma2vectorkinfitloop = Fillsigma2vectorkinfit(Row, Vkinfitloop, Col, Row); //sigma2vectorkinfitloop.Print();
	    pullsvectorloop = Fillpullsvector(Row, sigma2vector, sigma2vectorkinfitloop, inputvect, etakinfitloop_temp); //pullsvectorloop.Print();
	    
	    step ++;
	  }

	  lagvalue_min_7C = lagravaluekinfitloop_temp;
	  inputvect_7C = inputvect;
	  sigma2vect_7C = sigma2vector;
	  etakinfit_min_7C = etakinfitloop_temp; 
	  pvalue=TMath::Prob(lagvalue_min_7C, 7);
	  sigma2vectorkinfit_min_7C = sigma2vectorkinfitloop;
	  pullkinfit = pullsvectorloop;
	  
	  
	  //cout << "lagvalue_min_7C = " <<  lagvalue_min_7C << "\n";
	  //etakinfit_min_7C.Print();
	}

      }

    }

    TLorentzVector TLVector_pmi = Gettrack4vectorkinfit(trkindx1);
    TLorentzVector TLVector_ppl = Gettrack4vectorkinfit(trkindx2);

    // select pi0 photons and fill all photon final states
    TLorentzVector pionphoton1_tmp, pionphoton2_tmp;
    //cout << "TLVector_pmi info:" << "\n"; TLVector_pmi.Print();
    
    TVectorD inputvect_permut(Row), inputvectfit_permut(Row);
    TVectorD sigma2vect_permut(Row);
    TVectorD inputvect_ordered(Row), inputvect_fitted_ordered(Row);
    
    double mgg_tmp = 0., mgg_tmp_checked = 0., mgg_diff_tmp = 0.;
    double error_mggsq = 0., chi2mgg_tmp = 0., chi2mgg_min = 1e14;
    double m3pi_tmp = 0., m3pi_diff_tmp = 0.;
    double error_m3pisq = 0., chi2m3pi_tmp = 0., chi2m3pi_min = 1e14;
    int nr1 = 0, nr2 = 0, nr3 = 0;
    int isrgam_indx = 0, pi0gam1_indx = 0, pi0gam2_indx = 0;
    const int error_type_indx = 0; // 1 for kin.fitted error

    for (int i = 0; i < 3; i ++) {// loop over number of photon candidates

      // determine permutation indices
      nr1 = i % 3; nr2 = (i + 1) % 3; nr3 = (i + 2) % 3;
      
      // input vectors
      inputvect_permut = Fillpermutvector(Row, inputvect_7C, nr1, nr2, nr3); //inputvect_permut.Print();
      sigma2vect_permut = Fillpermutvector(Row, sigma2vect_7C, nr1, nr2, nr3); //sigma2vect_permut.Print();

      //inputvect_permut = Fillpermutvector(Row, etakinfit_min_7C, nr1, nr2, nr3); //inputvect_permut.Print();
      //sigma2vect_permut = Fillpermutvector(Row, sigma2vectorkinfit_min_7C, nr1, nr2, nr3); //sigma2vect_permut.Print();

      // filling pi0 photon candidates in permutation
      pionphoton1_tmp = Getphoton4vector(inputvect_permut(5), inputvect_permut(6), inputvect_permut(7), inputvect_permut(8));
      pionphoton2_tmp = Getphoton4vector(inputvect_permut(10),inputvect_permut(11),inputvect_permut(12),inputvect_permut(13));
      //pionphoton1_tmp.Print();
      //pionphoton2_tmp.Print();

      // calculate mgg and chi2mgg
      mgg_tmp = TMath::Sqrt((pionphoton1_tmp + pionphoton2_tmp).M2());
      mgg_tmp_checked = TMath::Sqrt(2. * pionphoton1_tmp.E() * pionphoton2_tmp.E() * (1. - TMath::Cos(pionphoton1_tmp.Angle(pionphoton2_tmp.Vect()))));
      mgg_diff_tmp = mgg_tmp - massneupion;
      error_mggsq = TMath::Power(mgg_tmp / 2, 2) * (sigma2vect_permut(5) / TMath::Power(inputvect_permut(5), 2) + sigma2vect_permut(10) / TMath::Power(inputvect_permut(10), 2));
      //error_mggsq = 2.855 * 2.855;
      chi2mgg_tmp = TMath::Power(mgg_diff_tmp, 2) / error_mggsq;
      
      m3pi_tmp = inputvect_permut(0);
      m3pi_diff_tmp = m3pi_tmp - 209.;
      //error_m3pisq = TMath::Power(m3pi_tmp / 2, 2) * (sigma2vect_permut(0)); //2.63 * 2.63;
      chi2m3pi_tmp = chi2mgg_tmp + TMath::Power(m3pi_diff_tmp, 2.) / sigma2vect_permut(0);
      
      //cout << "(nr1, nr2, nr3) = (" << nr1 << ", " << nr2 << ", " << nr3 << "), mgg squared = " << mgg_tmp << "+/-" << error_mggsq << ", checked = " << mgg_tmp_checked << ", mgg_diff = " << mgg_diff_tmp << ", chi2mgg_tmp = " << chi2mgg_tmp << ", error type = " << error_type_indx << endl;
      //cout << "(nr1, nr2, nr3) = (" << nr1 << ", " << nr2 << ", " << nr3 << "), m3pi squared = " << m3pi_tmp << "+/-" << error_m3pisq << ", m3pi_diff = " << m3pi_diff_tmp << ", chi2m3pi_tmp = " << chi2m3pi_tmp << endl;

      //cout << "(nr1, nr2, nr3) = (" << nr1 << ", " << nr2 << ", " << nr3 << "), chi2mgg_tmp = " << chi2mgg_tmp << endl;
      // fill kin.fit values
      inputvectfit_permut = Fillpermutvector(Row, etakinfit_min_7C, nr1, nr2, nr3);

      //if (chi2m3pi_tmp < chi2m3pi_min) {
      if (chi2mgg_tmp < chi2mgg_min) {
	chi2mgg_min = chi2mgg_tmp;
	inputvect_ordered = inputvect_permut;
	inputvect_fitted_ordered = inputvectfit_permut;
	isrgam_indx = nr1;
	pi0gam1_indx = nr2;
	pi0gam2_indx = nr3;
	//cout << "min (nr1, nr2, nr3) = (" << nr1 << ", " << nr2 << ", " << nr3 << "), chi2mgg = " << chi2mgg_tmp << "\n";
      }
	
    }// end loop

    chi2mgg_min_7C = chi2mgg_min;

    // without kin.fit
    TLorentzVector TLVector_isrpho = Getphoton4vector(inputvect_ordered(0),inputvect_ordered(1), inputvect_ordered(2), inputvect_ordered(3));
    TLorentzVector TLVector_pi0pho1 = Getphoton4vector(inputvect_ordered(5), inputvect_ordered(6), inputvect_ordered(7), inputvect_ordered(8));
    TLorentzVector TLVector_pi0pho2 = Getphoton4vector(inputvect_ordered(10), inputvect_ordered(11), inputvect_ordered(12), inputvect_ordered(13));
    TLorentzVector TLvector_isrpho_miss = Beam - (TLVector_pi0pho1 + TLVector_pi0pho2 + TLVector_ppl + TLVector_pmi);

    IMisrpho_miss = TLvector_isrpho_miss.M2();

    // masses
    MASSLIST[0] = (TLVector_pi0pho1 + TLVector_pi0pho2).M(); // pi0 mass no fit
    MASSLIST[1] = (TLVector_pi0pho1 + TLVector_pi0pho2 + TLVector_ppl + TLVector_pmi).M(); // 3pi mass no fit

    // fill 7C kin.fitted values
    //inputvect_fitted_ordered.Print();
    TLorentzVector TLVector_isrpho_kinfit7C = Getphoton4vector(inputvect_fitted_ordered(0),inputvect_fitted_ordered(1), inputvect_fitted_ordered(2), inputvect_fitted_ordered(3));
    TLorentzVector TLVector_pi0pho1_kinfit7C = Getphoton4vector(inputvect_fitted_ordered(5), inputvect_fitted_ordered(6), inputvect_fitted_ordered(7), inputvect_fitted_ordered(8));
    TLorentzVector TLVector_pi0pho2_kinfit7C = Getphoton4vector(inputvect_fitted_ordered(10), inputvect_fitted_ordered(11), inputvect_fitted_ordered(12), inputvect_fitted_ordered(13));
    TLorentzVector TLVector_pi0gg12_kinfit7C = TLVector_pi0pho1_kinfit7C + TLVector_pi0pho2_kinfit7C;
    TLorentzVector TLVector_3pi_kinfit7C = TLVector_pi0pho1_kinfit7C + TLVector_pi0pho2_kinfit7C + TLVector_ppl + TLVector_pmi;
    TLorentzVector TLVector_2pi_kinfit7C = TLVector_ppl + TLVector_pmi;

    IM3pi_7C = (TLVector_pi0pho1_kinfit7C + TLVector_pi0pho2_kinfit7C + TLVector_ppl + TLVector_pmi).M();
    IM_3pi_nofit = (TLVector_pi0pho1 + TLVector_pi0pho2 + TLVector_ppl + TLVector_pmi).M();

    // boost 3pi system
    TVector3 Boost3vector_3pi = -TLVector_3pi_kinfit7C.BoostVector();
    TVector3 Boost3vector_2pi = -TLVector_2pi_kinfit7C.BoostVector();
    
    TLorentzVector TLVector_ppl_boost = TLVector_ppl;
    TLorentzVector TLVector_pmi_boost = TLVector_pmi;
    TLVector_ppl_boost.Boost(Boost3vector_2pi);
    TLVector_pmi_boost.Boost(Boost3vector_2pi);
    
    angle_ppl_3piboost = TLVector_ppl_boost.Angle(TLVector_pi0gg12_kinfit7C.Vect())*TMath::RadToDeg(); // angle between pi+ in the rest 3pi system with repect to isr photon cadidate

    betapi0 = (TLVector_pi0gg12_kinfit7C.Vect()).Mag() / TLVector_pi0gg12_kinfit7C.E();

    // define 4-mom for pho1, pho2 and pho3
    TVector3 pho1_vect, pho2_vect, pho3_vect; 
    TLorentzVector TLVector_pho1, TLVector_pho2, TLVector_pho3;
    TLorentzVector TLVector_3pi12, TLVector_pho13, TLVector_pho23;

    // list of 4_mom of all three kin.fit photons
    TLorentzVector TLVector_pho_list[3] = {TLVector_isrpho_kinfit7C, TLVector_pi0pho1_kinfit7C, TLVector_pi0pho2_kinfit7C};
    
    double pho1E = 0., pho3E = 1E6, pho2E = 0.;

    for (int i = 0; i < 3; i ++) {// find min and max

      if (TLVector_pho_list[i].E() > pho1E) {

	pho1E = TLVector_pho_list[i].E();
	
	TLVector_pho1 = TLVector_pho_list[i];
	
      }

      if (TLVector_pho_list[i].E() < pho3E) {

	pho3E = TLVector_pho_list[i].E();
	
	TLVector_pho3 = TLVector_pho_list[i];
	
      }

      //cout << "phoE_list[" << i << "] = " << TLVector_pho_list[i].E() << ", pho1E = " << TLVector_pho1.E()  << ", pho3E = "<< TLVector_pho3.E() << "\n\n";

      
    }

    for (int i = 0; i < 3; i ++) {//find middle

      if (TLVector_pho_list[i].E() != pho1E && TLVector_pho_list[i].E() != pho3E ) {

	TLVector_pho2 = TLVector_pho_list[i];

	//cout << "pho2E = " << TLVector_pho2.E() << endl;
	
      }
      
      
    }

    IM3pi_pi12 = (TLVector_pho1 + TLVector_pho2 + TLVector_ppl + TLVector_pmi).M();
    IM3pi_pi13 = (TLVector_pho1 + TLVector_pho3 + TLVector_ppl + TLVector_pmi).M();
    IM3pi_pi23 = (TLVector_pho2 + TLVector_pho3 + TLVector_ppl + TLVector_pmi).M();

    /*
    cout << "pho1E = " << TLVector_pho1.E() << ", pho2E = " << TLVector_pho2.E() << ", pho3E = "<< TLVector_pho3.E() << "\n"
	 << "IM3pi_pi12 = " << IM3pi_pi12 << ", IM3pi_pi13 = " << IM3pi_pi13 << ", IM3pi_pi23 = " << IM3pi_pi23 << ", IM3pi_7C = " << IM3pi_7C << "\n\n";
    */

    TVector3 pi0gam1_vect, pi0gam2_vect; // selected pi0 photon, 3pi mc true four vectors
    TLorentzVector TLVector_pi0pho1_true, TLVector_pi0pho2_true, TLVector_3pi_true;

    if (kineid != -999) {
      pi0gam1_vect.SetXYZ(pxmc[pho_indx[pi0gam1_indx]], pymc[pho_indx[pi0gam1_indx]], pzmc[pho_indx[pi0gam1_indx]]);
      pi0gam2_vect.SetXYZ(pxmc[pho_indx[pi0gam2_indx]], pymc[pho_indx[pi0gam2_indx]], pzmc[pho_indx[pi0gam2_indx]]);
    }
    
    TLVector_pi0pho1_true = GetLorentzVector(pi0gam1_vect, 0.);
    TLVector_pi0pho2_true = GetLorentzVector(pi0gam2_vect, 0.);

    TLVector_3pi_true = piminusMC_TLvect + piplusMC_TLvect + TLVector_pi0pho1_true + TLVector_pi0pho2_true;
    IM3pi_true = TLVector_3pi_true.M();

    // energy
    ENERGYLIST[0] = TLVector_isrpho_kinfit7C.E();
    ENERGYLIST[1] = TLVector_pi0pho1_kinfit7C.E(); 
    ENERGYLIST[2] = DeltaE(TLVector_ppl, TLVector_pmi);//(TLVector_ppl_boost.Vect() + TLVector_pmi_boost.Vect()).Mag() - (Beam.M() - TMath::Sqrt(masschpion * masschpion + pplmomsqr_boost) - TMath::Sqrt(masschpion * masschpion + pmimomsqr_boost));
    ENERGYLIST[3] = TLVector_pi0pho2_kinfit7C.E(); 
    ENERGYLIST[4] = TLVector_isrpho.E();

    // photons 4 vectors
    // pi0 photon1
    MOM4PHO1[0] = TLVector_pi0pho1_kinfit7C.E();
    MOM4PHO1[1] = TLVector_pi0pho1_kinfit7C.X();
    MOM4PHO1[2] = TLVector_pi0pho1_kinfit7C.Y();
    MOM4PHO1[3] = TLVector_pi0pho1_kinfit7C.Z();

    // pi0 photon2
    MOM4PHO2[0] = TLVector_pi0pho2_kinfit7C.E();
    MOM4PHO2[1] = TLVector_pi0pho2_kinfit7C.X();
    MOM4PHO2[2] = TLVector_pi0pho2_kinfit7C.Y();
    MOM4PHO2[3] = TLVector_pi0pho2_kinfit7C.Z();

    // pi0 photon3
    MOM4PHO3[0] = TLVector_isrpho_kinfit7C.E();
    MOM4PHO3[1] = TLVector_isrpho_kinfit7C.X();
    MOM4PHO3[2] = TLVector_isrpho_kinfit7C.Y();
    MOM4PHO3[3] = TLVector_isrpho_kinfit7C.Z();
    
    // track 4 vectors
    MOM4TRKPLUS[0] = TLVector_ppl.E();
    MOM4TRKPLUS[1] = TLVector_ppl.X();
    MOM4TRKPLUS[2] = TLVector_ppl.Y();
    MOM4TRKPLUS[3] = TLVector_ppl.Z();

    MOM4TRKMINS[0] = TLVector_pmi.E();
    MOM4TRKMINS[1] = TLVector_pmi.X();
    MOM4TRKMINS[2] = TLVector_pmi.Y();
    MOM4TRKMINS[3] = TLVector_pmi.Z();

    // masses
    trkmass = Trkmass(TLVector_ppl, TLVector_pmi);
    ppIM = (TLVector_ppl + TLVector_pmi).M();
    //betapi0 = (TLVector_pi0gg12_kinfit7C.Vect()).Mag() / TLVector_pi0gg12_kinfit7C.E();
    ppIM_beta = ((TLVector_ppl + TLVector_pmi).Vect()).Mag() / (TLVector_ppl + TLVector_pmi).E();
    //(TLVector_pi0pho1_kinfit7C + TLVector_pi0pho2_kinfit7C + TLVector_ppl + TLVector_pmi).M();
    IM3pi_beta = ((TLVector_pi0pho1_kinfit7C + TLVector_pi0pho2_kinfit7C + TLVector_ppl + TLVector_pmi).Vect()).Mag() / (TLVector_pi0pho1_kinfit7C + TLVector_pi0pho2_kinfit7C + TLVector_ppl + TLVector_pmi).E();
    
    IM_pi0_7C = TLVector_pi0gg12_kinfit7C.M();
    IM_pi0_nofit = (TLVector_pi0pho1 + TLVector_pi0pho2).M();
    MASSLIST[2] = IM_pi0_7C; //TLVector_pi0gg12_kinfit7C.M(); // pi0 mass
    MASSLIST[3] = IM3pi_7C; //(TLVector_pi0pho1_kinfit7C + TLVector_pi0pho2_kinfit7C + TLVector_ppl + TLVector_pmi).M(); // 3pi mass 
    MASSLIST[4] = trkmass; //Trkmass(TLVector_ppl, TLVector_pmi); // track mass
    MASSLIST[5] = ppIM; //(TLVector_ppl + TLVector_pmi).M(); // pi+ pi- invaraint mass
    MASSLIST[6] = (TLVector_pi0pho1_kinfit7C + TLVector_isrpho_kinfit7C).M(); // pi0 mass, gamma 13
    MASSLIST[7] = (TLVector_pi0pho2_kinfit7C + TLVector_isrpho_kinfit7C).M(); // pi0 mass, gamma 23
    MASSLIST[8] = (TLVector_pi0pho1_kinfit7C + TLVector_isrpho_kinfit7C + TLVector_ppl + TLVector_pmi).M(); // 3pi mass, gamma13 
    MASSLIST[9] = (TLVector_pi0pho2_kinfit7C + TLVector_isrpho_kinfit7C + TLVector_ppl + TLVector_pmi).M(); // 3pi mass, gamma23 
    MASSLIST[10] = (TLVector_ppl + TLVector_pmi).M2() * 1e-6; //invariant mass square of pi+ pi- in GeV^{2}
    MASSLIST[11] = (TLVector_ppl + TLVector_pi0gg12_kinfit7C).M2() * 1e-6; //invariant mass square of pi+ pi0 in GeV^{2}
    MASSLIST[12] = IMisrpho_miss; // invariant mass of recoil mass against final hadron system 
    MASSLIST[13] = ppIM_beta; // ppIM beta 
    MASSLIST[14] = IM3pi_beta; // IM3pi_beta
    //cout << "isrpho missing mass = " << MASSLIST[12] << endl;
    
    // angle
    ANGLELIST[0] = TLVector_pi0pho1_kinfit7C.Angle(TLVector_pi0pho2_kinfit7C.Vect())*TMath::RadToDeg(); // gamma12
    ANGLELIST[1] = TLVector_pi0pho2_kinfit7C.Angle(TLVector_isrpho_kinfit7C.Vect())*TMath::RadToDeg(); // gamma23
    ANGLELIST[2] = TMath::Cos(TLVector_isrpho_kinfit7C.Theta()); // cos isr
    ANGLELIST[3] = TLVector_pi0pho1_kinfit7C.Angle(TLVector_isrpho_kinfit7C.Vect())*TMath::RadToDeg(); // gamma13
    ANGLELIST[4] = TLVector_pmi.Angle(TLVector_ppl.Vect())*TMath::RadToDeg(); // pi+ pi-
    ANGLELIST[5] = angle_ppl_3piboost;

    // resolutions
    RESOLIST[0] = TLVector_pi0gg12_kinfit7C.M() - massneupion; // Mpi
    RESOLIST[1] = IM_3pi; // M3pi

    // pull
    PULLIST[0] = pullkinfit(0);
    PULLIST[1] = pullkinfit(1);
    PULLIST[2] = pullkinfit(2);
    PULLIST[3] = pullkinfit(3);
    PULLIST[4] = pullkinfit(4);
    PULLIST[5] = pullkinfit(5);
    PULLIST[6] = pullkinfit(6);
    PULLIST[7] = pullkinfit(7);
    PULLIST[8] = pullkinfit(8);
    PULLIST[9] = pullkinfit(9);
    PULLIST[10] = pullkinfit(10);
    PULLIST[11] = pullkinfit(11);
    PULLIST[12] = pullkinfit(12);
    PULLIST[13] = pullkinfit(13);
    PULLIST[14] = pullkinfit(14);
    
    //cout << Row << pullkinfit(0) << endl;
    test_value = ENERGYLIST[2];//DeltaE(TLVector_ppl_boost, TLVector_pmi_boost);

    /// fill trees
    if (lagvalue_min_7C > 70.) continue; 
    evnt_final ++;

    //cout << "SELECTED PHOTONS\n";
    double clu_Emin = 1e13;
    int clu_indx_Emin = 0;
    int bkg_counter = 1;
    for (int i = 0; i < 3; i ++) {
      if (pho_indx[i] == -1) {
	bkg_counter = bkg_counter + 1;
      }
      
      if (E_list[i] < clu_Emin) {
	clu_Emin = E_list[i];
	clu_indx_Emin = pho_indx[i];
      }
    }

    //int clust_mult = 1; // cluster multiplicity
    if (pho_indx[0] == pho_indx[1] || pho_indx[0] == pho_indx[2] || pho_indx[1] == pho_indx[2]) {
      bkg_counter = bkg_counter + 1;
      //clust_mult ++;
    }

    //bkg_indx = clust_mult;
    bkg_indx = bkg_counter;

    double E_ntmc;
    TString pid_type = "";
    for (int i = 0; i < ntmc; i ++) {
      if (pidmc[i] == 1) {
	pid_type = "gamma";
      }
      else if (pidmc[i] == 2) {
	pid_type = "e+";
      }
      else if (pidmc[i] == 3) {
	pid_type = "e-";
      }
      else if (pidmc[i] == 5) {
	pid_type = "mu+";
      }
      else if (pidmc[i] == 6) {
	pid_type = "mu-";
      }
      else if (pidmc[i] == 7) {
	pid_type = "pi0";
      }
      else if (pidmc[i] == 8) {
	pid_type = "pi+";
      }
      else if (pidmc[i] == 9) {
	pid_type = "pi-";
      }
      else {
	pid_type = "unknown";
      }
      
      E_ntmc = TMath::Sqrt(pxmc[i] * pxmc[i] + pymc[i] * pymc[i] + pzmc[i] * pzmc[i]); 

      //cout << "ntmc=" << i << ", pidmc(pid_type)=" << pidmc[i] << "(" << pid_type << "), pxmc=" << pxmc[i] << ", pymc=" << pymc[i] << ", pzmc=" << pzmc[i] << ", E= " << E_ntmc << "\n";
    }

    //cout << "\nTRUE\n";
    int Pnum_sum = 0;
    for (int i = 0; i < nclumc; i ++) {
      //if (i < 3 && Pnum1[i] - 1 == -1) {
      //if (i == 0 && Pnum1[i] - 1 == -1) {
      if (Pnum1[i] - 1 == -1) {
	Pnum_sum ++;
	//cout << Pnum1[i] - 1 << endl;
      }
      if (Pid1[i] == 1) {
	pid_type = "gamma";
      }
      else {
	pid_type = "others";
      }
      
      //cout << "nclumc="<< i << ", Npar=" << Npar[i] << ", Pnum1=" << Pnum1[i] - 1 << ", Pid1(pid_type)=" << Pid1[i] << "(" << pid_type << "), pxmc=" << pxmc[Pnum1[i]-1] << ", pymc=" << pymc[Pnum1[i]-1] << ", pzmc=" << pzmc[Pnum1[i]-1] << ", E= " << TMath::Sqrt(pxmc[Pnum1[i]-1]*pxmc[Pnum1[i]-1]+pymc[Pnum1[i]-1]*pymc[Pnum1[i]-1]+pzmc[Pnum1[i]-1]*pzmc[Pnum1[i]-1]) << "\n";
      //Pnum2="<< Pnum2[i] - 1 << ", Pid2=" << Pid2[i] << ", Pnum3=" << Pnum3[i] - 1 << ", Pid3=" << Pid3[i] << "\n";
    }

    //cout << "\nRECON\n";
    TLorentzVector TVect_clust_test; //Getphoton4vector(double E, double x, double y, double z)
    for (int i = 0; i < nclu; i ++) {
      TVect_clust_test = Getphoton4vector(enecl[i], xcl[i] - bx, ycl[i] - by, zcl[i] - bz); //SetPxPyPzE(xcl[i]-bx, ycl[i]-by, zcl[i]-bz, enecl[i]);
      //cout << xcl[i] - bx << ", " << ycl[i] - by << ", " <<  zcl[i] - bz << ", " << TVect_clust_test.M() << endl;
      
      //cout << "nclu=" << i << ", px=" << TVect_clust_test.X() << ", py=" << TVect_clust_test.Y() << ", pz=" << TVect_clust_test.Z() << ", E=" << TVect_clust_test.E() << "\n";

    }

    //cout << "\nCheck pi0 photons" << endl;
    int recon_indx_tmp = 0;

    /*
    cout << "pi0gam1_indx = " << pi0gam1_indx << ", pi0gam2_indx = " << pi0gam2_indx << "\n"
	 << "E_true_list[pi0gam1_indx] = " << E_true_list[pi0gam1_indx] << ", (px, py, pz) = (" << pxmc[pho_indx[pi0gam1_indx]] << ", " << pymc[pho_indx[pi0gam1_indx]] << ", " << pzmc[pho_indx[pi0gam1_indx]] << "), pho_indx = " << pho_indx[pi0gam1_indx] << "\n"
	 << "E_true_list[pi0gam2_indx] = " << E_true_list[pi0gam2_indx] << ", (px, py, pz) = (" << pxmc[pho_indx[pi0gam2_indx]] << ", " << pymc[pho_indx[pi0gam2_indx]] << ", " << pzmc[pho_indx[pi0gam2_indx]] << "), pho_indx = " << pho_indx[pi0gam2_indx] << "\n";
    */

    PI0PHORESD[0] = TLVector_pi0pho1_true.E() - TLVector_pi0pho1.E();
    PI0PHORESD[1] = TLVector_pi0pho2_true.E() - TLVector_pi0pho2.E();
    PI0PHORESD[2] = TLVector_pi0pho1_true.E() - TLVector_pi0pho1_kinfit7C.E();
    PI0PHORESD[3] = TLVector_pi0pho2_true.E() - TLVector_pi0pho2_kinfit7C.E();

    Bool_t checked[2] = {kFALSE, kFALSE};

    for (int i = 0; i < 2; i ++) {// find the first match

      //cout << "EPI0GAM[" << i << "] = " << EPI0GAM[i] << ", EPI0NTMC[" << i << "] = "<< EPI0NTMC[i] << endl;

      //if (TMath::Abs(EPI0GAM[i] - E_true_list[pi0gam1_indx]) < 1e-3) {
      if (pho_indx[pi0gam1_indx] == EPI0NTMC[i]) {
	recon_indx_tmp ++;
	checked[i] = kTRUE;

	//cout << E_true_list[pi0gam1_indx] << ", pi0gam1_indx = " << pi0gam1_indx << ", recon1_indx = " << recon_indx_tmp << ", first mathched at EPI0NTMC=" << EPI0NTMC[i] << ", pho_indx[pi0gam1_indx]=" << pho_indx[pi0gam1_indx] << "\n";
      }
     
      /*
      //if (TMath::Abs(EPI0GAM[i] - E_true_list[pi0gam2_indx]) < 1e-3) {
      if (pho_indx[pi0gam2_indx] == EPI0NTMC[i]) {
	recon_indx_tmp ++;
	checked[i] = kTRUE;
	cout << E_true_list[pi0gam2_indx] << ", pi0gam2_indx = " << pi0gam2_indx << ", recon2_indx = " << recon_indx_tmp << ", mathched " << checked[i] << " at " << i << "\n";
      }
      */

    }

    for (int i = 0; i < 2; i ++) {// find the first match

      //cout << checked[i] << endl;
	
      if (pho_indx[pi0gam2_indx] == EPI0NTMC[i] && checked[i]==kFALSE) {
	recon_indx_tmp ++;
	checked[i] = kTRUE;

	//cout << E_true_list[pi0gam2_indx] << ", pi0gam2_indx = " << pi0gam2_indx << ", recon2_indx = " << recon_indx_tmp << ", first mathched at EPI0NTMC=" << EPI0NTMC[i] << ", pho_indx[pi0gam2_indx]=" << pho_indx[pi0gam2_indx] << "\n";
      }
      
    }
    
    recon_indx = recon_indx_tmp;
    //cout << "Matched recon pi0 photons = " << recon_indx << endl;

    if (recon_indx == 0) {
      evnt_recon_type0 ++;
    }
    else if (recon_indx == 1) {
      evnt_recon_type1 ++;
    }
    else {
      evnt_recon_type2 ++;
          
    }

    herror_type -> Fill(recon_indx, bkg_indx);

    ALLCHAIN_CUT.Fill();
    
  }// end entry loop

  // save trees
  ALLCHAIN_GEN.Write();
  ALLCHAIN_STR2.Write();
  ALLCHAIN_CUT.Write();
  
  // summary
  cout << "SUMMARY" << "\n"
       << "Total number of events: " << evnt_tot << " (evnt_sig = " << evnt_sig << ", evnt_bkg = " << evnt_bkg << ") \n"
       << "Trigger:" << " evnt_trigger (trigger_indx = 1) = " << evnt_trigger << ", evnt_no_trigger (trigger_indx = 999) = " << evnt_no_trigger << ", sum = " << evnt_no_trigger + evnt_trigger << "\n"
       << "EvntClass:" << " evnt str2 (evtcls_indx = 1) = " << evnt_str2 << ", evnt_no_str2 (evtcls_indx = 999) = " << evnt_no_str2 << ", sum = " << evnt_no_str2 + evnt_str2 << "\n"
       << "FILFO:" << " evnt_filfo (filfo_indx = 1) = " << evnt_filfo << ", evnt_no_filfo (filfo_indx = 999) = " << evnt_no_filfo << ", sum = " << evnt_no_filfo + evnt_filfo << "\n"
       << "2 tracks and 3 prompt photons: " << " evnt_fstate (fstate_indx = 1) = " << evnt_fstate << ", evnt_no_fstate (fstate_indx = 999) = " << evnt_no_fstate << ", sum = " << evnt_no_fstate + evnt_fstate << "\n"
       << "Number of events after all cuts: " << evnt_cutted << endl;

  cout << "PARMETERS" << "\n"
       << "egammamin = " << egammamin << "\n"
       << "Zvmin = " << Zvmin << "\n"
       << "Rhovmin = " << Rhovmin << "\n";
    
}// end Main()

