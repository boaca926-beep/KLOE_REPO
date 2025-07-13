#include "sm_para.h"
#include "path.h"
#include "sfw2d.txt"
#include "sfw1d.txt"
#include "method.h"
#include "result.h"
#include "cut_para.h"

int result(){

  gErrorIgnoreLevel = kInfo;

  cout << "Extract omega parameters ..." << endl;

  // data amd bkg histos
  getObj(f_cut);
  
  //checkList(HIM3pi_fit);
  checkList(HSIG);

  // data
  cout << "gsf = " << gsf << ", exp_type = " << exp_type << endl;
  hdata = (TH1D *) HIM3pi_fit -> FindObject("h1d_IM3pi_" + exp_type) -> Clone();
  //hdata -> Draw();
  hdata -> SetName("hdata");
  binsize = hdata -> GetNbinsX();
  hmin = hdata -> GetXaxis() -> GetXmin();
  hmax = hdata -> GetXaxis() -> GetXmax();

  // bkg MC
  eeg_hist = (TH1D *) HIM3pi_fit -> FindObject("h1d_IM3pi_TEEG") -> Clone();
  //eeg_hist -> Draw();
  omegapi_hist = (TH1D *) HIM3pi_fit -> FindObject("h1d_IM3pi_TOMEGAPI") -> Clone();
  ksl_hist = (TH1D *) HIM3pi_fit -> FindObject("h1d_IM3pi_TKSL") -> Clone();
  etagam_hist = (TH1D *) HIM3pi_fit -> FindObject("h1d_IM3pi_TETAGAM") -> Clone();
  hkpm = (TH1D *) HIM3pi_fit -> FindObject("h1d_IM3pi_TKPM") -> Clone();
  hrhopi = (TH1D *) HIM3pi_fit -> FindObject("h1d_IM3pi_TRHOPI") -> Clone();
  
  mcrest_hist = (TH1D *) HIM3pi_fit -> FindObject("h1d_IM3pi_TBKGREST") -> Clone(); // bkgsum
  mcrest_hist -> Add(hkpm, 1.);
  mcrest_hist -> Add(hrhopi, 1.);

  // efficiency
  get_efficy();
  
  // smearing matrix and smearing signal IM3pi MC true distr.
  TRandom *generator = new TRandom();
  rnd=new TRandom3();

  
  //cout << "Indx = " << type_indx << ", name = " << hdata -> GetName() << endl;
  cout << "nb_data = " << hdata -> Integral(1, binsize) << endl;

  hcorrmatrix = getCorrMatrix(binsize, hmin, hmax);
  hcorrmatrix -> SetName("hcorrmatrix");
  hcorrmatrix -> Sumw2();
  
  hsmearmatr = getSmearMatrix(hcorrmatrix, binsize, hmin, hmax);
  hsmearmatr -> SetName("hsmearmatr");
  hsmearmatr -> Sumw2();

  hsmearmatr_err = getSmearMatrixErr(hcorrmatrix, binsize, hmin, hmax);
  hsmearmatr_err -> SetName("hsmearmatr_err");
  hsmearmatr_err -> Sumw2();

  hsmearmatr_efficy = getSmearMatrixEfficy(hcorrmatrix, binsize, hmin, hmax);
  hsmearmatr_efficy -> SetName("hsmearmatr_efficy");
  hsmearmatr_efficy -> Sumw2();

  hsmearmatr_efficy_err = getSmearMatrixEfficy_Err(hcorrmatrix, binsize, hmin, hmax);
  hsmearmatr_efficy_err -> SetName("hsmearmatr_efficy_err");
  hsmearmatr_efficy_err -> Sumw2();
  
  getSmearSig(); // smeared signal

  // MC normalization
  MCNorm(); 
  scaleGSF(); 

  // calculate MC yields from sfw2d
  
  // fill the output tree
  binsize = hdata -> GetNbinsX();
  double xmin = hdata -> GetXaxis() -> GetXmin();
  double xmax = hdata -> GetXaxis() -> GetXmax();

  cout << "binszie = " << binsize << ", x range(" << xmin << ", " << xmax << ") MeV/c^2\n";
  
  double m3pi = 0., m3pi_lower = 0.,  m3pi_upper = 0.;
  double evnt_data = 0., evnt_data_err = 0.;
  double evnt_eeg = 0., evnt_eeg_err = 0.;
  double evnt_ksl = 0., evnt_ksl_err = 0.;
  double evnt_omegapi = 0., evnt_omegapi_err = 0.;
  double evnt_etagam = 0., evnt_etagam_err = 0.;
  double evnt_mcrest = 0., evnt_mcrest_err = 0.;
  double evnt_bkgsum = 0.;
  double nb_isr3pi_obs = 0., nb_isr3pi_obs_err = 0.;
  double efficy = 0., efficy_err = 0.;
  double W0_full = 0.;
  double isrlumi_apprx = 0.;

  TTree* TINPUT = new TTree("TINPUT", "recreate");
  TINPUT -> SetAutoSave(0);

  TINPUT -> Branch("Br_m3pi", &m3pi, "Br_m3pi/D");
  TINPUT -> Branch("Br_efficy", &efficy, "Br_efficy/D");
  TINPUT -> Branch("Br_efficy_err", &efficy_err, "Br_efficy_err/D");
  TINPUT -> Branch("Br_isrlumi_apprx", &isrlumi_apprx, "Br_isrlumi_apprx/D");
  TINPUT -> Branch("Br_nb_isr3pi_obs", &nb_isr3pi_obs, "Br_nb_isr3pi_obs/D");
  TINPUT -> Branch("Br_nb_isr3pi_obs_err", &nb_isr3pi_obs_err, "Br_nb_isr3pi_obs_err/D");

  //hefficy -> Draw();

  for (int i = 1; i <= binsize; i ++ ) {

    m3pi = hdata -> GetBinCenter(i);
    m3pi_lower = hdata -> GetBinLowEdge(i);
    m3pi_upper = hdata -> GetBinLowEdge(i + 1);
    
    // data
    evnt_data = hdata -> GetBinContent(i);
    evnt_data_err = hdata -> GetBinError(i);

    // eeg
    evnt_eeg = heeg -> GetBinContent(i);
    //evnt_eeg_err = GetBkgErr(nb_data_sum, eeg_hist -> GetBinContent(i), nb_eeg_sum, feeg, feeg_err);
    evnt_eeg_err = heeg -> GetBinError(i);
    
    // ksl
    evnt_ksl = hksl -> GetBinContent(i);
    //evnt_ksl_err = GetBkgErr(nb_data_sum, ksl_hist -> GetBinContent(i), nb_ksl_sum, fksl, fksl_err);
    evnt_ksl_err = hksl -> GetBinError(i);

    // omegapi
    evnt_omegapi = homegapi -> GetBinContent(i);
    //evnt_omegapi_err = GetBkgErr(nb_data_sum, omegapi_hist -> GetBinContent(i), nb_omegapi_sum, fomegapi, fomegapi_err);
    evnt_omegapi_err = homegapi -> GetBinError(i);

    // etagam
    evnt_etagam = hetagam -> GetBinContent(i);
    //evnt_etagam_err = GetBkgErr(nb_data_sum, etagam_hist -> GetBinContent(i), nb_etagam_sum, fetagam, fetagam_err);
    evnt_etagam_err = hetagam -> GetBinError(i);

    // mcrest
    evnt_mcrest = hmcrest -> GetBinContent(i);
    //evnt_mcrest_err = GetBkgErr(nb_data_sum, mcrest_hist -> GetBinContent(i), nb_mcrest_sum, fmcrest, fmcrest_err);
    evnt_mcrest_err = hmcrest -> GetBinError(i);

    // bkgsum
    evnt_bkgsum = evnt_eeg + evnt_ksl + evnt_omegapi + evnt_etagam + evnt_mcrest;

    // isr3pi_obs
    nb_isr3pi_obs = evnt_data - evnt_bkgsum; //h1d_IM3pi_T3PIOBS -> GetBinContent(i);
    nb_isr3pi_obs_err = TMath::Sqrt(evnt_data + evnt_eeg_err * evnt_eeg_err + evnt_ksl_err * evnt_ksl_err + evnt_omegapi_err * evnt_omegapi_err + evnt_etagam_err * evnt_etagam_err + evnt_mcrest_err * evnt_mcrest_err);

    // efficiency
    efficy = hefficy -> GetBinContent(i);
    efficy_err = hefficy -> GetBinError(i);

    // isr lumi
    W0_full = Get_W0_full(&m3pi);
    
    isrlumi_apprx = GetISRLumi_apprx(m3pi, m3pi_lower, m3pi_upper, W0_full);

    if (i == 119) {

      // summary
      
      cout << "bin indx = " << i << ", m3pi = " << m3pi << ", m3pi_lower = " << m3pi_lower << ", m3pi_upper = " << m3pi_upper << "\n"
	   << "evnt_eeg_err = " << evnt_eeg_err << ", nb_data_sum = " << nb_data_sum << ", event_eeg = " << eeg_hist -> GetBinContent(i) << ", nb_eeg_sum = " << nb_eeg_sum << ", feeg = " << feeg << "+/-" << feeg_err << "\n"
	   << "evnt_ksl_err = " << evnt_ksl_err << ", nb_data_sum = " << nb_data_sum << ", event_ksl = " << ksl_hist -> GetBinContent(i) << ", nb_ksl_sum = " << nb_ksl_sum << ", fksl = " << fksl << "+/-" << fksl_err << "\n"
	   << "evnt_omegapi_err = " << evnt_omegapi_err << ", nb_data_sum = " << nb_data_sum << ", event_omegapi = " << omegapi_hist -> GetBinContent(i) << ", nb_omegapi_sum = " << nb_omegapi_sum << ", fomegapi = " << fomegapi << "+/-" << fomegapi_err << "\n"
	   << "evnt_etagam_err = " << evnt_etagam_err << ", nb_data_sum = " << nb_data_sum << ", event_etagam = " << etagam_hist -> GetBinContent(i) << ", nb_etagam_sum = " << nb_etagam_sum << ", fetagam = " << fetagam << "+/-" << fetagam_err << "\n"
	   << "evnt_mcrest_err = " << evnt_mcrest_err << ", nb_data_sum = " << nb_data_sum << ", event_mcrest = " << mcrest_hist -> GetBinContent(i) << ", nb_mcrest_sum = " << nb_mcrest_sum << ", fmcrest = " << fmcrest << "+/-" << fmcrest_err << "\n\n"

      	   << "\tevnt_data = " << evnt_data << " +/- " << evnt_data_err << "\n"
	   << "\tevnt_eeg = " << evnt_eeg << " +/- " << evnt_eeg_err << "\n"
	   << "\tevnt_ksl = " << evnt_ksl << " +/- " << evnt_ksl_err << "\n"
	   << "\tevnt_omgapi = " << evnt_omegapi << " +/- " << evnt_omegapi_err << "\n"
	   << "\tevnt_etagam = " << evnt_etagam << " +/- " << evnt_etagam_err << "\n"
	   << "\tevnt_mcrest = " << evnt_mcrest << " +/- " << evnt_mcrest_err << "\n"
	   << "\tevnt_bkgsum = " << evnt_bkgsum << "\n"
	   << "\tnb_isr3pi_obs = " << nb_isr3pi_obs << "+/-" << nb_isr3pi_obs_err << "\n"
	   << "\tefficy = " << efficy << "+/-" << efficy_err << "\n"
	   << "\tLumi_tot " << Lumi_tot << ", gsf = " << gsf << ", Lumi_int = " << Lumi_int << "\n"
	   << "\tisrlumi_apprx = " << isrlumi_apprx << "\n";

      
    }    
    
    TINPUT -> Fill();

    //if (i > 10) break;
    
  }

  TFile *f_result = new TFile("result.root", "recreate");
  
  hcorrmatrix -> Write();
  hsmearmatr -> Write();
  hsmearmatr_err -> Write();
  hsmearmatr_efficy -> Write();
  hsmearmatr_efficy_err -> Write();
  hefficy -> Write();

  hsig -> Write("hsig"); // SIG
  hdata -> Write("hdata"); // DATA
  heeg -> Write("heeg"); // scale EEG
  homegapi -> Write("homegapi"); // scale OMEGAPI
  hetagam -> Write("hetagam"); // scale ETAGAM
  hksl -> Write("hksl"); // scale KSL 
  hmcrest -> Write("hmcrest"); // scale MCREST

  hsig_true -> Write();
  hsig_gen -> Write();
  
  TINPUT -> Write();
  
  f_result -> Close();

  return 0;
  
}
