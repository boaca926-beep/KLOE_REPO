#include "sm_para.h"
#include "method.h"
#include "plot.h"
#include "path.h"
#include "fit_func.h"

int fit_func() {

  getObj(f_in);


  double sigma_max = 1606.53; //nb
  double BB = GetBB_new(sigma_max, mass_omega); //sigma_max * 1e-9 * (mass_omega * mass_omega) / 12. / pi;

  cout << "mass_omega = " << mass_omega << ", Gam_omega = " << Gam_omega << ", BB = " << BB << endl;
  
  TF1 *fun_bw = new TF1("func_bw", Func_bw, 0., 2000., 3);
  fun_bw -> SetParameters(mass_omega, Gam_omega, BB);
  
  int counter = 0;
  const int size_limit = 10;
  const int bin_indx = 119;
  
  for (Int_t irow = 0; irow < TINPUT -> GetEntries(); irow++) {// begin for loop

    TINPUT -> GetEntry(irow);

    M3PI[counter] = TINPUT -> GetLeaf("Br_m3pi") -> GetValue(0);
    M3PI_ERR[counter] = 0.;
    
    EFFICY[counter] = TINPUT -> GetLeaf("Br_efficy") -> GetValue(0);
    EFFICY_ERR[counter] = TINPUT -> GetLeaf("Br_efficy_err") -> GetValue(0);
    ISRLUMI[counter] = TINPUT -> GetLeaf("Br_isrlumi_apprx") -> GetValue(0);
    CRX3PI_BW[counter] = fun_bw -> Eval(M3PI[counter]);

    //N3PI_PRE[counter] = ISRLUMI[counter] * CRX3PI_BW[counter] * EFFICY[counter];
    N3PI_PRE[counter] = ISRLUMI[counter] * CRX3PI_BW[counter];

    if (irow + 1 == bin_indx) {
      cout << irow + 1 << ", m3pi = " << M3PI[counter] << ", EFFICY [%] = " << EFFICY[counter] * 100. << " +/- " << EFFICY_ERR[counter] * 100. << ", ISRLUMI = " << ISRLUMI[counter] << ", CRX3PI_BW = " << CRX3PI_BW[counter] << ", N3PI_PRE = " << N3PI_PRE[counter] << "\n";
    
    }
    
    //if (irow > size_limit) break;
    counter ++;
    
  }

  // smearing matrix
  const int binsize = hsmearmatr -> ProjectionY() -> GetNbinsX();
  const double hmin = hsig_true -> GetXaxis() -> GetXmin();
  const double hmax = hsig_true -> GetXaxis() -> GetXmax();

  cout << binsize << endl;

  // efficiency combined smearing matrix
  
  for (int j = 1; j <= binsize; j ++ ) {

    double X_sum = 0., X_mean = 0.;
    double w = 0., wX_sum = 0., w_sum = 0.;
    double err2_sum = 0.;
    for (int k = 1; k <= binsize; k ++ ) {

      X_sum += hsmearmatr_efficy -> GetBinContent(k, j) * N3PI_PRE[k - 1];
      err2_sum = get_sigma2_y(hsmearmatr_efficy -> GetBinContent(k, j), N3PI_PRE[k - 1], hsmearmatr_efficy_err -> GetBinContent(k, j), 0.);
      //X_sum += hsmearmatr_efficy -> GetBinContent(k, j);
      if (err2_sum !=0){
	w = 1 / err2_sum;
      }
      else {
	w = 0.;
      }
      wX_sum += w * hsmearmatr_efficy -> GetBinContent(k, j) * N3PI_PRE[k - 1];
      w_sum += w;
    }

    X_mean = wX_sum / w_sum;
    //X_mean = X_sum /  binsize;
    //cout << X_sum << ", " << X_mean << ", " << X_sum / binsize << endl;
      
    double prob_tmp = 0., prob_err_tmp = 0.;
    double prob_efficy_tmp = 0., prob_efficy_err_tmp = 0.;
    double nb_sig = 0., nb_sig_err = 0.;
    double nb_sig_true = 0., nb_sig_true_err = 0.;
    double nb_sig_gen = 0., nb_sig_gen_err = 0.;
    double nb_sig_smeared = 0.; // use semarmatr
    double nb_sig_smeared_efficy = 0., nb_sig_smeared_efficy_err = 0.; // use smearmatr_efficy
    double nb_rec_tmp = 0.;
    double nb_rec_tmp_err2 = 0., nb_rec_tmp_err2_sum = 0.;
    double efficy = 0.;
    double sigma2 = 0., sigma2_sum = 0.;
    double var_X = 0., varsum_X = 0.;
    double s2 = 0.;
    for (int i = 1; i <= binsize; i ++ ) {

      prob_tmp = hsmearmatr -> GetBinContent(i, j);
      prob_err_tmp = hsmearmatr_err -> GetBinError(i, j);

      prob_efficy_tmp = hsmearmatr_efficy -> GetBinContent(i, j);
      prob_efficy_err_tmp = hsmearmatr_efficy_err -> GetBinContent(i, j);
      
      efficy = hefficy -> GetBinContent(i);
      
      nb_sig_true = hsig_true -> GetBinContent(i);
      nb_sig_true_err = hsig_true -> GetBinError(i);
	
      nb_sig_gen = hsig_gen -> GetBinContent(i);
      nb_sig_gen_err = hsig_gen -> GetBinError(i);
      
      nb_sig_smeared += prob_tmp * nb_sig_true;
      nb_sig_smeared_efficy += prob_efficy_tmp * nb_sig_gen;

      nb_rec_tmp += prob_efficy_tmp * N3PI_PRE[i - 1];//prob_tmp * N3PI_PRE[i - 1] * EFFICY[i - 1];
      //nb_rec_tmp_err2_sum = 0.;

      //var_X = (prob_efficy_tmp * N3PI_PRE[i - 1] - X_mean) * (prob_efficy_tmp * N3PI_PRE[i - 1] - X_mean) / binsize; 
      //cout << X_mean << endl;
      
      if (prob_efficy_tmp != 0.) {
	sigma2 = get_sigma2_y(prob_efficy_tmp, nb_sig_gen, prob_efficy_err_tmp, nb_sig_gen_err);
	nb_rec_tmp_err2 = get_sigma2_y(prob_efficy_tmp, N3PI_PRE[i - 1], prob_efficy_err_tmp, 0.);
	w = 1 / nb_rec_tmp_err2;
	var_X = w * (prob_efficy_tmp * N3PI_PRE[i - 1] - X_mean) * (prob_efficy_tmp * N3PI_PRE[i - 1] - X_mean);
      }
      else {
	sigma2 = 0.;
	nb_rec_tmp_err2 = 0.;
	w = 0.;
	var_X = 0.;
      }

      // sum error2
      sigma2_sum += sigma2;
      s2 += var_X;
      w_sum += w;
      //nb_rec_tmp_err2_sum += nb_rec_tmp_err2 + var_X;
      //nb_rec_tmp_err2_sum += nb_rec_tmp_err2 + s2 / w_sum; //+ s2 / (w_sum - 1);
      nb_rec_tmp_err2_sum += s2 / (w_sum - 1);
     
      if (j == bin_indx) {
	//cout << "(rec_indx, true_indx) = (" << j << ", " << i << "), nb_sig_true = " << nb_sig_true << ", prob_tmp = " << prob_tmp << "+/-"<< prob_err_tmp << ", nb_sig_gen = " << nb_sig_gen << ", prob_efficy_tmp = " << prob_efficy_tmp << "+/-" << prob_efficy_err_tmp << ", N3PI_PRE = " << N3PI_PRE[i - 1] << ", nb_rec_tmp_err2_sum = " <<  nb_rec_tmp_err2_sum << endl;
      }
      
    }

    nb_sig = hcorrmatrix_projY -> GetBinContent(j);
    nb_sig_err = hcorrmatrix_projY -> GetBinError(j);
    nb_sig_smeared_efficy_err = TMath::Sqrt(sigma2_sum);

    N3PI_SMEAR[j -1] = nb_rec_tmp;
    N3PI_SMEAR_ERR[j -1] = TMath::Sqrt(nb_rec_tmp_err2_sum);

    //if (j == bin_indx) {
      cout << "rec bin = " << j << ", m3pi = " << hsig_true -> GetBinCenter(j) << ", nb_sig = " << nb_sig << "+/-" << nb_sig_err << ", nb_sig_smeared = " << nb_sig_smeared << ", nb_sig_smeared_efficy = " << nb_sig_smeared_efficy << "+/-" << nb_sig_smeared_efficy_err <<  ", nb_sig_true = " << hsig_true -> GetBinContent(j) << ", nb_sig_gen = " << hsig_gen -> GetBinContent(j) << ", N3PI_PRE = " << N3PI_PRE[j - 1] << ", N3PI_SMEAR = " << N3PI_SMEAR[j -1] << "+/-" << N3PI_SMEAR_ERR[j -1] << "\n\n";
      //}
    
    //cout << "rec bin = " << j << ", m3pi = " << hsig_true -> GetBinCenter(j) << ", N3PI_PRE = " << N3PI_PRE[j - 1] << ", N3PI_SMEAR = " << N3PI_SMEAR[j -1] << "+/-" << N3PI_SMEAR_ERR[j -1] << "\n\n";
      
  }

  TGraphErrors *gf = new TGraphErrors(binsize, M3PI, N3PI_SMEAR, M3PI_ERR, N3PI_SMEAR_ERR);

  TCanvas *cv = new TCanvas("cv", "n3pi smear", 0, 0, 1000, 800);

  gf -> GetXaxis() -> SetRangeUser(750., 810.);
  
  gf -> Draw("AP");
  
  return 0;
  
}
