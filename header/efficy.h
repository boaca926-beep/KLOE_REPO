const TString inputFile = "../../input_small_evtcls_TUFO/cut/tree_pre";
const TString outputCut = "../../efficy_evtcls/";
const TString treeType = "TUFO";

/*
const double mass_sigma_nb = 1;
const double IM3pi_min = 380;
const double IM3pi_max = 1020;
const double IM3pi_sigma = 2.65;
const int IM3pi_bin = TMath::Nint((IM3pi_max - IM3pi_min) / mass_sigma_nb / IM3pi_sigma);
*/

const int list_size = 2;
const TString CUT_TYPE[list_size] = {"sel", "evtcls"};

TH1D *H1DLIST[list_size];

TList *HIM3PI = new TList();

