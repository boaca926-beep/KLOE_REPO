// BB
const int list_size = 5;
const TString para_nm = "BB";  
const double value_nomi = 5.07; // nominal value (efficiency corrected)
const double error_stat = 6.99985; // statistical error

const TString CUT_NM[list_size] = {"presel",
				   "bkg",
				   "sfw2d",
				   "fit",
				   "lumi"
};

const double SYSTERR_plus[list_size] = {0.0843402,
					0.0224535,
					0,
					0.0577114,
					0.0192809
};

const double SYSTERR_neg[list_size] = {-0.0696828,
				       -0.0224535,
				       0,
				       0,
				       -0.0191654
};
  
