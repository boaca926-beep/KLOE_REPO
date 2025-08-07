// width
const int list_size = 5;
const TString para_nm = "mass_omega";
const double value_nomi = 782.730; // nominal value (efficiency corrected)
const double error_stat = 0.0398315; // statistical error
  
const TString CUT_NM[list_size] = {"presel",
				   "bkg",
				   "sfw2d",
				   "fit",
				   "lumi"
};

const double SYSTERR_plus[list_size] = {0.0453377,
					0.022596,
					0,
					0,
					0
};

const double SYSTERR_neg[list_size] = {-0.0627435,
				       -0.0280213,
				       0,
				       0,
				       0
};
  
