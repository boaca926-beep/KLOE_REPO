// width
const int list_size = 5;
const TString para_nm = "Gamma_omega"; 
const double value_nomi = 8.67507; // nominal value
const double error_stat = 0.120853; // statistical error
  
const TString CUT_NM[list_size] = {"presel",
				   "bkg",
				   "sfw2d",
				   "fit",
				   "lumi"
};

const double SYSTERR_plus[list_size] = {0.1137,
					0.0546038,
					0,
					0,
					0.0192809
};

const double SYSTERR_neg[list_size] = {-0.113141,
				       -0.0801576,
				       0,
				       -0.0855686,
				       0
};
  
