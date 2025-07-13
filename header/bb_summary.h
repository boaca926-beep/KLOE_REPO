// BB
const int list_size = 13;
const TString para_nm = "BB"; // 6.38161 +/- 0.0630261 
const double value_nomi = 6.38161; // nominal value
const double error_stat = 0.0630261; // statistical error

const TString CUT_NM[list_size] = {"Minimum cluster energy",
				   "Chi-square",
				   "Missing photon energy",
				   "Opening angle",
				   "pi0 velocity",
				   "Normalization",
				   "Mass bin width",
				   "Fitting range",
				   "Cross section models",
				   "Integrated luminosity",
				   "Triggers",
				   "FILFO",
				   "Event classification"};


const double SYSTERR_plus[list_size] = {0.0353243,
					0.0231834,
					0.0161245,
					0.0192525,
					0.0,
					0.00617284,
					0.0292947,
					0.0126963,
					0.0230429,
					0.0192021,
					0.0124216,
					0.0195723,
					0.};

const double SYSTERR_neg[list_size] = {-0.0,
				       -0.0231834,
				       -0.0,
				       -0.00325477,
				       -0.0,
				       -0.00617284,
				       -0.0292947,
				       -0.0126963,
				       -0.0230429,
				       -0.0190873,
				       -0.0124216,
				       -0.0195723,
				       -0.};
  
