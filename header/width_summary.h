// width
const int list_size = 13;
const TString para_nm = "Gamma_omega"; // 8.72995 +/- 0.112604
const double value_nomi = 8.72995; // nominal value
const double error_stat = 0.112604; // statistical error
  
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


const double SYSTERR_plus[list_size] = {0.0453582,
					0.0245415,
					0.,
					0.00617825,
					0.0580603,
					0.0285264,
					0.0484302,
					0.0296298,
					0.0977192,
					0.0,
					0.0270241,
					0.0340537,
					0.};

const double SYSTERR_neg[list_size] = {-0.136389,
				       -0.0245415,
				       -0.0289126,
				       -0.0382834,
				       -0.0581105,
				       -0.0285264,
				       -0.0484302,
				       -0.0296298,
				       -0.0977192,
				       -0.0,
				       -0.0270241,
				       -0.0340537,
				       -0.};
  
