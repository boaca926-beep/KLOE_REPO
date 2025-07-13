// omega mass
const int list_size = 13;
const TString para_nm = "M_omega"; //  782.732 +/- 0.0376286
const double value_nomi = 782.732; // nominal value
const double error_stat = 0.0376286; // statistical error

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


const double SYSTERR_plus[list_size] = {0.0322383,
					0.0171958,
					0.00809502,
					0.00278623,
					0.0241274,
					0.00256418,
					0.00736121,
					0.00664431,
					0.,
					0.,
					0.00589199,
					0.0110924,
					0.};

const double SYSTERR_neg[list_size] = {-0.0207026,
				       -0.0171958,
				       -0.,
				       -0.,
				       -0.0162089,
				       -0.00256418,
				       -0.00736121,
				       -0.00664431,
				       -0.,
				       -0.,
				       -0.00589199,
				       -0.0110924,
				       -0.};
  
