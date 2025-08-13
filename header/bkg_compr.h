const TString output_folder = "../../bkg_compr_Eprompt_max";
const TString var_nm = "Eprompt_max";
const TString unit = "[MeV]";
const TString var_symb = "E_{#gamma}^{max}";

const int binsize = 300;
const double var_min = 150;
const double var_max = 450;
const double sfw1d_isr3pi = 4.60022e-02;

void checkList(TList *list_tmp){

  TIter next(list_tmp);
  TObject* object = 0;
  int obj_indx = 0;
  while ((object = next()))
    {

      cout << "[" << obj_indx << "]: " << object -> GetName() << endl;
      obj_indx ++;
      
    }

}

