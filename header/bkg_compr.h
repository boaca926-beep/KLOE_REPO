const TString output_folder = "../../bkg_compr_pull_E1";
const TString var_nm = "pull_x1";
const TString unit = "";
const TString var_symb = "pull x1";

const int binsize = 200;
const double var_min = -5;
const double var_max = 5;
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

