#include "../header/path.h"

//int getchains(TString inputFile = "../../input_norm_TDATA/input/sig") {
int getchains(TString inputFile, TString outputPath) {

  //
  TChain *CUT_CHAIN = new TChain("ALLCHAIN_CUT");
  TChain *GEN_CHAIN = new TChain("ALLCHAIN_GEN");

  cout << "Merge chains of MC type: " << " from " << inputFile << ", output "<< outputPath << endl;

  TString file_path = "";

  // 
  TSystemDirectory dir(inputFile, inputFile);
  TList *files = dir.GetListOfFiles();

  if (files) {
    TSystemFile *file;
    TString fname;
    TIter next(files);

    int f_indx = 0;
    while ((file=(TSystemFile*)next())) {
      
      fname = file -> GetName();
      f_indx ++;

      
      if (fname.Contains(".root")) {

	file_path = inputFile + "/" + fname;
	
	cout << "Adding file: " << file_path << " to the chain of files" << endl;

	GEN_CHAIN -> Add(file_path);
	CUT_CHAIN -> Add(file_path);
	
      }

    }

  }

  TFile *output = new TFile(outputPath + data_type + ".root", "RECREATE");

  TTree *ALLCHAIN_CUT = CUT_CHAIN -> CloneTree();
  TTree *ALLCHAIN_GEN = GEN_CHAIN -> CloneTree();
  
  ALLCHAIN_CUT->Write();
  ALLCHAIN_GEN->Write();
  
  output->Close();
  
  return 0;
  
}
