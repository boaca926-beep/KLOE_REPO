#include <TString.h>
#include <TChain.h>
#include <iostream>
#include <fstream>
#include <string>

#include "../header/tree.h"

using namespace std;

int Analys(TString inputFile, TString outputFile){

  cout << inputFile << endl;
  
  Int_t RUNNR=0;
  TChain *tree = new TChain("ETAPPG/h1");
  // reading list of ROOT files from the file list
  // files can be commented with '!'
  string line;
  ifstream filelist(inputFile); // = delete
  //  ifstream filelist("filelist.txt");
  if (filelist.is_open()) {
    while (!filelist.eof()) {
      if (getline(filelist, line, '\n'))
        if (line[0] != '!') {
          tree->Add(line.data());
	  //cout << "Adding file: " << line << " to the chain of files" << endl;
          //tree->SetBranchAddress("run_nr",&RUNNR); 
          /*for (Int_t irow=0;irow<tree->GetEntries();irow++) {
      	 	tree->GetEntry(irow);      			
      	 }*/
	  RUNNR ++;
	  //cout<<RUNNR<<endl;
        }
    }
    filelist.close();
  } else {
    cout << "Unable to open filelist" << endl;
    return 0;
  }

  cout << inputFile << "\n"
       << "Submitted runs: " << RUNNR << endl;

  TFile *myfile;
  myfile = new TFile(outputFile,"RECREATE");
  
  MyClass *analysis = new MyClass(tree);
  analysis->Main();

  myfile->Write();
  myfile->Close();

  delete analysis;

  return 0;
  
}
