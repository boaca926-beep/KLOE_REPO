//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Apr  2 13:02:19 2025 by ROOT version 6.32.10
// from TTree h1/etappg
// found on file: exp41902.root
//////////////////////////////////////////////////////////

#ifndef MyClass_h
#define MyClass_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class MyClass {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           run_nr;
   Int_t           ev_nr;
   Int_t           mcflag;
   Int_t           phid;
   Int_t           atyp[3];
   Int_t           btyp[3];
   Int_t           t3down;
   Int_t           t3flag;
   Int_t           ecltrgw;
   Int_t           eclfilfo;
   Int_t           necls;
   Int_t           eclword[8];   //[necls]
   Int_t           eclstream[8];   //[necls]
   Int_t           ecltagnum[8];   //[necls]
   Int_t           kineid;
   Int_t           bit_select;
   Int_t           trgtype;
   Int_t           bphi;
   Int_t           ephi;
   Int_t           wphi;
   Int_t           t1d;
   Int_t           t2d;
   Float_t         bx;
   Float_t         by;
   Float_t         bz;
   Float_t         bpx;
   Float_t         bpy;
   Float_t         bpz;
   Float_t         bene;
   Int_t           nv;
   Int_t           vtx[20];   //[nv]
   Float_t         xv[20];   //[nv]
   Float_t         yv[20];   //[nv]
   Float_t         zv[20];   //[nv]
   Int_t           ntv;
   Int_t           iv[30];   //[ntv]
   Int_t           trknumv[30];   //[ntv]
   Float_t         curv[30];   //[ntv]
   Float_t         phiv[30];   //[ntv]
   Float_t         cotv[30];   //[ntv]
   Float_t         pxtv[30];   //[ntv]
   Float_t         pytv[30];   //[ntv]
   Float_t         pztv[30];   //[ntv]
   Float_t         pmodv[30];   //[ntv]
   Float_t         lenvtx[30];   //[ntv]
   Float_t         vtxchi2[30];   //[ntv]
   Int_t           pidtv[30];   //[ntv]
   Int_t           nt;
   Int_t           trkind[100];   //[nt]
   Int_t           trkver[100];   //[nt]
   Float_t         cur[100];   //[nt]
   Float_t         phi[100];   //[nt]
   Float_t         cot[100];   //[nt]
   Float_t         pxt[100];   //[nt]
   Float_t         pyt[100];   //[nt]
   Float_t         pzt[100];   //[nt]
   Float_t         pmod[100];   //[nt]
   Float_t         xfirst[100];   //[nt]
   Float_t         yfirst[100];   //[nt]
   Float_t         zfirst[100];   //[nt]
   Float_t         length[100];   //[nt]
   Float_t         curla[100];   //[nt]
   Float_t         phila[100];   //[nt]
   Float_t         cotla[100];   //[nt]
   Float_t         pxtla[100];   //[nt]
   Float_t         pytla[100];   //[nt]
   Float_t         pztla[100];   //[nt]
   Float_t         pmodla[100];   //[nt]
   Float_t         xlast[100];   //[nt]
   Float_t         ylast[100];   //[nt]
   Float_t         zlast[100];   //[nt]
   Float_t         xpca[100];   //[nt]
   Float_t         ypca[100];   //[nt]
   Float_t         zpca[100];   //[nt]
   Float_t         curpca[100];   //[nt]
   Float_t         cotpca[100];   //[nt]
   Float_t         phipca[100];   //[nt]
   Int_t           nhit[100];   //[nt]
   Int_t           nfithit[100];   //[nt]
   Float_t         sigcurv[100];   //[nt]
   Float_t         sigcot[100];   //[nt]
   Float_t         sigphi[100];   //[nt]
   Float_t         xp_dcw[100][7][2];   //[nt]
   Float_t         xp_bp[100][7][2];   //[nt]
   Float_t         xp_ip[100][7][2];   //[nt]
   Int_t           n_seg[100][2];   //[nt]
   Float_t         len_seg[100][3][2];   //[nt]
   Float_t         dist_ip[100][2];   //[nt]
   Int_t           good_extrap[100][2];   //[nt]
   Int_t           mc_corr[100];   //[nt]
   Int_t           ntfmc;
   Int_t           ncontr[100];   //[ntfmc]
   Int_t           trkine1[100];   //[ntfmc]
   Int_t           trtype1[100];   //[ntfmc]
   Int_t           trhits1[100];   //[ntfmc]
   Int_t           trkine2[100];   //[ntfmc]
   Int_t           trtype2[100];   //[ntfmc]
   Int_t           trhits2[100];   //[ntfmc]
   Int_t           trkine3[100];   //[ntfmc]
   Int_t           trtype3[100];   //[ntfmc]
   Int_t           trhits3[100];   //[ntfmc]
   Int_t           nclu;
   Float_t         enecl[100];   //[nclu]
   Float_t         tcl[100];   //[nclu]
   Float_t         xcl[100];   //[nclu]
   Float_t         ycl[100];   //[nclu]
   Float_t         zcl[100];   //[nclu]
   Float_t         mmclu[100];   //[nclu]
   Int_t           charged_old[100];   //[nclu]
   Int_t           charged_new[100];   //[nclu]
   Float_t         ene_plane[100][5];   //[nclu]
   Int_t           intime[100];   //[nclu]
   Int_t           nclumc;
   Int_t           Npar[100];   //[nclumc]
   Int_t           Pnum1[100];   //[nclumc]
   Int_t           Pid1[100];   //[nclumc]
   Int_t           Pnum2[100];   //[nclumc]
   Int_t           Pid2[100];   //[nclumc]
   Int_t           Pnum3[100];   //[nclumc]
   Int_t           Pid3[100];   //[nclumc]
   Int_t           ntmc;
   Int_t           kine[50];   //[ntmc]
   Int_t           pidmc[50];   //[ntmc]
   Int_t           virmom[50];   //[ntmc]
   Int_t           indv[50];   //[ntmc]
   Float_t         pxmc[50];   //[ntmc]
   Float_t         pymc[50];   //[ntmc]
   Float_t         pzmc[50];   //[ntmc]
   Int_t           nvtxmc;
   Int_t           kinmom[30];   //[nvtxmc]
   Int_t           mother[30];   //[nvtxmc]
   Float_t         tofvmc[30];   //[nvtxmc]
   Float_t         xvmc[30];   //[nvtxmc]
   Float_t         yvmc[30];   //[nvtxmc]
   Float_t         zvmc[30];   //[nvtxmc]
   Float_t         trkvtxmc[30];   //[nvtxmc]
   Int_t           ntclo;
   Float_t         charge[100];   //[ntclo]
   Int_t           trk_quality[100];   //[ntclo]
   Int_t           extr_zone[100];   //[ntclo]
   Float_t         extrapolation[100][9];   //[ntclo]
   Int_t           best_clu[100];   //[ntclo]
   Int_t           associated[100];   //[ntclo]
   Int_t           nassclu[100];   //[ntclo]
   Float_t         trk_distalong[100];   //[ntclo]
   Float_t         trk_dista_tra[100];   //[ntclo]
   Int_t           n_parts[100];   //[ntclo]
   Float_t         leng_parts[100][3];   //[ntclo]
   Float_t         mome_parts[100][3];   //[ntclo]
   Float_t         mome_partsend[100][3];   //[ntclo]
   Float_t         beta_parts[100][3];   //[ntclo]

   // List of branches
   TBranch        *b_run_nr;   //!
   TBranch        *b_ev_nr;   //!
   TBranch        *b_mcflag;   //!
   TBranch        *b_phid;   //!
   TBranch        *b_atyp;   //!
   TBranch        *b_btyp;   //!
   TBranch        *b_t3down;   //!
   TBranch        *b_t3flag;   //!
   TBranch        *b_ecltrgw;   //!
   TBranch        *b_eclfilfo;   //!
   TBranch        *b_necls;   //!
   TBranch        *b_eclword;   //!
   TBranch        *b_eclstream;   //!
   TBranch        *b_ecltagnum;   //!
   TBranch        *b_kineid;   //!
   TBranch        *b_bit_select;   //!
   TBranch        *b_trgtype;   //!
   TBranch        *b_bphi;   //!
   TBranch        *b_ephi;   //!
   TBranch        *b_wphi;   //!
   TBranch        *b_t1d;   //!
   TBranch        *b_t2d;   //!
   TBranch        *b_bx;   //!
   TBranch        *b_by;   //!
   TBranch        *b_bz;   //!
   TBranch        *b_bpx;   //!
   TBranch        *b_bpy;   //!
   TBranch        *b_bpz;   //!
   TBranch        *b_bene;   //!
   TBranch        *b_nv;   //!
   TBranch        *b_vtx;   //!
   TBranch        *b_xv;   //!
   TBranch        *b_yv;   //!
   TBranch        *b_zv;   //!
   TBranch        *b_ntv;   //!
   TBranch        *b_iv;   //!
   TBranch        *b_trknumv;   //!
   TBranch        *b_curv;   //!
   TBranch        *b_phiv;   //!
   TBranch        *b_cotv;   //!
   TBranch        *b_pxtv;   //!
   TBranch        *b_pytv;   //!
   TBranch        *b_pztv;   //!
   TBranch        *b_pmodv;   //!
   TBranch        *b_lenvtx;   //!
   TBranch        *b_vtxchi2;   //!
   TBranch        *b_pidtv;   //!
   TBranch        *b_nt;   //!
   TBranch        *b_trkind;   //!
   TBranch        *b_trkver;   //!
   TBranch        *b_cur;   //!
   TBranch        *b_phi;   //!
   TBranch        *b_cot;   //!
   TBranch        *b_pxt;   //!
   TBranch        *b_pyt;   //!
   TBranch        *b_pzt;   //!
   TBranch        *b_pmod;   //!
   TBranch        *b_xfirst;   //!
   TBranch        *b_yfirst;   //!
   TBranch        *b_zfirst;   //!
   TBranch        *b_length;   //!
   TBranch        *b_curla;   //!
   TBranch        *b_phila;   //!
   TBranch        *b_cotla;   //!
   TBranch        *b_pxtla;   //!
   TBranch        *b_pytla;   //!
   TBranch        *b_pztla;   //!
   TBranch        *b_pmodla;   //!
   TBranch        *b_xlast;   //!
   TBranch        *b_ylast;   //!
   TBranch        *b_zlast;   //!
   TBranch        *b_xpca;   //!
   TBranch        *b_ypca;   //!
   TBranch        *b_zpca;   //!
   TBranch        *b_curpca;   //!
   TBranch        *b_cotpca;   //!
   TBranch        *b_phipca;   //!
   TBranch        *b_nhit;   //!
   TBranch        *b_nfithit;   //!
   TBranch        *b_sigcurv;   //!
   TBranch        *b_sigcot;   //!
   TBranch        *b_sigphi;   //!
   TBranch        *b_xp_dcw;   //!
   TBranch        *b_xp_bp;   //!
   TBranch        *b_xp_ip;   //!
   TBranch        *b_n_seg;   //!
   TBranch        *b_len_seg;   //!
   TBranch        *b_dist_ip;   //!
   TBranch        *b_good_extrap;   //!
   TBranch        *b_mc_corr;   //!
   TBranch        *b_ntfmc;   //!
   TBranch        *b_ncontr;   //!
   TBranch        *b_trkine1;   //!
   TBranch        *b_trtype1;   //!
   TBranch        *b_trhits1;   //!
   TBranch        *b_trkine2;   //!
   TBranch        *b_trtype2;   //!
   TBranch        *b_trhits2;   //!
   TBranch        *b_trkine3;   //!
   TBranch        *b_trtype3;   //!
   TBranch        *b_trhits3;   //!
   TBranch        *b_nclu;   //!
   TBranch        *b_enecl;   //!
   TBranch        *b_tcl;   //!
   TBranch        *b_xcl;   //!
   TBranch        *b_ycl;   //!
   TBranch        *b_zcl;   //!
   TBranch        *b_mmclu;   //!
   TBranch        *b_charged_old;   //!
   TBranch        *b_charged_new;   //!
   TBranch        *b_ene_plane;   //!
   TBranch        *b_intime;   //!
   TBranch        *b_nclumc;   //!
   TBranch        *b_Npar;   //!
   TBranch        *b_Pnum1;   //!
   TBranch        *b_Pid1;   //!
   TBranch        *b_Pnum2;   //!
   TBranch        *b_Pid2;   //!
   TBranch        *b_Pnum3;   //!
   TBranch        *b_Pid3;   //!
   TBranch        *b_ntmc;   //!
   TBranch        *b_kine;   //!
   TBranch        *b_pidmc;   //!
   TBranch        *b_virmom;   //!
   TBranch        *b_indv;   //!
   TBranch        *b_pxmc;   //!
   TBranch        *b_pymc;   //!
   TBranch        *b_pzmc;   //!
   TBranch        *b_nvtxmc;   //!
   TBranch        *b_kinmom;   //!
   TBranch        *b_mother;   //!
   TBranch        *b_tofvmc;   //!
   TBranch        *b_xvmc;   //!
   TBranch        *b_yvmc;   //!
   TBranch        *b_zvmc;   //!
   TBranch        *b_trkvtxmc;   //!
   TBranch        *b_ntclo;   //!
   TBranch        *b_charge;   //!
   TBranch        *b_trk_quality;   //!
   TBranch        *b_extr_zone;   //!
   TBranch        *b_extrapolation;   //!
   TBranch        *b_best_clu;   //!
   TBranch        *b_associated;   //!
   TBranch        *b_nassclu;   //!
   TBranch        *b_trk_distalong;   //!
   TBranch        *b_trk_dista_tra;   //!
   TBranch        *b_n_parts;   //!
   TBranch        *b_leng_parts;   //!
   TBranch        *b_mome_parts;   //!
   TBranch        *b_mome_partsend;   //!
   TBranch        *b_beta_parts;   //!

   int pstrnb;
   double speedc;
   double massomega;
   double masschpion;
   double massneupion;
   double egammamin;
   double Zv_trk1;
   double Zv_trk2;
   double Zvmin;
   double Rhov_trk1;
   double Rhov_trk2;
   double Rhovmin;   
   
   //for smearing the MC tracks to be more like data
   TRandom *generator;
   
   MyClass(TTree *tree=0);
   virtual ~MyClass();
   virtual TVectorD Getpiontrnb();
   virtual double   Trkmass(TLorentzVector bestppl, TLorentzVector bestpmi);
   virtual double   DeltaE(TLorentzVector bestppl, TLorentzVector bestpmi);
   virtual TLorentzVector GetLorentzVector(TVector3 vector, double mass);

   // kinematic fit 7C
   virtual TVectorD FillInputvector_7C(int size, int index1, int index2, int index3, int indtr1, int indtr2);
   virtual TVectorD FillSigma2vector_7C(int size, int index1, int index2, int index3, int indtr1, int indtr2);
   virtual TVectorD Gfunc_7C(TVectorD etavector, int rownb, TLorentzVector beam, int indtr1, int indtr2);
   virtual TVectorD Lambdavector(int size, TMatrixD inveSmatrix, TMatrixD dgmatrix, TVectorD constr, TVectorD diffv);
   virtual TVectorD etakinfitfunc(TMatrixD Vmatrix, TMatrixD dgmatrixTrans, TVectorD lambdavector, TVectorD etatilde, int rownb);
   virtual TVectorD Fillsigma2vectorkinfit(int size, TMatrixD matrix, int matrixcol, int matrixrow);
   virtual TVectorD Fillpullsvector(int size, TVectorD sigma2vector_old, TVectorD sigma2vector_new, TVectorD inputvector_old, TVectorD inputvector_new);
   virtual TVectorD Fillpermutvector(int size, TVectorD input, int index1, int index2, int index3);
   virtual TMatrixD CovMatrix(TVectorD vector, int row, int col);
   virtual TMatrixD Mtrans(TMatrixD ma);
   virtual TMatrixD MInvert(TMatrixD ma);
   virtual TMatrixD Getafunc_7C(TVectorD etavector, int rownb, int colnb);
   virtual TLorentzVector Gettrack4vectorkinfit(int index);
   virtual TLorentzVector Getphoton4vector(double E, double x, double y, double z);
  
   
   // cut 
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Main();
   virtual bool     Notify();
   virtual bool     IfBroken(int idx1, int idx2);
   virtual bool     IfTriggered();
   virtual bool     IfFilfoed();
   virtual bool     IfStreamed(int pstrnb);
   virtual void     Show(Long64_t entry = -1);

   int pho_indx[100], clu_indx[100], pid_indx[100];
   
   double E_list[100], E_true_list[100], X_list[100], Y_list[100], Z_list[100], T_list[100];
   double Sigma2_E_list[100], Sigma2_X_list[100], Sigma2_Y_list[100], Sigma2_Z_list[100], Sigma2_T_list[100];
   double etakinfit_min_7C[15];   

   TLorentzVector Beam;

  TH1D *hist1d, *hrecon_typeI, *hIM3pi, *hIM3pi_bkg_indx_rest, *hIM3pi_bkg_indx0, *hprompt_distr, *hstr_distr;
  TH2D *hsmearmatr_clust, *hsmearmatr_trk, *hsmearmatr_typeI, *hsmearmatr_typeII, *herror_type;
   
   /*
   bool IfTriggered_EMC();
   bool IfTriggered_DC();
   bool IfTriggered_EMCandDC();
   bool IfFilfoed();
   bool IfFilfoed_28();
   
   
   TVectorD GetTVectorD(int size);
   
   
   
   */
  
};

#endif

#ifdef MyClass_cxx
MyClass::MyClass(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("exp41902.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("exp41902.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("exp41902.root:/ETAPPG");
      dir->GetObject("h1",tree);

   }
   Init(tree);

   // initializing vaiables
   pstrnb = 2;
   generator = new TRandom();
   speedc = 29.9792458;
   massomega = 782.65;
   masschpion = 139.57; // MeV
   massneupion = 134.977; // MeV
   Zv_trk1 = 0;
   Zv_trk2 = 0;
   Rhov_trk1 = 0;
   Rhov_trk2 = 0;

   egammamin = 15;
   Rhovmin = 4;
   Zvmin = 10;
   
   for (int i = 0; i < 100; i ++) {
     E_true_list[i] = 0.;
     E_list[i] = 0.;
     X_list[i] = 0.; 
     Y_list[i] = 0.;
     Z_list[i] = 0.; 
     T_list[i] = 0.;

     Sigma2_E_list[i] = 0.;
     Sigma2_X_list[i] = 0.;
     Sigma2_Y_list[i] = 0.;
     Sigma2_Z_list[i] = 0.;
     Sigma2_T_list[i] = 0.;

     pho_indx[i] = 0;
     clu_indx[i] = 0;
     pid_indx[i] = 0;
     
   }

   for (int i = 0; i < 15; i ++) {
     etakinfit_min_7C[i] = 0.;
   }

   hprompt_distr = new TH1D("prompt_distr", "number of prompt photon distribution", 10, 0, 10);
   hstr_distr = new TH1D("str_distr", "stream distribution", 10, 0, 10);
   hist1d = new TH1D("hist1d", "hist1d", 200, 0., 1000.);
   hsmearmatr_clust = new TH2D("hsmearmatr_clust", "hsmearmatr_clust", 100, 400., 1000., 100, 400., 1000.);
   hsmearmatr_typeI = new TH2D("hsmearmatr_typeI", "IM3pi true v.s. recon box", 100, 400., 1000., 100, 400., 1000.);
   hsmearmatr_typeII = new TH2D("hsmearmatr_typeII", "IM3pi true v.s. recon typeII", 100, 400., 1000., 100, 400., 1000.);
   
   hsmearmatr_trk = new TH2D("hsmearmatr_trk", "hsmearmatr_trk", 200, 0., 1000., 200, 0., 1000.);
   hrecon_typeI = new TH1D("hrecon_typeI", "hrecon_typeI", 100, 400., 1000.);
   hIM3pi = new TH1D("hIM3pi", "hIM3pi", 100, 400., 1000.);
   
   hIM3pi_bkg_indx_rest = new TH1D("hIM3pi_bkg_indx_rest", "hIM3pi_bkg_indx_rest", 100, 400., 1000.);
   hIM3pi_bkg_indx0 = new TH1D("hIM3pi_bkg_indx0", "hIM3pi_bkg_indx0", 100, 0., 1000.);

   herror_type = new TH2D("herror_type", "herror_type", 3, 0, 3, 3, 0, 3);
   
}

MyClass::~MyClass()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t MyClass::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t MyClass::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void MyClass::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run_nr", &run_nr, &b_run_nr);
   fChain->SetBranchAddress("ev_nr", &ev_nr, &b_ev_nr);
   fChain->SetBranchAddress("mcflag", &mcflag, &b_mcflag);
   fChain->SetBranchAddress("phid", &phid, &b_phid);
   fChain->SetBranchAddress("atyp", atyp, &b_atyp);
   fChain->SetBranchAddress("btyp", btyp, &b_btyp);
   fChain->SetBranchAddress("t3down", &t3down, &b_t3down);
   fChain->SetBranchAddress("t3flag", &t3flag, &b_t3flag);
   fChain->SetBranchAddress("ecltrgw", &ecltrgw, &b_ecltrgw);
   fChain->SetBranchAddress("eclfilfo", &eclfilfo, &b_eclfilfo);
   fChain->SetBranchAddress("necls", &necls, &b_necls);
   fChain->SetBranchAddress("eclword", eclword, &b_eclword);
   fChain->SetBranchAddress("eclstream", eclstream, &b_eclstream);
   fChain->SetBranchAddress("ecltagnum", ecltagnum, &b_ecltagnum);
   fChain->SetBranchAddress("kineid", &kineid, &b_kineid);
   fChain->SetBranchAddress("bit_select", &bit_select, &b_bit_select);
   fChain->SetBranchAddress("trgtype", &trgtype, &b_trgtype);
   fChain->SetBranchAddress("bphi", &bphi, &b_bphi);
   fChain->SetBranchAddress("ephi", &ephi, &b_ephi);
   fChain->SetBranchAddress("wphi", &wphi, &b_wphi);
   fChain->SetBranchAddress("t1d", &t1d, &b_t1d);
   fChain->SetBranchAddress("t2d", &t2d, &b_t2d);
   fChain->SetBranchAddress("bx", &bx, &b_bx);
   fChain->SetBranchAddress("by", &by, &b_by);
   fChain->SetBranchAddress("bz", &bz, &b_bz);
   fChain->SetBranchAddress("bpx", &bpx, &b_bpx);
   fChain->SetBranchAddress("bpy", &bpy, &b_bpy);
   fChain->SetBranchAddress("bpz", &bpz, &b_bpz);
   fChain->SetBranchAddress("bene", &bene, &b_bene);
   fChain->SetBranchAddress("nv", &nv, &b_nv);
   fChain->SetBranchAddress("vtx", vtx, &b_vtx);
   fChain->SetBranchAddress("xv", xv, &b_xv);
   fChain->SetBranchAddress("yv", yv, &b_yv);
   fChain->SetBranchAddress("zv", zv, &b_zv);
   fChain->SetBranchAddress("ntv", &ntv, &b_ntv);
   fChain->SetBranchAddress("iv", iv, &b_iv);
   fChain->SetBranchAddress("trknumv", trknumv, &b_trknumv);
   fChain->SetBranchAddress("curv", curv, &b_curv);
   fChain->SetBranchAddress("phiv", phiv, &b_phiv);
   fChain->SetBranchAddress("cotv", cotv, &b_cotv);
   fChain->SetBranchAddress("pxtv", pxtv, &b_pxtv);
   fChain->SetBranchAddress("pytv", pytv, &b_pytv);
   fChain->SetBranchAddress("pztv", pztv, &b_pztv);
   fChain->SetBranchAddress("pmodv", pmodv, &b_pmodv);
   fChain->SetBranchAddress("lenvtx", lenvtx, &b_lenvtx);
   fChain->SetBranchAddress("vtxchi2", vtxchi2, &b_vtxchi2);
   fChain->SetBranchAddress("pidtv", pidtv, &b_pidtv);
   fChain->SetBranchAddress("nt", &nt, &b_nt);
   fChain->SetBranchAddress("trkind", trkind, &b_trkind);
   fChain->SetBranchAddress("trkver", trkver, &b_trkver);
   fChain->SetBranchAddress("cur", cur, &b_cur);
   fChain->SetBranchAddress("phi", phi, &b_phi);
   fChain->SetBranchAddress("cot", cot, &b_cot);
   fChain->SetBranchAddress("pxt", pxt, &b_pxt);
   fChain->SetBranchAddress("pyt", pyt, &b_pyt);
   fChain->SetBranchAddress("pzt", pzt, &b_pzt);
   fChain->SetBranchAddress("pmod", pmod, &b_pmod);
   fChain->SetBranchAddress("xfirst", xfirst, &b_xfirst);
   fChain->SetBranchAddress("yfirst", yfirst, &b_yfirst);
   fChain->SetBranchAddress("zfirst", zfirst, &b_zfirst);
   fChain->SetBranchAddress("length", length, &b_length);
   fChain->SetBranchAddress("curla", curla, &b_curla);
   fChain->SetBranchAddress("phila", phila, &b_phila);
   fChain->SetBranchAddress("cotla", cotla, &b_cotla);
   fChain->SetBranchAddress("pxtla", pxtla, &b_pxtla);
   fChain->SetBranchAddress("pytla", pytla, &b_pytla);
   fChain->SetBranchAddress("pztla", pztla, &b_pztla);
   fChain->SetBranchAddress("pmodla", pmodla, &b_pmodla);
   fChain->SetBranchAddress("xlast", xlast, &b_xlast);
   fChain->SetBranchAddress("ylast", ylast, &b_ylast);
   fChain->SetBranchAddress("zlast", zlast, &b_zlast);
   fChain->SetBranchAddress("xpca", xpca, &b_xpca);
   fChain->SetBranchAddress("ypca", ypca, &b_ypca);
   fChain->SetBranchAddress("zpca", zpca, &b_zpca);
   fChain->SetBranchAddress("curpca", curpca, &b_curpca);
   fChain->SetBranchAddress("cotpca", cotpca, &b_cotpca);
   fChain->SetBranchAddress("phipca", phipca, &b_phipca);
   fChain->SetBranchAddress("nhit", nhit, &b_nhit);
   fChain->SetBranchAddress("nfithit", nfithit, &b_nfithit);
   fChain->SetBranchAddress("sigcurv", sigcurv, &b_sigcurv);
   fChain->SetBranchAddress("sigcot", sigcot, &b_sigcot);
   fChain->SetBranchAddress("sigphi", sigphi, &b_sigphi);
   fChain->SetBranchAddress("xp_dcw", xp_dcw, &b_xp_dcw);
   fChain->SetBranchAddress("xp_bp", xp_bp, &b_xp_bp);
   fChain->SetBranchAddress("xp_ip", xp_ip, &b_xp_ip);
   fChain->SetBranchAddress("n_seg", n_seg, &b_n_seg);
   fChain->SetBranchAddress("len_seg", len_seg, &b_len_seg);
   fChain->SetBranchAddress("dist_ip", dist_ip, &b_dist_ip);
   fChain->SetBranchAddress("good_extrap", good_extrap, &b_good_extrap);
   fChain->SetBranchAddress("mc_corr", mc_corr, &b_mc_corr);
   fChain->SetBranchAddress("ntfmc", &ntfmc, &b_ntfmc);
   fChain->SetBranchAddress("ncontr", ncontr, &b_ncontr);
   fChain->SetBranchAddress("trkine1", trkine1, &b_trkine1);
   fChain->SetBranchAddress("trtype1", trtype1, &b_trtype1);
   fChain->SetBranchAddress("trhits1", trhits1, &b_trhits1);
   fChain->SetBranchAddress("trkine2", trkine2, &b_trkine2);
   fChain->SetBranchAddress("trtype2", trtype2, &b_trtype2);
   fChain->SetBranchAddress("trhits2", trhits2, &b_trhits2);
   fChain->SetBranchAddress("trkine3", trkine3, &b_trkine3);
   fChain->SetBranchAddress("trtype3", trtype3, &b_trtype3);
   fChain->SetBranchAddress("trhits3", trhits3, &b_trhits3);
   fChain->SetBranchAddress("nclu", &nclu, &b_nclu);
   fChain->SetBranchAddress("enecl", enecl, &b_enecl);
   fChain->SetBranchAddress("tcl", tcl, &b_tcl);
   fChain->SetBranchAddress("xcl", xcl, &b_xcl);
   fChain->SetBranchAddress("ycl", ycl, &b_ycl);
   fChain->SetBranchAddress("zcl", zcl, &b_zcl);
   fChain->SetBranchAddress("mmclu", mmclu, &b_mmclu);
   fChain->SetBranchAddress("charged_old", charged_old, &b_charged_old);
   fChain->SetBranchAddress("charged_new", charged_new, &b_charged_new);
   fChain->SetBranchAddress("ene_plane", ene_plane, &b_ene_plane);
   fChain->SetBranchAddress("intime", intime, &b_intime);
   fChain->SetBranchAddress("nclumc", &nclumc, &b_nclumc);
   fChain->SetBranchAddress("Npar", Npar, &b_Npar);
   fChain->SetBranchAddress("Pnum1", Pnum1, &b_Pnum1);
   fChain->SetBranchAddress("Pid1", Pid1, &b_Pid1);
   fChain->SetBranchAddress("Pnum2", Pnum2, &b_Pnum2);
   fChain->SetBranchAddress("Pid2", Pid2, &b_Pid2);
   fChain->SetBranchAddress("Pnum3", Pnum3, &b_Pnum3);
   fChain->SetBranchAddress("Pid3", Pid3, &b_Pid3);
   fChain->SetBranchAddress("ntmc", &ntmc, &b_ntmc);
   fChain->SetBranchAddress("kine", kine, &b_kine);
   fChain->SetBranchAddress("pidmc", pidmc, &b_pidmc);
   fChain->SetBranchAddress("virmom", virmom, &b_virmom);
   fChain->SetBranchAddress("indv", indv, &b_indv);
   fChain->SetBranchAddress("pxmc", pxmc, &b_pxmc);
   fChain->SetBranchAddress("pymc", pymc, &b_pymc);
   fChain->SetBranchAddress("pzmc", pzmc, &b_pzmc);
   fChain->SetBranchAddress("nvtxmc", &nvtxmc, &b_nvtxmc);
   fChain->SetBranchAddress("kinmom", kinmom, &b_kinmom);
   fChain->SetBranchAddress("mother", mother, &b_mother);
   fChain->SetBranchAddress("tofvmc", tofvmc, &b_tofvmc);
   fChain->SetBranchAddress("xvmc", xvmc, &b_xvmc);
   fChain->SetBranchAddress("yvmc", yvmc, &b_yvmc);
   fChain->SetBranchAddress("zvmc", zvmc, &b_zvmc);
   fChain->SetBranchAddress("trkvtxmc", trkvtxmc, &b_trkvtxmc);
   fChain->SetBranchAddress("ntclo", &ntclo, &b_ntclo);
   fChain->SetBranchAddress("charge", charge, &b_charge);
   fChain->SetBranchAddress("trk_quality", trk_quality, &b_trk_quality);
   fChain->SetBranchAddress("extr_zone", extr_zone, &b_extr_zone);
   fChain->SetBranchAddress("extrapolation", extrapolation, &b_extrapolation);
   fChain->SetBranchAddress("best_clu", best_clu, &b_best_clu);
   fChain->SetBranchAddress("associated", associated, &b_associated);
   fChain->SetBranchAddress("nassclu", nassclu, &b_nassclu);
   fChain->SetBranchAddress("trk_distalong", trk_distalong, &b_trk_distalong);
   fChain->SetBranchAddress("trk_dista_tra", trk_dista_tra, &b_trk_dista_tra);
   fChain->SetBranchAddress("n_parts", n_parts, &b_n_parts);
   fChain->SetBranchAddress("leng_parts", leng_parts, &b_leng_parts);
   fChain->SetBranchAddress("mome_parts", mome_parts, &b_mome_parts);
   fChain->SetBranchAddress("mome_partsend", mome_partsend, &b_mome_partsend);
   fChain->SetBranchAddress("beta_parts", beta_parts, &b_beta_parts);
   Notify();
}

bool MyClass::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return true;
}

void MyClass::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t MyClass::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}

TLorentzVector MyClass::GetLorentzVector(TVector3 vector, double mass) {

  TLorentzVector tvector(0.,0.,0.,0.);
  double E = TMath::Sqrt(TMath::Power(mass, 2) + vector.Mag2());
  tvector.SetPxPyPzE(vector(0), vector(1), vector(2), E);
  //cout << tvector.M() << endl;
  
  return tvector;

}

double MyClass::DeltaE(TLorentzVector bestppl, TLorentzVector bestpmi) {
  double deltaE = 0.;
  
  TVector3 Boost3vector = -Beam.BoostVector();
    
  TLorentzVector bestppl_boost = bestppl;
  bestppl_boost.Boost(Boost3vector);

  TLorentzVector bestpmi_boost = bestpmi;
  bestpmi_boost.Boost(Boost3vector);

  double pplmomsqr_boost = (bestppl_boost.Vect()).Mag2();
  double pmimomsqr_boost = (bestpmi_boost.Vect()).Mag2();

  deltaE = (bestppl_boost.Vect() + bestpmi_boost.Vect()).Mag() - (Beam.M() - TMath::Sqrt(masschpion * masschpion + pplmomsqr_boost) - TMath::Sqrt(masschpion * masschpion + pmimomsqr_boost));

  //cout << "pplmomsqr_boost = " << pplmomsqr_boost << "\n"; 
  
  return deltaE;
}

double MyClass::Trkmass(TLorentzVector bestppl, TLorentzVector bestpmi) {
  TVector3 ppl3mom, pmi3mom, p3mom_miss;
  double ppl3mom_mag = 0., pmi3mom_mag = 0., p3mom_mag_miss = 0.;
  double sqrtS = 0., mtrk2 = 0., mtrk = 0.;
  double sqr2_denom = 0., sqr2_nom = 0.;

  sqrtS = Beam.E();
  //cout << "Beam.E = " << Beam.E() << ", Beam.rho = " << Beam.Rho() << endl;
  //Beam.Print();
  //cout << "Beam trans = " << Beam.Px() << Beam. endl;

  ppl3mom = bestppl.Vect();
  ppl3mom_mag = ppl3mom.Mag();

  pmi3mom = bestpmi.Vect();
  pmi3mom_mag = pmi3mom.Mag();

  p3mom_miss = Beam.Vect() - ppl3mom - pmi3mom;
  p3mom_mag_miss = p3mom_miss.Mag();

  //cout << "pmi3mom_mag = " << pmi3mom_mag << ", bestpmi.rho = " << bestpmi.Rho() << endl;

  sqr2_denom = TMath::Power((ppl3mom_mag * ppl3mom_mag - pmi3mom_mag * pmi3mom_mag + (sqrtS - p3mom_mag_miss) * (sqrtS - p3mom_mag_miss)),2);
  sqr2_nom = (p3mom_mag_miss - sqrtS) * (p3mom_mag_miss - sqrtS);

  mtrk2 = sqr2_denom / (4 * sqr2_nom) - ppl3mom_mag * ppl3mom_mag;

  //if(mtrk2 < 0.){
  //mtrk2 = 0.;
    //cout << "mtrk2 = " << mtrk2 << endl;
  //}
  
  mtrk = TMath::Sqrt(mtrk2);

  //cout << "sqrt(s) = " << Beam.E() << endl;
  //bestppl.Print();
  //bestpmi.Print();
  //cout << "p+ = " << ppl3mom_mag << ", p- = " << pmi3mom_mag  << endl;
  //cout << "p3mom_mag_miss = " << p3mom_mag_miss << endl;
  //cout << "trkmass = " << sqr2_denom << ", sqr2_nom = " << sqr2_nom << endl;
  //cout << "TRKMASS2 = " << mtrk2 << ", TRKMASS = " << mtrk << endl;

  //cout << "p3mom_mag_miss = " << p3mom_miss.Mag() << endl;
  //cout << "p3mom_mag_miss check = " << TMath::Sqrt(p3mom_miss.Mag2()) << endl;

  return mtrk;
}

TVectorD MyClass::Getpiontrnb() {
  TVectorD indexvector(3);
  int svar_pos = 0, svar_neg = 0, index1 = 0, index2 = 0;
  Float_t distance1 = 1000000.0, distance2 = 100000.0;
  Float_t dist = 0.0;

  for (int k = 0; k < nt; k++) {// begin loop
    dist = sqrt((bx-xpca[k])*(bx-xpca[k]) + (by-ypca[k])*(by-ypca[k]) +(bz-zpca[k])*(bz-zpca[k]));
    if (curpca[k] < 0 && dist < distance1) {
      distance1 = dist;
      index1 = k;
      svar_neg=1;
    }
    else if (curpca[k] > 0 && dist < distance2) {
      distance2 = dist;
      index2 = k;
      svar_pos=1;
    }
  }// end loop

  indexvector(0) = svar_neg+svar_pos;
  indexvector(1) = index1;
  indexvector(2) = index2; //indexvector.Print();

  return indexvector;
}

bool MyClass::IfBroken(int idx1, int idx2) {
  if (TMath::Abs(cot[idx1] + cot[idx2]) < 0.1 && TMath::Abs((cur[idx1] + cur[idx2]) / cur[idx1]) < 0.2) {
    //cout << TMath::Abs(cot[idx1] + cot[idx2]) << endl;
    return kTRUE;
  }
  else return kFALSE;
}

bool MyClass::IfTriggered() {
  bool passtrigger=kFALSE;

  if (trgtype == 6 || trgtype == 2 || trgtype == 4) passtrigger = kTRUE;
  //if (trgtype == 6) passtrigger = kTRUE;

  if (passtrigger) return kFALSE;

  else return kTRUE;
}

bool MyClass::IfFilfoed() {
  bool passfilfo=kFALSE;

  if (((eclfilfo & ( 1 << 20 )) >> 20) == 1) passfilfo=kTRUE;

  if (passfilfo) return kFALSE;

  else return kTRUE;

}

bool MyClass::IfStreamed(int pstrnb) {
  Bool_t bTagged = kFALSE, chrad = kFALSE;

  //cout << pstrnb << endl;
  for (int NrEC = 0; NrEC < necls; NrEC ++){

    if (eclstream[NrEC] == pstrnb) {

      bTagged = kTRUE;

      if ((ecltagnum[NrEC] &  1 ) == 1) chrad = kTRUE;

    }

  }

  //if (bTagged) return kTRUE;
  if (bTagged) return kTRUE;
  else return kFALSE;
}

TVectorD MyClass::FillInputvector_7C(int size, int index1, int index2, int index3, int indtr1, int indtr2) {
  TVectorD vector(size);
  double inputarray[15] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
  double r1sq=0., r2sq=0., r3sq=0., pp1sq=0., pp2sq=0.;
  double r1=0., r2=0., r3=0., px1=0., px2=0., px3=0., py1=0., py2=0., py3=0., pz1=0., pz2=0., pz3=0.;
  double Epp1=0., Epp2=0., pp1=0., pp2=0., ppx1=0., ppx2=0., ppy1=0., ppy2=0., ppz1=0., ppz2=0., pp1tran=0., pp2tran=0.;
  double x1=0., x2=0., x3=0., y1=0., y2=0., y3=0., z1=0., z2=0., z3=0., E1=0., E2=0., E3=0., t1=0., t2=0., t3=0.;

  // calulate x1, x2, x3, y1, y2, y3, z1, z2, z3
  // calulate E1, E2, E3, t1, t2, t3
  // calulate px1, px2, px3
  x1=X_list[index1], x2=X_list[index2], x3=X_list[index3];
  y1=Y_list[index1], y2=Y_list[index2], y3=Y_list[index3];
  z1=Z_list[index1], z2=Z_list[index2], z3=Z_list[index3];

  E1=E_list[index1], E2=E_list[index2], E3=E_list[index3];
  t1=T_list[index1], t2=T_list[index2], t3=T_list[index3];

  // measurements array
  inputarray[0]=E1, inputarray[1]=x1, inputarray[2]=y1, inputarray[3]=z1, inputarray[4]=t1; //cout<<inputarray[4]<<endl;
  inputarray[5]=E2, inputarray[6]=x2, inputarray[7]=y2, inputarray[8]=z2, inputarray[9]=t2;
  inputarray[10]=E3, inputarray[11]=x3, inputarray[12]=y3, inputarray[13]=z3, inputarray[14]=t3;
  return vector.Use(size,inputarray);
}

TVectorD MyClass::FillSigma2vector_7C(int size, int index1, int index2, int index3, int indtr1, int indtr2) {
  TVectorD vector(size);
  double sigma2array[15]={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};

  // fill sigma2 array
  sigma2array[0]=Sigma2_E_list[index1], sigma2array[1]=Sigma2_X_list[index1], sigma2array[2]=Sigma2_Y_list[index1], sigma2array[3]=Sigma2_Z_list[index1], sigma2array[4]=Sigma2_T_list[index1];
  sigma2array[5]=Sigma2_E_list[index2], sigma2array[6]=Sigma2_X_list[index2], sigma2array[7]=Sigma2_Y_list[index2], sigma2array[8]=Sigma2_Z_list[index2], sigma2array[9]=Sigma2_T_list[index2];
  sigma2array[10]=Sigma2_E_list[index3], sigma2array[11]=Sigma2_X_list[index3], sigma2array[12]=Sigma2_Y_list[index3], sigma2array[13]=Sigma2_Z_list[index3], sigma2array[14]=Sigma2_T_list[index3];

  //sigma2array[15]=sigcurv[indtr1], sigma2array[16]=sigcot[indtr1], sigma2array[17]=sigphi[indtr1];
  //sigma2array[18]=sigcurv[indtr2], sigma2array[19]=sigcot[indtr2], sigma2array[20]=sigphi[indtr2];

  return vector.Use(size,sigma2array);
}

TMatrixD MyClass::CovMatrix(TVectorD vector, int row, int col) {
  TMatrixD covmatrix_temp(row,col);

  for (int i = 0; i < row; i++) {
    for (int j = 0; j < col; j++) {
      if (i==j) covmatrix_temp(i,j)=vector(i);
    }
  }
  return covmatrix_temp;
}

TMatrixD MyClass::Mtrans(TMatrixD ma) {
  TMatrixD ma_temp = ma;
  return ma_temp.T();
}

TMatrixD MyClass::MInvert(TMatrixD ma) {
  TMatrixD ma_temp = ma;
  return ma_temp.Invert();
}

TVectorD MyClass::Gfunc_7C(TVectorD etavector, int rownb, TLorentzVector beam, int indtr1, int indtr2) {
  TVectorD vector(rownb);
  //double Beam_E=benergy;
  //beam.Print();

  double r1sq=0., r2sq=0., r3sq=0., pp1sq=0., pp2sq=0.;
  double r1=0., r2=0., r3=0., px1=0., px2=0., px3=0., py1=0., py2=0., py3=0., pz1=0., pz2=0., pz3=0.;
  double Epp1=0., Epp2=0., pp1=0., pp2=0., ppx1=0., ppx2=0., ppy1=0., ppy2=0., ppz1=0., ppz2=0., pp1tran=0., pp2tran=0.;
  double curv1=0., cotan1=0., phi1=0, curv2=0., cotan2=0., phi2=0;
  double x1=0., x2=0., x3=0., y1=0., y2=0., y3=0., z1=0., z2=0., z3=0., E1=0., E2=0., E3=0., t1=0., t2=0., t3=0.;
  double garray[7] = {0.,0.,0.,0.,0.,0.,0.};
  double costheta = 0., massfactor = 0.;

  E1=etavector(0), x1=etavector(1), y1=etavector(2), z1=etavector(3), t1=etavector(4);
  E2=etavector(5), x2=etavector(6), y2=etavector(7), z2=etavector(8), t2=etavector(9);
  E3=etavector(10), x3=etavector(11), y3=etavector(12), z3=etavector(13), t3=etavector(14);

  // calculate squared r1sq, r2sq, r3sq and calculate r1, r2, r3
  r1sq=TMath::Power(x1,2.)+TMath::Power(y1,2.)+TMath::Power(z1,2.);
  r2sq=TMath::Power(x2,2.)+TMath::Power(y2,2.)+TMath::Power(z2,2.);
  r3sq=TMath::Power(x3,2.)+TMath::Power(y3,2.)+TMath::Power(z3,2.);

  r1=TMath::Sqrt(r1sq), r2=TMath::Sqrt(r2sq), r3=TMath::Sqrt(r3sq);

  //
  costheta = (x1*x2+y1*y2+z1*z2)/(r1*r2);
  //massfactor = (massneupion*massneupion)/(2*(1-costheta));
  //std::cout<<"E1 = "<<E1<<", E2 = "<<etavector(5)<<std::endl;
  //std::cout<<"costheta = "<<costheta<<", massfactor = "<<massfactor<<std::endl;
  //std::cout<<"constraint E2 = "<<massfactor/E1<<std::endl;

  //E2=massfactor/E1;


  px1=E1*x1/r1, px2=E2*x2/r2, px3=E3*x3/r3;
  py1=E1*y1/r1, py2=E2*y2/r2, py3=E3*y3/r3;
  pz1=E1*z1/r1, pz2=E2*z2/r2, pz3=E3*z3/r3;

  px1=E1*x1/r1, px2=E2*x2/r2, px3=E3*x3/r3;
  py1=E1*y1/r1, py2=E2*y2/r2, py3=E3*y3/r3;
  pz1=E1*z1/r1, pz2=E2*z2/r2, pz3=E3*z3/r3;

  // calculate squared ppx1sq, ppx2s, pp1tran, pp2tran, pp1, pp2, Epp1, Epp2
  // ppi1.SetX(ptransvSmeared * TMath::Cos(phipca[index]));
  //ppi1.SetY(ptransvSmeared * TMath::Sin(phipca[index]));
  //ppi1.SetZ(ptransvSmeared * cotpca[index]);

  curv1=curpca[indtr1], cotan1=cotpca[indtr1], phi1=phipca[indtr1];
  curv2=curpca[indtr2], cotan2=cotpca[indtr2], phi2=phipca[indtr2];
  //cout << "curv1 = " << curv1 << ", cotan1 = " << cotan1 << endl;
  //cout << "curv2 = " << curv2 << ", cotan2 = " << cotan2 << endl;
  //cout<<cotan2-cotpca[indtr2]<<endl;
  //cout<<phi2-phipca[indtr2]<<endl;
  pp1tran=1000.*1./TMath::Abs(curv1);
  pp2tran=1000.*1./TMath::Abs(curv2);

  ppx1=pp1tran*TMath::Cos(phi1);
  ppx2=pp2tran*TMath::Cos(phi2);

  ppy1=pp1tran*TMath::Sin(phi1);
  ppy2=pp2tran*TMath::Sin(phi2);

  ppz1=pp1tran*cotan1;
  ppz2=pp2tran*cotan2;

  pp1sq=TMath::Power(ppx1,2.)+TMath::Power(ppy1,2.)+TMath::Power(ppz1,2.);
  pp2sq=TMath::Power(ppx2,2.)+TMath::Power(ppy2,2.)+TMath::Power(ppz2,2.);

  Epp1=TMath::Sqrt(masschpion*masschpion+pp1sq);
  Epp2=TMath::Sqrt(masschpion*masschpion+pp2sq);

  pp1=TMath::Sqrt(pp1sq);
  pp2=TMath::Sqrt(pp2sq);

  TLorentzVector tvector1(0.,0.,0.,0.), tvector2(0.,0.,0.,0.);
  tvector1.SetPxPyPzE(px1,py1,pz1,E1);
  tvector2.SetPxPyPzE(px2,py2,pz2,E2);

  double openangle_loop=tvector1.Angle(tvector2.Vect());
  double cosangle_loop=TMath::Cos(openangle_loop);

  garray[0]=E1+E2+E3+Epp2+Epp1-beam.E();
  //cout<<"!!!!!!!!!!!!!!Beam E="<<Beam.X()<<"!!!!!!!!!!"<<endl;
  garray[1]=px1+px2+px3+ppx1+ppx2-beam.X();
  garray[2]=py1+py2+py3+ppy1+ppy2-beam.Y();
  garray[3]=pz1+pz2+pz3+ppz1+ppz2-beam.Z();
  garray[4]=r1-speedc*t1;
  garray[5]=r2-speedc*t2; //cout<<garray[5]<<endl;
  garray[6]=r3-speedc*t3;

  //cout << "E1 = " << E1 << ", E2 = " << E2 << ", E3 = " << E3 << endl;
  //cout << "Epp2 = " << Epp2 << ", Epp1 = " << Epp1 << ", beamE = " << beam.E() << endl;
  //cout << "g_new" << endl;
  //cout << "g1 = " << garray[0] << ", g2 = " << garray[1] << ", g3 = " << garray[2] << ", g4 = " << garray[3] << ", g5 = " << garray[4] << ", g6 = " << garray[5] << ", g7 = " << garray[6] << endl;
  //cout << "r1 = " << r1 << ", t1 = " << t1 << endl;

  vector.Use(rownb,garray);

  return vector;
}

TMatrixD MyClass::Getafunc_7C(TVectorD etavector, int rownb, int colnb) {
  TMatrixD matrix_temp(rownb, colnb);
  //cout << rownb << endl;
  double dg1deta[15] = {1.,0.,0.,0.,0.,1.,0.,0.,0.,0.,1.,0.,0.,0.,0.};
  double dg2deta[15] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
  double dg3deta[15] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
  double dg4deta[15] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
  double dg5deta[15] = {0.,0.,0.,0.,-speedc,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
  double dg6deta[15] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,-speedc,0.,0.,0.,0.,0.};
  double dg7deta[15] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,-speedc};
  double r1sq=0., r2sq=0., r3sq=0., pp1sq=0., pp2sq=0.;
  double r1=0., r2=0., r3=0., px1=0., px2=0., px3=0., py1=0., py2=0., py3=0., pz1=0., pz2=0., pz3=0.;
  double Epp1=0., Epp2=0., pp1=0., pp2=0., ppx1=0., ppx2=0., ppy1=0., ppy2=0., ppz1=0., ppz2=0., pp1tran=0., pp2tran=0.;
  double curv1=0., cotan1=0., phi1=0, curv2=0., cotan2=0., phi2=0;
  double x1=0., x2=0., x3=0., y1=0., y2=0., y3=0., z1=0., z2=0., z3=0., E1=0., E2=0., E3=0., t1=0., t2=0., t3=0.;
  TVectorD dg1(colnb), dg2(colnb), dg3(colnb), dg4(colnb), dg5(colnb), dg6(colnb), dg7(colnb);
  TVectorD dgarray[7]={dg1, dg2, dg3, dg4, dg5, dg6, dg7};
  double hatE1=0., hatE2=0., hatE3=0, hatx1=0., hatx2=0., hatx3=0., haty1=0., haty2=0., haty3=0., hatz1=0., hatz2=0., hatz3=0., costheta = 0., massfactor = 0.;


  // calulate x1, x2, x3, y1, y2, y3, z1, z2, z3
  // calulate E1, E2, E3, t1, t2, t3
  // calulate px1, px2, px3
  // calculate squared r1sq, r2sq, r3sq and calculate r1, r2, r3
  x1=etavector(1), x2=etavector(6), x3=etavector(11);
  y1=etavector(2), y2=etavector(7), y3=etavector(12);
  z1=etavector(3), z2=etavector(8), z3=etavector(13);

  t1=etavector(4), t2=etavector(9), t3=etavector(14);

  r1sq=TMath::Power(x1,2.)+TMath::Power(y1,2.)+TMath::Power(z1,2.);
  r2sq=TMath::Power(x2,2.)+TMath::Power(y2,2.)+TMath::Power(z2,2.);
  r3sq=TMath::Power(x3,2.)+TMath::Power(y3,2.)+TMath::Power(z3,2.);

  r1=TMath::Sqrt(r1sq), r2=TMath::Sqrt(r2sq), r3=TMath::Sqrt(r3sq);

  E1=etavector(0);
  E2=etavector(5);
  E3=etavector(10);
  //
  costheta = (x1*x2+y1*y2+z1*z2)/(r1*r2);
  massfactor = (massneupion*massneupion)/(2*(1-costheta));
  //std::cout<<"E1 = "<<E1<<", E2 = "<<etavector(5)<<std::endl;
  //std::cout<<"costheta = "<<costheta<<", massfactor = "<<massfactor<<std::endl;
  //std::cout<<"constraint E2 = "<<massfactor/E1<<std::endl;

  //E2=massfactor/E1;


  px1=E1*(x1)/r1, px2=E2*(x2)/r2, px3=E3*(x3)/r3;
  py1=E1*(y1)/r1, py2=E2*(y2)/r2, py3=E3*(y3)/r3;
  pz1=E1*(z1)/r1, pz2=E2*(z2)/r2, pz3=E3*(z3)/r3;


  // calculate more variables
  hatE1=E1/r1, hatE2=E2/r2, hatE3=E3/r3, hatx1=(x1)/r1, hatx2=(x2)/r2, hatx3=(x3)/r3, haty1=(y1)/r1, haty2=(y2)/r2, haty3=(y3)/r3, hatz1=(z1)/r1, hatz2=(z2)/r2, hatz3=(z3)/r3;


  // fill dgdeta array

  dg2deta[0]=hatx1, dg2deta[1]=hatE1*(TMath::Power(haty1,2)+TMath::Power(hatz1,2)), dg2deta[2]=-hatE1*hatx1*haty1, dg2deta[3]=-hatE1*hatz1*hatx1;
  dg2deta[5]=hatx2, dg2deta[6]=hatE2*(TMath::Power(haty2,2)+TMath::Power(hatz2,2)), dg2deta[7]=-hatE2*hatx2*haty2, dg2deta[8]=-hatE2*hatz2*hatx2;
  dg2deta[10]=hatx3, dg2deta[11]=hatE3*(TMath::Power(haty3,2)+TMath::Power(hatz3,2)), dg2deta[12]=-hatE3*hatx3*haty3, dg2deta[13]=-hatE3*hatz3*hatx3;


  dg3deta[0]=haty1, dg3deta[1]=dg2deta[2], dg3deta[2]=hatE1*(TMath::Power(hatx1,2)+TMath::Power(hatz1,2)), dg3deta[3]=-hatE1*hatz1*haty1;
  dg3deta[5]=haty2, dg3deta[6]=dg2deta[7], dg3deta[7]=hatE2*(TMath::Power(hatx2,2)+TMath::Power(hatz2,2)), dg3deta[8]=-hatE2*hatz2*haty2;
  dg3deta[10]=haty3, dg3deta[11]=dg2deta[12], dg3deta[12]=hatE3*(TMath::Power(hatx3,2)+TMath::Power(hatz3,2)), dg3deta[13]=-hatE3*hatz3*haty3;


  dg4deta[0]=hatz1, dg4deta[1]=dg2deta[3], dg4deta[2]=dg3deta[3], dg4deta[3]=hatE1*(TMath::Power(hatx1,2)+TMath::Power(haty1,2));
  dg4deta[5]=hatz2, dg4deta[6]=dg2deta[8], dg4deta[7]=dg3deta[8], dg4deta[8]=hatE2*(TMath::Power(hatx2,2)+TMath::Power(haty2,2));
  dg4deta[10]=hatz3, dg4deta[11]=dg2deta[13], dg4deta[12]=dg3deta[13], dg4deta[13]=hatE3*(TMath::Power(hatx3,2)+TMath::Power(haty3,2));


  dg5deta[1]=hatx1, dg5deta[2]=haty1, dg5deta[3]=hatz1;

  dg6deta[6]=hatx2, dg6deta[7]=haty2, dg6deta[8]=hatz2;

  dg7deta[11]=hatx3, dg7deta[12]=haty3, dg7deta[13]=hatz3;


  // fill dgarray
  dgarray[0]=dg1.Use(15,dg1deta);
  dgarray[1]=dg2.Use(15,dg2deta);
  dgarray[2]=dg3.Use(15,dg3deta);
  dgarray[3]=dg4.Use(15,dg4deta); //cout << "error" << endl;
  dgarray[4]=dg5.Use(15,dg5deta);
  dgarray[5]=dg6.Use(15,dg6deta);
  dgarray[6]=dg7.Use(15,dg7deta);

  // fill G matrix
  TMatrixDRow(matrix_temp,0)=dg1, TMatrixDRow(matrix_temp,1)=dg2, TMatrixDRow(matrix_temp,2)=dg3; TMatrixDRow(matrix_temp,3)=dg4, TMatrixDRow(matrix_temp,4)=dg5, TMatrixDRow(matrix_temp,5)=dg6, TMatrixDRow(matrix_temp,6)=dg7;


  return matrix_temp;
}

TVectorD MyClass::Lambdavector(int size, TMatrixD inveSmatrix, TMatrixD dgmatrix, TVectorD constr, TVectorD diffv) {
  TVectorD vector(size), rvector(size); //rvector.Print();
  rvector=constr+dgmatrix*diffv; //rvector.Print();
  vector=inveSmatrix*rvector;

  return vector;

}

TVectorD MyClass::etakinfitfunc(TMatrixD Vmatrix, TMatrixD dgmatrixTrans, TVectorD lambdavector, TVectorD etatilde, int rownb) {
  TVectorD vector(rownb);
  vector=etatilde-(Vmatrix*dgmatrixTrans)*lambdavector;
  return vector;
}

TVectorD MyClass::Fillsigma2vectorkinfit(int size, TMatrixD matrix, int matrixcol, int matrixrow) {
  TVectorD vector(size);

  for (int i = 0; i < matrixrow; i++) {
    for (int j = 0; j< matrixcol; j++) {
      if (i==j) vector[i]=matrix[i][j];
    }
  }

  return vector;
}

TVectorD MyClass::Fillpullsvector(int size, TVectorD sigma2vector_old, TVectorD sigma2vector_new, TVectorD inputvector_old, TVectorD inputvector_new) {
  double Array[21] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
  TVectorD vector(size); //vector.Print();

  for (int i=0; i<size; i++) {
    //if((sigma2vector_old-sigma2vector_new)(i)<=0) printf("<=0!!!=%lf\n",(sigma2vector_old-sigma2vector_new)(i));
    Array[i]=((inputvector_old-inputvector_new)(i))/TMath::Sqrt((sigma2vector_old-sigma2vector_new)(i));
    //if (TMath::Abs(vector(i)) < 0.001) cout<<vector(i)<<endl;
    //vector(i)=(inputvector_old-inputvector_new)(i);
  }
  vector.Use(size,Array);
  return vector;
}

TLorentzVector MyClass::Gettrack4vectorkinfit(int index) {

  TLorentzVector tvector_smeared(0.,0.,0.,0.);
  TLorentzVector tvector1(0.,0.,0.,0.);
  double E=0., psq=0., ptran=0., px=0., py=0., pz=0.;
  Double_t smearing_factor = 0.1501;
  Double_t scale_factor = 0.000367;
  double curv = curpca[index];
  double phi = phipca[index];
  double cotan = cotpca[index];
  
  ptran = 1000. * 1. / TMath::Abs(curv);
  Double_t ptransvSmeared = ptran*(1.+scale_factor);
  //ptransvSmeared = generator->Gaus(ptransvSmeared,(TMath::Sqrt(sigcurv[index])/1000.)*smearing_factor*ptransvSmeared*ptransvSmeared);
  ptransvSmeared = generator->Gaus(ptransvSmeared,(TMath::Sqrt(sigcurv[index]))*smearing_factor*ptransvSmeared);
  
  /*
  cout << "scale_factor = " << scale_factor << ", smearing_factor = " << smearing_factor << "\n"
       << "ptran = " << ptran << ", ptran*(1.+scale_factor) = " << ptran*(1.+scale_factor) << ", ptransvSmeared = " << ptransvSmeared << "\n"
       << "TMath::Sqrt(sigcurv[index])) = " << TMath::Sqrt(sigcurv[index]) << ", smearing = " << smearing_factor << "\n";
  */
  
  //px=ptran*TMath::Cos(phi); py=ptran*TMath::Sin(phi), pz=ptran*cotan;
  //psq=TMath::Power(ptran,2)+TMath::Power(pz,2);
  //E=TMath::Sqrt(masschpion*masschpion+psq);
  //cout<<TMath::Abs(curv)<<endl;
  //std::cout<<"E = "<<E<<endl;

  //tvector.SetPxPyPzE(px,py,pz,E);

  tvector1.SetX(ptran * TMath::Cos(phi));
  tvector1.SetY(ptran * TMath::Sin(phi));
  tvector1.SetZ(ptran * cotan);
  tvector1.SetE(sqrt(masschpion*masschpion + tvector1.Vect().Mag2()));

  tvector_smeared.SetX(ptransvSmeared * TMath::Cos(phi));
  tvector_smeared.SetY(ptransvSmeared * TMath::Sin(phi));
  tvector_smeared.SetZ(ptransvSmeared * cotan);
  tvector_smeared.SetE(sqrt(masschpion*masschpion + tvector_smeared.Vect().Mag2()));
  
  //std::cout << "check ptran = " << (tvector1.Vect()).Mag() << endl;
  //std::cout<<"===================="<<endl;
  //cout << "ptran = " << ptran << endl;

  /*
    cout << "no smearing (px, py, pz, E) = (" << tvector1.X() << "," << tvector1.Y() << ", " << tvector1.Z() << ", " << tvector1.E() << ")" << "\n"
    << "smeared (px, py, pz, E) = (" << tvector_smeared.X() << "," << tvector_smeared.Y() << ", " << tvector_smeared.Z() << ", " << tvector_smeared.E() << ")" << "\n\n";
  */    


  //return tvector1;
  return tvector_smeared;
}

TVectorD MyClass::Fillpermutvector(int size, TVectorD input, int index1, int index2, int index3) {
  TVectorD vector(size);
  TVectorD testvector(size);

  for (int i=0;i<15;i++) {
    testvector(i)=i;
  }
  //input.Print();
  //index1=1, index2=0, index3=2;
  //printf("%d,%d,%d \n",5*index1,5*index2,5*index3);
  //cout<<testvector(index1)<<endl;
  vector(0)=input(5*index1), vector(1)=input(5*index1+1), vector(2)=input(5*index1+2), vector(3)=input(5*index1+3), vector(4)=input(5*index1+4);
  vector(5)=input(5*index2), vector(6)=input(5*index2+1), vector(7)=input(5*index2+2), vector(8)=input(5*index2+3), vector(9)=input(5*index2+4);
  vector(10)=input(5*index3), vector(11)=input(5*index3+1), vector(12)=input(5*index3+2), vector(13)=input(5*index3+3), vector(14)=input(5*index3+4);


  //testvector.Print();
  //input.Print();
  //vector.Print();
  //printf("=========================\n");

  return vector;
}

TLorentzVector MyClass::Getphoton4vector(double E, double x, double y, double z) {
  //given a cluster index returns the 4-mom of a photon
  TVector3 gamma(x-bx,y-by,z-bz);
  Double_t scale1;
  scale1=E/gamma.Mag();
  TLorentzVector gamma4mom(scale1*gamma, E);
  //cout << gamma4mom.M() << endl;
  return gamma4mom;

}

#endif // #ifdef MyClass_cxx
