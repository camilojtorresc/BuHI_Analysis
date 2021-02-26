//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Dec 22 15:10:22 2020 by ROOT version 6.16/01
// from TTree ntuple/Bc ntuple
// found on file: dataMC/Rootuple_Bu_JpsiK_2016HI-MCAOD_bacht1.root
//////////////////////////////////////////////////////////

#ifndef classntupleMC_h
#define classntupleMC_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"
#include "vector"
#include "TLorentzVector.h"
#include "TVector3.h"

class classntupleMC {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   UInt_t          nB;
   UInt_t          nMu;
   vector<float>   *J_px1;
   vector<float>   *J_py1;
   vector<float>   *J_pz1;
   vector<float>   *J_px2;
   vector<float>   *J_py2;
   vector<float>   *J_pz2;
   vector<float>   *B_k_px;
   vector<float>   *B_k_py;
   vector<float>   *B_k_pz;
   vector<int>     *pi1Q;
   vector<int>     *pi1_trackerhits;
   vector<int>     *pi1_pixelhits;
   vector<float>   *pi1dz;
   vector<float>   *pi1dzE;
   vector<float>   *pi1dxy;
   vector<float>   *pi1dxyE;
   vector<float>   *B_charge;
   vector<float>   *B_mass;
   vector<float>   *B_px;
   vector<float>   *B_py;
   vector<float>   *B_pz;
   vector<int>     *B_k_charge1;
   vector<float>   *B_k_px_track;
   vector<float>   *B_k_py_track;
   vector<float>   *B_k_pz_track;
   vector<float>   *B_J_mass;
   vector<float>   *B_J_px;
   vector<float>   *B_J_py;
   vector<float>   *B_J_pz;
   vector<float>   *B_J_px1;
   vector<float>   *B_J_py1;
   vector<float>   *B_J_pz1;
   vector<int>     *B_J_charge1;
   vector<float>   *B_J_px2;
   vector<float>   *B_J_py2;
   vector<float>   *B_J_pz2;
   vector<int>     *B_J_charge2;
   vector<int>     *muon1Q;
   vector<int>     *muon2Q;
   vector<int>     *muon1Soft;
   vector<int>     *muon2Soft;
   vector<int>     *muon1StationTight;
   vector<int>     *muon2StationTight;
   vector<int>     *mu1_trackerhits;
   vector<int>     *mu2_trackerhits;
   vector<int>     *mu1_pixelhits;
   vector<int>     *mu2_pixelhits;
   vector<float>   *mu1dz;
   vector<float>   *mu1dzE;
   vector<float>   *mu1dxy;
   vector<float>   *mu1dxyE;
   vector<float>   *mu2dz;
   vector<float>   *mu2dzE;
   vector<float>   *mu2dxy;
   vector<float>   *mu2dxyE;
   vector<float>   *B_Prob;
   vector<float>   *B_J_Prob;
   vector<float>   *B_DecayVtxX;
   vector<float>   *B_DecayVtxY;
   vector<float>   *B_DecayVtxZ;
   vector<double>  *B_DecayVtxXE;
   vector<double>  *B_DecayVtxYE;
   vector<double>  *B_DecayVtxZE;
   vector<double>  *B_DecayVtxXYE;
   vector<double>  *B_DecayVtxXZE;
   vector<double>  *B_DecayVtxYZE;
   vector<float>   *B_J_DecayVtxX;
   vector<float>   *B_J_DecayVtxY;
   vector<float>   *B_J_DecayVtxZ;
   vector<float>   *B_J_DecayVtxXE;
   vector<float>   *B_J_DecayVtxYE;
   vector<float>   *B_J_DecayVtxZE;
   vector<float>   *B_J_DecayVtxXYE;
   vector<float>   *B_J_DecayVtxXZE;
   vector<float>   *B_J_DecayVtxYZE;
   Float_t         priVtxX;
   Float_t         priVtxY;
   Float_t         priVtxZ;
   Float_t         priVtxXE;
   Float_t         priVtxYE;
   Float_t         priVtxZE;
   Float_t         priVtxXYE;
   Float_t         priVtxXZE;
   Float_t         priVtxYZE;
   Float_t         priVtxCL;
   vector<float>   *pVtxIPX;
   vector<float>   *pVtxIPY;
   vector<float>   *pVtxIPZ;
   vector<float>   *pVtxIPXE;
   vector<float>   *pVtxIPYE;
   vector<float>   *pVtxIPZE;
   vector<float>   *pVtxIPXYE;
   vector<float>   *pVtxIPXZE;
   vector<float>   *pVtxIPYZE;
   vector<float>   *pVtxIPCL;
   vector<float>   *priRfVtxX;
   vector<float>   *priRfVtxY;
   vector<float>   *priRfVtxZ;
   vector<float>   *priRfVtxXE;
   vector<float>   *priRfVtxYE;
   vector<float>   *priRfVtxZE;
   vector<float>   *priRfVtxXYE;
   vector<float>   *priRfVtxXZE;
   vector<float>   *priRfVtxYZE;
   vector<float>   *priRfVtxCL;
   vector<int>     *priRfNTrkDif;
   UInt_t          nVtx;
   Int_t           nTrk;
   Int_t           nTrk_Vtx;
   vector<int>     *nTrk_VtxIP;
   vector<int>     *nTrk_VtxIP_Q;
   Int_t           run;
   Int_t           event;
   Int_t           lumiblock;
   UInt_t          trigger;
   vector<int>     *tri_PAL2DoubleMu0;
   vector<int>     *tri_Multiplicity150;
   vector<int>     *tri_Multiplicity185;
   vector<int>     *tri_Multiplicity220;
   TLorentzVector  *gen_bc_p4;
   TLorentzVector  *gen_jpsi_p4;
   TLorentzVector  *gen_pion3_p4;
   TLorentzVector  *gen_muon1_p4;
   TLorentzVector  *gen_muon2_p4;
   TVector3        *gen_bc_vtx;
   TVector3        *gen_jpsi_vtx;
   Float_t         gen_bc_ct;

   // List of branches
   TBranch        *b_nB;   //!
   TBranch        *b_nMu;   //!
   TBranch        *b_J_px1;   //!
   TBranch        *b_J_py1;   //!
   TBranch        *b_J_pz1;   //!
   TBranch        *b_J_px2;   //!
   TBranch        *b_J_py2;   //!
   TBranch        *b_J_pz2;   //!
   TBranch        *b_B_k_px;   //!
   TBranch        *b_B_k_py;   //!
   TBranch        *b_B_k_pz;   //!
   TBranch        *b_pi1Q;   //!
   TBranch        *b_pi1_trackerhits;   //!
   TBranch        *b_pi1_pixelhits;   //!
   TBranch        *b_pi1dz;   //!
   TBranch        *b_pi1dzE;   //!
   TBranch        *b_pi1dxy;   //!
   TBranch        *b_pi1dxyE;   //!
   TBranch        *b_B_charge;   //!
   TBranch        *b_B_mass;   //!
   TBranch        *b_B_px;   //!
   TBranch        *b_B_py;   //!
   TBranch        *b_B_pz;   //!
   TBranch        *b_B_k_charge1;   //!
   TBranch        *b_B_k_px_track;   //!
   TBranch        *b_B_k_py_track;   //!
   TBranch        *b_B_k_pz_track;   //!
   TBranch        *b_B_J_mass;   //!
   TBranch        *b_B_J_px;   //!
   TBranch        *b_B_J_py;   //!
   TBranch        *b_B_J_pz;   //!
   TBranch        *b_B_J_px1;   //!
   TBranch        *b_B_J_py1;   //!
   TBranch        *b_B_J_pz1;   //!
   TBranch        *b_B_J_charge1;   //!
   TBranch        *b_B_J_px2;   //!
   TBranch        *b_B_J_py2;   //!
   TBranch        *b_B_J_pz2;   //!
   TBranch        *b_B_J_charge2;   //!
   TBranch        *b_muon1Q;   //!
   TBranch        *b_muon2Q;   //!
   TBranch        *b_muon1Soft;   //!
   TBranch        *b_muon2Soft;   //!
   TBranch        *b_muon1StationTight;   //!
   TBranch        *b_muon2StationTight;   //!
   TBranch        *b_mu1_trackerhits;   //!
   TBranch        *b_mu2_trackerhits;   //!
   TBranch        *b_mu1_pixelhits;   //!
   TBranch        *b_mu2_pixelhits;   //!
   TBranch        *b_mu1dz;   //!
   TBranch        *b_mu1dzE;   //!
   TBranch        *b_mu1dxy;   //!
   TBranch        *b_mu1dxyE;   //!
   TBranch        *b_mu2dz;   //!
   TBranch        *b_mu2dzE;   //!
   TBranch        *b_mu2dxy;   //!
   TBranch        *b_mu2dxyE;   //!
   TBranch        *b_B_Prob;   //!
   TBranch        *b_B_J_Prob;   //!
   TBranch        *b_B_DecayVtxX;   //!
   TBranch        *b_B_DecayVtxY;   //!
   TBranch        *b_B_DecayVtxZ;   //!
   TBranch        *b_B_DecayVtxXE;   //!
   TBranch        *b_B_DecayVtxYE;   //!
   TBranch        *b_B_DecayVtxZE;   //!
   TBranch        *b_B_DecayVtxXYE;   //!
   TBranch        *b_B_DecayVtxXZE;   //!
   TBranch        *b_B_DecayVtxYZE;   //!
   TBranch        *b_B_J_DecayVtxX;   //!
   TBranch        *b_B_J_DecayVtxY;   //!
   TBranch        *b_B_J_DecayVtxZ;   //!
   TBranch        *b_B_J_DecayVtxXE;   //!
   TBranch        *b_B_J_DecayVtxYE;   //!
   TBranch        *b_B_J_DecayVtxZE;   //!
   TBranch        *b_B_J_DecayVtxXYE;   //!
   TBranch        *b_B_J_DecayVtxXZE;   //!
   TBranch        *b_B_J_DecayVtxYZE;   //!
   TBranch        *b_priVtxX;   //!
   TBranch        *b_priVtxY;   //!
   TBranch        *b_priVtxZ;   //!
   TBranch        *b_priVtxXE;   //!
   TBranch        *b_priVtxYE;   //!
   TBranch        *b_priVtxZE;   //!
   TBranch        *b_priVtxXYE;   //!
   TBranch        *b_priVtxXZE;   //!
   TBranch        *b_priVtxYZE;   //!
   TBranch        *b_priVtxCL;   //!
   TBranch        *b_pVtxIPX;   //!
   TBranch        *b_pVtxIPY;   //!
   TBranch        *b_pVtxIPZ;   //!
   TBranch        *b_pVtxIPXE;   //!
   TBranch        *b_pVtxIPYE;   //!
   TBranch        *b_pVtxIPZE;   //!
   TBranch        *b_pVtxIPXYE;   //!
   TBranch        *b_pVtxIPXZE;   //!
   TBranch        *b_pVtxIPYZE;   //!
   TBranch        *b_pVtxIPCL;   //!
   TBranch        *b_priRfVtxX;   //!
   TBranch        *b_priRfVtxY;   //!
   TBranch        *b_priRfVtxZ;   //!
   TBranch        *b_priRfVtxXE;   //!
   TBranch        *b_priRfVtxYE;   //!
   TBranch        *b_priRfVtxZE;   //!
   TBranch        *b_priRfVtxXYE;   //!
   TBranch        *b_priRfVtxXZE;   //!
   TBranch        *b_priRfVtxYZE;   //!
   TBranch        *b_priRfVtxCL;   //!
   TBranch        *b_priRfNTrkDif;   //!
   TBranch        *b_nVtx;   //!
   TBranch        *b_nTrk;   //!
   TBranch        *b_nTrk_Vtx;   //!
   TBranch        *b_nTrk_VtxIP;   //!
   TBranch        *b_nTrk_VtxIP_Q;   //!
   TBranch        *b_run;   //!
   TBranch        *b_event;   //!
   TBranch        *b_lumiblock;   //!
   TBranch        *b_trigger;   //!
   TBranch        *b_tri_PAL2DoubleMu0;   //!
   TBranch        *b_tri_Multiplicity150;   //!
   TBranch        *b_tri_Multiplicity185;   //!
   TBranch        *b_tri_Multiplicity220;   //!
   TBranch        *b_gen_bc_p4;   //!
   TBranch        *b_gen_jpsi_p4;   //!
   TBranch        *b_gen_pion3_p4;   //!
   TBranch        *b_gen_muon1_p4;   //!
   TBranch        *b_gen_muon2_p4;   //!
   TBranch        *b_gen_bc_vtx;   //!
   TBranch        *b_gen_jpsi_vtx;   //!
   TBranch        *b_gen_bc_ct;   //!

   classntupleMC(TTree *tree=0);
   virtual ~classntupleMC();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef classntupleMC_cxx
classntupleMC::classntupleMC(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("dataMC/Rootuple_Bu_JpsiK_2016HI-MCAOD_bacht1.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("dataMC/Rootuple_Bu_JpsiK_2016HI-MCAOD_bacht1.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("dataMC/Rootuple_Bu_JpsiK_2016HI-MCAOD_bacht1.root:/rootuple");
      dir->GetObject("ntuple",tree);

   }
   Init(tree);
}

classntupleMC::~classntupleMC()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t classntupleMC::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t classntupleMC::LoadTree(Long64_t entry)
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

void classntupleMC::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   J_px1 = 0;
   J_py1 = 0;
   J_pz1 = 0;
   J_px2 = 0;
   J_py2 = 0;
   J_pz2 = 0;
   B_k_px = 0;
   B_k_py = 0;
   B_k_pz = 0;
   pi1Q = 0;
   pi1_trackerhits = 0;
   pi1_pixelhits = 0;
   pi1dz = 0;
   pi1dzE = 0;
   pi1dxy = 0;
   pi1dxyE = 0;
   B_charge = 0;
   B_mass = 0;
   B_px = 0;
   B_py = 0;
   B_pz = 0;
   B_k_charge1 = 0;
   B_k_px_track = 0;
   B_k_py_track = 0;
   B_k_pz_track = 0;
   B_J_mass = 0;
   B_J_px = 0;
   B_J_py = 0;
   B_J_pz = 0;
   B_J_px1 = 0;
   B_J_py1 = 0;
   B_J_pz1 = 0;
   B_J_charge1 = 0;
   B_J_px2 = 0;
   B_J_py2 = 0;
   B_J_pz2 = 0;
   B_J_charge2 = 0;
   muon1Q = 0;
   muon2Q = 0;
   muon1Soft = 0;
   muon2Soft = 0;
   muon1StationTight = 0;
   muon2StationTight = 0;
   mu1_trackerhits = 0;
   mu2_trackerhits = 0;
   mu1_pixelhits = 0;
   mu2_pixelhits = 0;
   mu1dz = 0;
   mu1dzE = 0;
   mu1dxy = 0;
   mu1dxyE = 0;
   mu2dz = 0;
   mu2dzE = 0;
   mu2dxy = 0;
   mu2dxyE = 0;
   B_Prob = 0;
   B_J_Prob = 0;
   B_DecayVtxX = 0;
   B_DecayVtxY = 0;
   B_DecayVtxZ = 0;
   B_DecayVtxXE = 0;
   B_DecayVtxYE = 0;
   B_DecayVtxZE = 0;
   B_DecayVtxXYE = 0;
   B_DecayVtxXZE = 0;
   B_DecayVtxYZE = 0;
   B_J_DecayVtxX = 0;
   B_J_DecayVtxY = 0;
   B_J_DecayVtxZ = 0;
   B_J_DecayVtxXE = 0;
   B_J_DecayVtxYE = 0;
   B_J_DecayVtxZE = 0;
   B_J_DecayVtxXYE = 0;
   B_J_DecayVtxXZE = 0;
   B_J_DecayVtxYZE = 0;
   pVtxIPX = 0;
   pVtxIPY = 0;
   pVtxIPZ = 0;
   pVtxIPXE = 0;
   pVtxIPYE = 0;
   pVtxIPZE = 0;
   pVtxIPXYE = 0;
   pVtxIPXZE = 0;
   pVtxIPYZE = 0;
   pVtxIPCL = 0;
   priRfVtxX = 0;
   priRfVtxY = 0;
   priRfVtxZ = 0;
   priRfVtxXE = 0;
   priRfVtxYE = 0;
   priRfVtxZE = 0;
   priRfVtxXYE = 0;
   priRfVtxXZE = 0;
   priRfVtxYZE = 0;
   priRfVtxCL = 0;
   priRfNTrkDif = 0;
   nTrk_VtxIP = 0;
   nTrk_VtxIP_Q = 0;
   tri_PAL2DoubleMu0 = 0;
   tri_Multiplicity150 = 0;
   tri_Multiplicity185 = 0;
   tri_Multiplicity220 = 0;
   gen_bc_p4 = 0;
   gen_jpsi_p4 = 0;
   gen_pion3_p4 = 0;
   gen_muon1_p4 = 0;
   gen_muon2_p4 = 0;
   gen_bc_vtx = 0;
   gen_jpsi_vtx = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("nB", &nB, &b_nB);
   fChain->SetBranchAddress("nMu", &nMu, &b_nMu);
   fChain->SetBranchAddress("J_px1", &J_px1, &b_J_px1);
   fChain->SetBranchAddress("J_py1", &J_py1, &b_J_py1);
   fChain->SetBranchAddress("J_pz1", &J_pz1, &b_J_pz1);
   fChain->SetBranchAddress("J_px2", &J_px2, &b_J_px2);
   fChain->SetBranchAddress("J_py2", &J_py2, &b_J_py2);
   fChain->SetBranchAddress("J_pz2", &J_pz2, &b_J_pz2);
   fChain->SetBranchAddress("B_k_px", &B_k_px, &b_B_k_px);
   fChain->SetBranchAddress("B_k_py", &B_k_py, &b_B_k_py);
   fChain->SetBranchAddress("B_k_pz", &B_k_pz, &b_B_k_pz);
   fChain->SetBranchAddress("pi1Q", &pi1Q, &b_pi1Q);
   fChain->SetBranchAddress("pi1_trackerhits", &pi1_trackerhits, &b_pi1_trackerhits);
   fChain->SetBranchAddress("pi1_pixelhits", &pi1_pixelhits, &b_pi1_pixelhits);
   fChain->SetBranchAddress("pi1dz", &pi1dz, &b_pi1dz);
   fChain->SetBranchAddress("pi1dzE", &pi1dzE, &b_pi1dzE);
   fChain->SetBranchAddress("pi1dxy", &pi1dxy, &b_pi1dxy);
   fChain->SetBranchAddress("pi1dxyE", &pi1dxyE, &b_pi1dxyE);
   fChain->SetBranchAddress("B_charge", &B_charge, &b_B_charge);
   fChain->SetBranchAddress("B_mass", &B_mass, &b_B_mass);
   fChain->SetBranchAddress("B_px", &B_px, &b_B_px);
   fChain->SetBranchAddress("B_py", &B_py, &b_B_py);
   fChain->SetBranchAddress("B_pz", &B_pz, &b_B_pz);
   fChain->SetBranchAddress("B_k_charge1", &B_k_charge1, &b_B_k_charge1);
   fChain->SetBranchAddress("B_k_px_track", &B_k_px_track, &b_B_k_px_track);
   fChain->SetBranchAddress("B_k_py_track", &B_k_py_track, &b_B_k_py_track);
   fChain->SetBranchAddress("B_k_pz_track", &B_k_pz_track, &b_B_k_pz_track);
   fChain->SetBranchAddress("B_J_mass", &B_J_mass, &b_B_J_mass);
   fChain->SetBranchAddress("B_J_px", &B_J_px, &b_B_J_px);
   fChain->SetBranchAddress("B_J_py", &B_J_py, &b_B_J_py);
   fChain->SetBranchAddress("B_J_pz", &B_J_pz, &b_B_J_pz);
   fChain->SetBranchAddress("B_J_px1", &B_J_px1, &b_B_J_px1);
   fChain->SetBranchAddress("B_J_py1", &B_J_py1, &b_B_J_py1);
   fChain->SetBranchAddress("B_J_pz1", &B_J_pz1, &b_B_J_pz1);
   fChain->SetBranchAddress("B_J_charge1", &B_J_charge1, &b_B_J_charge1);
   fChain->SetBranchAddress("B_J_px2", &B_J_px2, &b_B_J_px2);
   fChain->SetBranchAddress("B_J_py2", &B_J_py2, &b_B_J_py2);
   fChain->SetBranchAddress("B_J_pz2", &B_J_pz2, &b_B_J_pz2);
   fChain->SetBranchAddress("B_J_charge2", &B_J_charge2, &b_B_J_charge2);
   fChain->SetBranchAddress("muon1Q", &muon1Q, &b_muon1Q);
   fChain->SetBranchAddress("muon2Q", &muon2Q, &b_muon2Q);
   fChain->SetBranchAddress("muon1Soft", &muon1Soft, &b_muon1Soft);
   fChain->SetBranchAddress("muon2Soft", &muon2Soft, &b_muon2Soft);
   fChain->SetBranchAddress("muon1StationTight", &muon1StationTight, &b_muon1StationTight);
   fChain->SetBranchAddress("muon2StationTight", &muon2StationTight, &b_muon2StationTight);
   fChain->SetBranchAddress("mu1_trackerhits", &mu1_trackerhits, &b_mu1_trackerhits);
   fChain->SetBranchAddress("mu2_trackerhits", &mu2_trackerhits, &b_mu2_trackerhits);
   fChain->SetBranchAddress("mu1_pixelhits", &mu1_pixelhits, &b_mu1_pixelhits);
   fChain->SetBranchAddress("mu2_pixelhits", &mu2_pixelhits, &b_mu2_pixelhits);
   fChain->SetBranchAddress("mu1dz", &mu1dz, &b_mu1dz);
   fChain->SetBranchAddress("mu1dzE", &mu1dzE, &b_mu1dzE);
   fChain->SetBranchAddress("mu1dxy", &mu1dxy, &b_mu1dxy);
   fChain->SetBranchAddress("mu1dxyE", &mu1dxyE, &b_mu1dxyE);
   fChain->SetBranchAddress("mu2dz", &mu2dz, &b_mu2dz);
   fChain->SetBranchAddress("mu2dzE", &mu2dzE, &b_mu2dzE);
   fChain->SetBranchAddress("mu2dxy", &mu2dxy, &b_mu2dxy);
   fChain->SetBranchAddress("mu2dxyE", &mu2dxyE, &b_mu2dxyE);
   fChain->SetBranchAddress("B_Prob", &B_Prob, &b_B_Prob);
   fChain->SetBranchAddress("B_J_Prob", &B_J_Prob, &b_B_J_Prob);
   fChain->SetBranchAddress("B_DecayVtxX", &B_DecayVtxX, &b_B_DecayVtxX);
   fChain->SetBranchAddress("B_DecayVtxY", &B_DecayVtxY, &b_B_DecayVtxY);
   fChain->SetBranchAddress("B_DecayVtxZ", &B_DecayVtxZ, &b_B_DecayVtxZ);
   fChain->SetBranchAddress("B_DecayVtxXE", &B_DecayVtxXE, &b_B_DecayVtxXE);
   fChain->SetBranchAddress("B_DecayVtxYE", &B_DecayVtxYE, &b_B_DecayVtxYE);
   fChain->SetBranchAddress("B_DecayVtxZE", &B_DecayVtxZE, &b_B_DecayVtxZE);
   fChain->SetBranchAddress("B_DecayVtxXYE", &B_DecayVtxXYE, &b_B_DecayVtxXYE);
   fChain->SetBranchAddress("B_DecayVtxXZE", &B_DecayVtxXZE, &b_B_DecayVtxXZE);
   fChain->SetBranchAddress("B_DecayVtxYZE", &B_DecayVtxYZE, &b_B_DecayVtxYZE);
   fChain->SetBranchAddress("B_J_DecayVtxX", &B_J_DecayVtxX, &b_B_J_DecayVtxX);
   fChain->SetBranchAddress("B_J_DecayVtxY", &B_J_DecayVtxY, &b_B_J_DecayVtxY);
   fChain->SetBranchAddress("B_J_DecayVtxZ", &B_J_DecayVtxZ, &b_B_J_DecayVtxZ);
   fChain->SetBranchAddress("B_J_DecayVtxXE", &B_J_DecayVtxXE, &b_B_J_DecayVtxXE);
   fChain->SetBranchAddress("B_J_DecayVtxYE", &B_J_DecayVtxYE, &b_B_J_DecayVtxYE);
   fChain->SetBranchAddress("B_J_DecayVtxZE", &B_J_DecayVtxZE, &b_B_J_DecayVtxZE);
   fChain->SetBranchAddress("B_J_DecayVtxXYE", &B_J_DecayVtxXYE, &b_B_J_DecayVtxXYE);
   fChain->SetBranchAddress("B_J_DecayVtxXZE", &B_J_DecayVtxXZE, &b_B_J_DecayVtxXZE);
   fChain->SetBranchAddress("B_J_DecayVtxYZE", &B_J_DecayVtxYZE, &b_B_J_DecayVtxYZE);
   fChain->SetBranchAddress("priVtxX", &priVtxX, &b_priVtxX);
   fChain->SetBranchAddress("priVtxY", &priVtxY, &b_priVtxY);
   fChain->SetBranchAddress("priVtxZ", &priVtxZ, &b_priVtxZ);
   fChain->SetBranchAddress("priVtxXE", &priVtxXE, &b_priVtxXE);
   fChain->SetBranchAddress("priVtxYE", &priVtxYE, &b_priVtxYE);
   fChain->SetBranchAddress("priVtxZE", &priVtxZE, &b_priVtxZE);
   fChain->SetBranchAddress("priVtxXYE", &priVtxXYE, &b_priVtxXYE);
   fChain->SetBranchAddress("priVtxXZE", &priVtxXZE, &b_priVtxXZE);
   fChain->SetBranchAddress("priVtxYZE", &priVtxYZE, &b_priVtxYZE);
   fChain->SetBranchAddress("priVtxCL", &priVtxCL, &b_priVtxCL);
   fChain->SetBranchAddress("pVtxIPX", &pVtxIPX, &b_pVtxIPX);
   fChain->SetBranchAddress("pVtxIPY", &pVtxIPY, &b_pVtxIPY);
   fChain->SetBranchAddress("pVtxIPZ", &pVtxIPZ, &b_pVtxIPZ);
   fChain->SetBranchAddress("pVtxIPXE", &pVtxIPXE, &b_pVtxIPXE);
   fChain->SetBranchAddress("pVtxIPYE", &pVtxIPYE, &b_pVtxIPYE);
   fChain->SetBranchAddress("pVtxIPZE", &pVtxIPZE, &b_pVtxIPZE);
   fChain->SetBranchAddress("pVtxIPXYE", &pVtxIPXYE, &b_pVtxIPXYE);
   fChain->SetBranchAddress("pVtxIPXZE", &pVtxIPXZE, &b_pVtxIPXZE);
   fChain->SetBranchAddress("pVtxIPYZE", &pVtxIPYZE, &b_pVtxIPYZE);
   fChain->SetBranchAddress("pVtxIPCL", &pVtxIPCL, &b_pVtxIPCL);
   fChain->SetBranchAddress("priRfVtxX", &priRfVtxX, &b_priRfVtxX);
   fChain->SetBranchAddress("priRfVtxY", &priRfVtxY, &b_priRfVtxY);
   fChain->SetBranchAddress("priRfVtxZ", &priRfVtxZ, &b_priRfVtxZ);
   fChain->SetBranchAddress("priRfVtxXE", &priRfVtxXE, &b_priRfVtxXE);
   fChain->SetBranchAddress("priRfVtxYE", &priRfVtxYE, &b_priRfVtxYE);
   fChain->SetBranchAddress("priRfVtxZE", &priRfVtxZE, &b_priRfVtxZE);
   fChain->SetBranchAddress("priRfVtxXYE", &priRfVtxXYE, &b_priRfVtxXYE);
   fChain->SetBranchAddress("priRfVtxXZE", &priRfVtxXZE, &b_priRfVtxXZE);
   fChain->SetBranchAddress("priRfVtxYZE", &priRfVtxYZE, &b_priRfVtxYZE);
   fChain->SetBranchAddress("priRfVtxCL", &priRfVtxCL, &b_priRfVtxCL);
   fChain->SetBranchAddress("priRfNTrkDif", &priRfNTrkDif, &b_priRfNTrkDif);
   fChain->SetBranchAddress("nVtx", &nVtx, &b_nVtx);
   fChain->SetBranchAddress("nTrk", &nTrk, &b_nTrk);
   fChain->SetBranchAddress("nTrk_Vtx", &nTrk_Vtx, &b_nTrk_Vtx);
   fChain->SetBranchAddress("nTrk_VtxIP", &nTrk_VtxIP, &b_nTrk_VtxIP);
   fChain->SetBranchAddress("nTrk_VtxIP_Q", &nTrk_VtxIP_Q, &b_nTrk_VtxIP_Q);
   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("lumiblock", &lumiblock, &b_lumiblock);
   fChain->SetBranchAddress("trigger", &trigger, &b_trigger);
   fChain->SetBranchAddress("tri_PAL2DoubleMu0", &tri_PAL2DoubleMu0, &b_tri_PAL2DoubleMu0);
   fChain->SetBranchAddress("tri_Multiplicity150", &tri_Multiplicity150, &b_tri_Multiplicity150);
   fChain->SetBranchAddress("tri_Multiplicity185", &tri_Multiplicity185, &b_tri_Multiplicity185);
   fChain->SetBranchAddress("tri_Multiplicity220", &tri_Multiplicity220, &b_tri_Multiplicity220);
   fChain->SetBranchAddress("gen_bc_p4", &gen_bc_p4, &b_gen_bc_p4);
   fChain->SetBranchAddress("gen_jpsi_p4", &gen_jpsi_p4, &b_gen_jpsi_p4);
   fChain->SetBranchAddress("gen_pion3_p4", &gen_pion3_p4, &b_gen_pion3_p4);
   fChain->SetBranchAddress("gen_muon1_p4", &gen_muon1_p4, &b_gen_muon1_p4);
   fChain->SetBranchAddress("gen_muon2_p4", &gen_muon2_p4, &b_gen_muon2_p4);
   fChain->SetBranchAddress("gen_bc_vtx", &gen_bc_vtx, &b_gen_bc_vtx);
   fChain->SetBranchAddress("gen_jpsi_vtx", &gen_jpsi_vtx, &b_gen_jpsi_vtx);
   fChain->SetBranchAddress("gen_bc_ct", &gen_bc_ct, &b_gen_bc_ct);
   Notify();
}

Bool_t classntupleMC::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void classntupleMC::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t classntupleMC::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef classntupleMC_cxx
