//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Nov 26 03:11:41 2020 by ROOT version 6.16/01
// from TTree treeS/signal
// found on file: ROOTS_Bujk_AOD_HI2016_Era1_best1.root
//////////////////////////////////////////////////////////

#ifndef classmake_h
#define classmake_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class classmake {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Double_t        massB;
   Double_t        Bpt;
   Double_t        Brapidity;
   Double_t        Bphi;
   Double_t        Bdl;
   Double_t        BdlE;
   Double_t        BdlIP;
   Double_t        BdlIPE;
   Double_t        massJ;
   Double_t        Jpsipt;
   Double_t        Jrapidity;
   Double_t        Jphi;
   Double_t        sigLxyJ;
   Double_t        cosalfaJ;
   Double_t        sigLxyB;
   Double_t        cosalfaB;
   Int_t           Kcharge;
   Double_t        pi1pt;
   Double_t        pi1eta;
   Double_t        pi1phi;
   Int_t           pi1Thits;
   Int_t           pi1Phits;
   Int_t           pi1q;
   Double_t        mu1pt;
   Double_t        mu2pt;
   Double_t        mu1phi;
   Double_t        mu2phi;
   Double_t        mu1eta;
   Double_t        mu2eta;
   Int_t           mu1soft;
   Int_t           mu2soft;
   Int_t           muon1StationTight;
   Int_t           muon2StationTight;
   Double_t        maxPtMu;
   Double_t        minPtMu;
   Double_t        Bpro;
   Double_t        Jpro;
   Double_t        dxypi1;
   Double_t        dxyEpi1;
   Double_t        DRmumu;
   Int_t           Ntrk;
   Int_t           Ntrkvip;
   Int_t           NtrkvipQ;
   UInt_t          nPV;
   Int_t           runn;
   Int_t           evtn;
   Int_t           Lumiblock;
   Int_t           Tripal1open;
   Int_t           Tripal1;
   Int_t           Tripal2;
   Int_t           Tripal3;
   UInt_t          Trigger;

   // List of branches
   TBranch        *b_massB;   //!
   TBranch        *b_Bpt;   //!
   TBranch        *b_Brapidity;   //!
   TBranch        *b_Bphi;   //!
   TBranch        *b_Bdl;   //!
   TBranch        *b_BdlE;   //!
   TBranch        *b_BdlIP;   //!
   TBranch        *b_BdlIPE;   //!
   TBranch        *b_massJ;   //!
   TBranch        *b_Jpsipt;   //!
   TBranch        *b_Jrapidity;   //!
   TBranch        *b_Jphi;   //!
   TBranch        *b_sigLxyJ;   //!
   TBranch        *b_cosalfaJ;   //!
   TBranch        *b_sigLxyB;   //!
   TBranch        *b_cosalfaB;   //!
   TBranch        *b_Kcharge;   //!
   TBranch        *b_pi1pt;   //!
   TBranch        *b_pi1eta;   //!
   TBranch        *b_pi1phi;   //!
   TBranch        *b_pi1Thits;   //!
   TBranch        *b_pi1Phits;   //!
   TBranch        *b_pi1q;   //!
   TBranch        *b_mu1pt;   //!
   TBranch        *b_mu2pt;   //!
   TBranch        *b_mu1phi;   //!
   TBranch        *b_mu2phi;   //!
   TBranch        *b_mu1eta;   //!
   TBranch        *b_mu2eta;   //!
   TBranch        *b_mu1soft;   //!
   TBranch        *b_mu2soft;   //!
   TBranch        *b_muon1StationTight;   //!
   TBranch        *b_muon2StationTight;   //!
   TBranch        *b_maxPtMu;   //!
   TBranch        *b_minPtMu;   //!
   TBranch        *b_Bpro;   //!
   TBranch        *b_Jpro;   //!
   TBranch        *b_dxypi1;   //!
   TBranch        *b_dxyEpi1;   //!
   TBranch        *b_DRmumu;   //!
   TBranch        *b_Ntrk;   //!
   TBranch        *b_Ntrkvip;   //!
   TBranch        *b_NtrkvipQ;   //!
   TBranch        *b_nPV;   //!
   TBranch        *b_runn;   //!
   TBranch        *b_evtn;   //!
   TBranch        *b_Lumiblock;   //!
   TBranch        *b_Tripal1open;   //!
   TBranch        *b_Tripal1;   //!
   TBranch        *b_Tripal2;   //!
   TBranch        *b_Tripal3;   //!
   TBranch        *b_Trigger;   //!

   classmake(TTree *tree=0);
   virtual ~classmake();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef classmake_cxx
classmake::classmake(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("ROOTS_Bujk_AOD_HI2016_Era1_best1.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("ROOTS_Bujk_AOD_HI2016_Era1_best1.root");
      }
      f->GetObject("treeS",tree);

   }
   Init(tree);
}

classmake::~classmake()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t classmake::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t classmake::LoadTree(Long64_t entry)
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

void classmake::Init(TTree *tree)
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

   fChain->SetBranchAddress("massB", &massB, &b_massB);
   fChain->SetBranchAddress("Bpt", &Bpt, &b_Bpt);
   fChain->SetBranchAddress("Brapidity", &Brapidity, &b_Brapidity);
   fChain->SetBranchAddress("Bphi", &Bphi, &b_Bphi);
   fChain->SetBranchAddress("Bdl", &Bdl, &b_Bdl);
   fChain->SetBranchAddress("BdlE", &BdlE, &b_BdlE);
   fChain->SetBranchAddress("BdlIP", &BdlIP, &b_BdlIP);
   fChain->SetBranchAddress("BdlIPE", &BdlIPE, &b_BdlIPE);
   fChain->SetBranchAddress("massJ", &massJ, &b_massJ);
   fChain->SetBranchAddress("Jpsipt", &Jpsipt, &b_Jpsipt);
   fChain->SetBranchAddress("Jrapidity", &Jrapidity, &b_Jrapidity);
   fChain->SetBranchAddress("Jphi", &Jphi, &b_Jphi);
   fChain->SetBranchAddress("sigLxyJ", &sigLxyJ, &b_sigLxyJ);
   fChain->SetBranchAddress("cosalfaJ", &cosalfaJ, &b_cosalfaJ);
   fChain->SetBranchAddress("sigLxyB", &sigLxyB, &b_sigLxyB);
   fChain->SetBranchAddress("cosalfaB", &cosalfaB, &b_cosalfaB);
   fChain->SetBranchAddress("Kcharge", &Kcharge, &b_Kcharge);
   fChain->SetBranchAddress("pi1pt", &pi1pt, &b_pi1pt);
   fChain->SetBranchAddress("pi1eta", &pi1eta, &b_pi1eta);
   fChain->SetBranchAddress("pi1phi", &pi1phi, &b_pi1phi);
   fChain->SetBranchAddress("pi1Thits", &pi1Thits, &b_pi1Thits);
   fChain->SetBranchAddress("pi1Phits", &pi1Phits, &b_pi1Phits);
   fChain->SetBranchAddress("pi1q", &pi1q, &b_pi1q);
   fChain->SetBranchAddress("mu1pt", &mu1pt, &b_mu1pt);
   fChain->SetBranchAddress("mu2pt", &mu2pt, &b_mu2pt);
   fChain->SetBranchAddress("mu1phi", &mu1phi, &b_mu1phi);
   fChain->SetBranchAddress("mu2phi", &mu2phi, &b_mu2phi);
   fChain->SetBranchAddress("mu1eta", &mu1eta, &b_mu1eta);
   fChain->SetBranchAddress("mu2eta", &mu2eta, &b_mu2eta);
   fChain->SetBranchAddress("mu1soft", &mu1soft, &b_mu1soft);
   fChain->SetBranchAddress("mu2soft", &mu2soft, &b_mu2soft);
   fChain->SetBranchAddress("muon1StationTight", &muon1StationTight, &b_muon1StationTight);
   fChain->SetBranchAddress("muon2StationTight", &muon2StationTight, &b_muon2StationTight);
   fChain->SetBranchAddress("maxPtMu", &maxPtMu, &b_maxPtMu);
   fChain->SetBranchAddress("minPtMu", &minPtMu, &b_minPtMu);
   fChain->SetBranchAddress("Bpro", &Bpro, &b_Bpro);
   fChain->SetBranchAddress("Jpro", &Jpro, &b_Jpro);
   fChain->SetBranchAddress("dxypi1", &dxypi1, &b_dxypi1);
   fChain->SetBranchAddress("dxyEpi1", &dxyEpi1, &b_dxyEpi1);
   fChain->SetBranchAddress("DRmumu", &DRmumu, &b_DRmumu);
   fChain->SetBranchAddress("Ntrk", &Ntrk, &b_Ntrk);
   fChain->SetBranchAddress("Ntrkvip", &Ntrkvip, &b_Ntrkvip);
   fChain->SetBranchAddress("NtrkvipQ", &NtrkvipQ, &b_NtrkvipQ);
   fChain->SetBranchAddress("nPV", &nPV, &b_nPV);
   fChain->SetBranchAddress("runn", &runn, &b_runn);
   fChain->SetBranchAddress("evtn", &evtn, &b_evtn);
   fChain->SetBranchAddress("Lumiblock", &Lumiblock, &b_Lumiblock);
   fChain->SetBranchAddress("Tripal1open", &Tripal1open, &b_Tripal1open);
   fChain->SetBranchAddress("Tripal1", &Tripal1, &b_Tripal1);
   fChain->SetBranchAddress("Tripal2", &Tripal2, &b_Tripal2);
   fChain->SetBranchAddress("Tripal3", &Tripal3, &b_Tripal3);
   fChain->SetBranchAddress("Trigger", &Trigger, &b_Trigger);
   Notify();
}

Bool_t classmake::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void classmake::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t classmake::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef classmake_cxx
