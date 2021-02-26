//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Feb 24 21:54:53 2021 by ROOT version 6.22/02
// from TTree butree/butree
// found on file: reducetree_Bujk_AOD_HI2016_Era1_best1.root
//////////////////////////////////////////////////////////

#ifndef classreduce_h
#define classreduce_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class classreduce {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Double_t        B_mass;
   Double_t        J_mass;
   Double_t        Bupt;
   Double_t        rapidityB;
   Double_t        pdl;
   Double_t        pdle;
   Double_t        pion1pt;
   Int_t           p1charg;
   Double_t        etapi1;
   Double_t        etamu1;
   Double_t        etamu2;
   Double_t        mu1pt;
   Double_t        mu2pt;
   Double_t        Bpro;
   Double_t        Jpro;
   Int_t           Ntrk;
   Int_t           Ntrkvip;
   Int_t           NtrkvipQ;
   Int_t           trig;
   UInt_t          Trigger;

   // List of branches
   TBranch        *b_B_mass;   //!
   TBranch        *b_J_mass;   //!
   TBranch        *b_Bupt;   //!
   TBranch        *b_rapidityB;   //!
   TBranch        *b_pdl;   //!
   TBranch        *b_pdle;   //!
   TBranch        *b_pion1pt;   //!
   TBranch        *b_p1charg;   //!
   TBranch        *b_etapi1;   //!
   TBranch        *b_etamu1;   //!
   TBranch        *b_etamu2;   //!
   TBranch        *b_mu1pt;   //!
   TBranch        *b_mu2pt;   //!
   TBranch        *b_Bpro;   //!
   TBranch        *b_Jpro;   //!
   TBranch        *b_Ntrk;   //!
   TBranch        *b_Ntrkvip;   //!
   TBranch        *b_NtrkvipQ;   //!
   TBranch        *b_trig;   //!
   TBranch        *b_Trigger;   //!

   classreduce(TTree *tree=0);
   virtual ~classreduce();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef classreduce_cxx
classreduce::classreduce(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("reducetree_Bujk_AOD_HI2016_Era1_best1.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("reducetree_Bujk_AOD_HI2016_Era1_best1.root");
      }
      f->GetObject("butree",tree);

   }
   Init(tree);
}

classreduce::~classreduce()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t classreduce::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t classreduce::LoadTree(Long64_t entry)
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

void classreduce::Init(TTree *tree)
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

   fChain->SetBranchAddress("B_mass", &B_mass, &b_B_mass);
   fChain->SetBranchAddress("J_mass", &J_mass, &b_J_mass);
   fChain->SetBranchAddress("Bupt", &Bupt, &b_Bupt);
   fChain->SetBranchAddress("rapidityB", &rapidityB, &b_rapidityB);
   fChain->SetBranchAddress("pdl", &pdl, &b_pdl);
   fChain->SetBranchAddress("pdle", &pdle, &b_pdle);
   fChain->SetBranchAddress("pion1pt", &pion1pt, &b_pion1pt);
   fChain->SetBranchAddress("p1charg", &p1charg, &b_p1charg);
   fChain->SetBranchAddress("etapi1", &etapi1, &b_etapi1);
   fChain->SetBranchAddress("etamu1", &etamu1, &b_etamu1);
   fChain->SetBranchAddress("etamu2", &etamu2, &b_etamu2);
   fChain->SetBranchAddress("mu1pt", &mu1pt, &b_mu1pt);
   fChain->SetBranchAddress("mu2pt", &mu2pt, &b_mu2pt);
   fChain->SetBranchAddress("Bpro", &Bpro, &b_Bpro);
   fChain->SetBranchAddress("Jpro", &Jpro, &b_Jpro);
   fChain->SetBranchAddress("Ntrk", &Ntrk, &b_Ntrk);
   fChain->SetBranchAddress("Ntrkvip", &Ntrkvip, &b_Ntrkvip);
   fChain->SetBranchAddress("NtrkvipQ", &NtrkvipQ, &b_NtrkvipQ);
   fChain->SetBranchAddress("trig", &trig, &b_trig);
   fChain->SetBranchAddress("Trigger", &Trigger, &b_Trigger);
   Notify();
}

Bool_t classreduce::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void classreduce::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t classreduce::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef classreduce_cxx
