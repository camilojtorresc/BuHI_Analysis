#include "classmakeMC.C"

#include <TROOT.h>
#include "TMath.h"
#include <iostream>
#include  <vector>
#include "TVector3.h"
#include "TMatrixD.h"
#include "TLorentzVector.h"
#include "TRandom3.h"
#include "TFile.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodBDT.h"
#include <TVector.h>

using namespace std;

void reduce_rootuples_Bujk(int vq)
{

TChain *ch = new TChain("treeS",""); 
ch->Add(Form("ROOTS_Bujk_AOD_HI2016_MC_best%1i.root/treeS",vq) ); 

TTree *tree = (TTree*) ch;
classmakeMC t(tree);
Long64_t nentries = t.fChain->GetEntries();
cout<<" Entries : "<<nentries<<endl;
cout<<"xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"<<endl;
cout<<" Reading data ..."<<nentries<<endl;

//------------------------------
TString outfileName( Form("reducetree_Bujk_AOD_HI2016_MC_best%1i.root",vq)  ); 
TFile* outputFile = TFile::Open( outfileName, "RECREATE" );
 
//*******************************
//*     VARIABLES TO SAVE       *
//*******************************

Double_t B_mass, J_mass;
Double_t Bupt, rapidityB, pdl, pdle;
Double_t pion1pt; 
Int_t    p1charg;
Double_t etapi1;  
Double_t mu1pt, mu2pt;
Double_t etamu1, etamu2;
Double_t Bpro, Jpro;
Int_t Ntrk, Ntrkvip, NtrkvipQ; 
Int_t trig;
UInt_t Trigger; 
Int_t runn, evtn;

Double_t pi1eta, pi1phi;
TLorentzVector  *gen_pion1_p4=0;
TLorentzVector  *gen_muon1_p4=0;
TLorentzVector  *gen_muon2_p4=0;
TLorentzVector  *gen_b_p4=0;
 
TTree *butree =  new TTree("butree","butree");
butree->Branch("B_mass", &B_mass,  "B_mass/D");
butree->Branch("J_mass", &J_mass,  "J_mass/D");
butree->Branch("Bupt", &Bupt,  "Bupt/D"); 
butree->Branch("rapidityB", &rapidityB,  "rapidityB/D");
butree->Branch("pdl", &pdl,  "pdl/D");
butree->Branch("pdle", &pdle,  "pdle/D");
butree->Branch("pion1pt", &pion1pt,  "pion1pt/D");
butree->Branch("p1charg", &p1charg,  "p1charg/I");
butree->Branch("etapi1",  &etapi1,    "etapi1/D");  
butree->Branch("etamu1",  &etamu1,    "etamu1/D");
butree->Branch("etamu2",  &etamu2,    "etamu2/D");
butree->Branch("mu1pt", &mu1pt,  "mu1pt/D");
butree->Branch("mu2pt", &mu2pt,  "mu2pt/D");
butree->Branch("Bpro",   &Bpro,    "Bpro/D");
butree->Branch("Jpro",   &Jpro,    "Jpro/D");
butree->Branch("Ntrk", &Ntrk, "Ntrk/I");
butree->Branch("Ntrkvip", &Ntrkvip, "Ntrkvip/I");
butree->Branch("NtrkvipQ", &NtrkvipQ, "NtrkvipQ/I"); 
butree->Branch("trig", &trig, "trig/I");
butree->Branch("Trigger", &Trigger, "Trigger/i");
butree->Branch("runn", &runn, "runn/I");
butree->Branch("evtn", &evtn, "evtn/I");

butree->Branch("pi1eta", &pi1eta,  "pi1eta/D");
butree->Branch("pi1phi", &pi1phi,  "pi1phi/D");
butree->Branch("gen_pion1_p4", &gen_pion1_p4);
butree->Branch("gen_muon1_p4", &gen_muon1_p4);
butree->Branch("gen_muon2_p4", &gen_muon2_p4);
butree->Branch("gen_b_p4", &gen_b_p4);

Int_t nTen = nentries/10;
Int_t k=0;
Int_t nbytes = 0, nb = 0;
for(Long64_t jentry=0; jentry<nentries;jentry++)
    {
      Long64_t ientry = t.LoadTree(jentry);
      if (ientry < 0) break;
      nb = t.fChain->GetEntry(jentry);   nbytes += nb;
      if(jentry%nTen==0) cout<<10*(jentry/nTen)<<"%-"<<flush;
      if(jentry==nentries-1) cout<<endl;

      //butree->Draw("B_mass","pion1pt>0.5 && mu1pt>3.0 && mu2pt>3.0 && (pdl/pdle)>5.0")
      //*** Mass windows cuts *****
      if(t.massB<=5.0 || t.massB>=5.6) continue;
      if(t.massJ<=2.9 || t.massJ>=3.3) continue;

      //*******  Jpsi cuts  **************
      //if(t.DRmumu>1.2) continue;
      if(t.mu1pt<1.5) continue;
      if(t.mu2pt<1.5) continue;
      if(t.Jpro<0.01)continue;   
      //if(t.Jpsipt<7.9) continue;
      //if(t.sigLxyJ<3.0)continue;
      //if(t.cosalfaJ<0.9)continue;

      if(t.mu1soft!=1 || t.mu2soft!=1)continue;
      //if(t.muon1StationTight!=1 || t.muon2StationTight!=1)continue;

      Int_t trigtmp = 1;
      //if(t.Tripal1open!=1 && t.Tripal1!=1 && t.Tripal2!=1 && t.Tripal3!=1)continue;
      //if(t.Tripal1open!=1 && t.Tripal1!=1 && t.Tripal2!=1 && t.Tripal3!=1){
      //if(t.Tripal2!=1 && t.Tripal3!=1){
      if(t.Tripal2!=1){		
	trigtmp=0;
      }
      
      //************  Bu cuts  ************
      if(t.Bpt<1.0) continue;
      //if(t.Bpro<0.05)continue;  
      //if((t.BdlIP/t.BdlIPE)<4.0) continue;
      //if(t.BdlIP<0.005) continue;      
      //if(t.cosalfaB<0.99)continue;

      //************** kaon cuts ***********
      if(t.pi1pt<0.5) continue;
      if(t.pi1Thits<6)continue;
      if(t.pi1Phits<1)continue;
      if(t.pi1q!=1)continue;


      //**************rapidity  cuts ***********
      if( abs(t.mu1eta)>2.4)continue;
      if( abs(t.mu2eta)>2.4)continue;
      if( abs(t.pi1eta)>2.4)continue;
      if( abs(t.Brapidity)>2.4) continue;

      //*******************
      //*   fill tree    *
      //*******************

      B_mass = t.massB;
      J_mass  = t.massJ;  
      Bupt    = t.Bpt;
      rapidityB = t.Brapidity;
      pdl = t.BdlIP;
      pdle = t.BdlIPE;
      pion1pt = t.pi1pt;
      p1charg = t.Kcharge;
      etapi1 = t.pi1eta;      
      etamu1  = t.mu1eta;
      etamu2  = t.mu2eta;
      mu1pt = t.mu1pt;
      mu2pt = t.mu2pt;
      Bpro  = t.Bpro;
      Jpro  = t.Jpro;
      Ntrk   = t.Ntrk;
      Ntrkvip = t.Ntrkvip;
      NtrkvipQ = t.NtrkvipQ;
      trig = trigtmp;
      Trigger = t.Trigger ;
      runn= t.runn;
      evtn =t.evtn;

      pi1eta = t.pi1eta;
      pi1phi = t.pi1phi;      
      gen_pion1_p4=t.gen_pion1_p4;
      gen_muon1_p4 = t.gen_muon1_p4;
      gen_muon2_p4 = t.gen_muon2_p4;
      gen_b_p4 = t.gen_b_p4;
      
      butree->Fill();

    }

 outputFile->Write("",TObject::kOverwrite);

}//End analysis
