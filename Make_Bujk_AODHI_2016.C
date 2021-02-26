// Primer Cambio


#include "classntuple.C"
#include <iostream>
#include  <vector>
#include "TVector3.h"
#include "TMatrixD.h"
#include "TLorentzVector.h"
#include "TRandom3.h"
#include "TFile.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodBDT.h"
#include <TMath.h>
#include <TVector.h>
using namespace std;


//float, 
Double_t V0_Lifetime(TVector3 pv, TVector3 sv, TMatrixD EPV, TMatrixD ESV, Double_t M, TVector3 pT, double &ct, double &ect)
{
  // NOTA1: Esta funcion calcula el tiempo de vida y su error usando la longitud propia de decaimiento transversal
  // Recuerde que estamos asumiendo que el error en el pt es despreciable por eso las matrices asociadas a este las definimos cero

  TVector3 svT(sv.X(),sv.Y(),0.0);
  TVector3 pvT(pv.X(),pv.Y(),0.0);
  TVector3 d = svT - pvT;

  TMatrixD VSV(2,2);
  VSV(0,0) = ESV(0,0);
  VSV(1,1) = ESV(1,1);
  VSV(1,0) = ESV(1,0);
  VSV(0,1) = VSV(1,0);

  TMatrixD VPV(2,2);
  VPV(0,0) = EPV(0,0);
  VPV(1,1) = EPV(1,1);
  VPV(1,0) = EPV(1,0);
  VPV(0,1) = VPV(1,0);

  TMatrixD VL(2,2); VL = VSV; VL+=VPV;

  TVector3 p = pT;

  TMatrixD VP(2,2);
  VP(0,0) = 0.0;
  VP(1,1) = 0.0;
  VP(0,1) = 0.0;
  VP(1,0) = 0.0;

  double Lxy = d.Dot(p)/p.Mag();
  double lf = Lxy*M/p.Mag();
  //cout<<" ---> "<<lf<<endl;
  ct = lf;
  
  //Ahora calaculamos el error en el tiempo de vida
  
  //computing Mass error
  //double sM2 = 0; //We assume 0 for now
  
  //computing Lxy error
  
  //Defining Matrix:
  TMatrixD A(2,2);
  TMatrixD B(2,2);
  TMatrixD C(2,2);
  TMatrixD EP(2,2);
  TMatrixD EL(2,2);
  
  //Aij = PiPj/p2
  //Bij = LiLj/Lxy2 (Li = SVi - PVi)
  //EPij = Vij(P)/p2
  //ELij = Vij(L)/Lxy^2;
  //Cij = LiPj/(pLxy)
  
  A(0,0) = p.X()*p.X()/p.Mag2();
  A(1,1) = p.Y()*p.Y()/p.Mag2();
  A(0,1) = p.X()*p.Y()/p.Mag2();
  A(1,0) = A(0,1);

  B(0,0) = d.X()*d.X()/(Lxy*Lxy);
  B(1,1) = d.Y()*d.Y()/(Lxy*Lxy);
  B(0,1) = d.X()*d.Y()/(Lxy*Lxy);
  B(1,0) = B(0,1);
  
  C(0,0) = d.X()*p.X()/(Lxy*p.Mag());
  C(1,1) = d.Y()*p.Y()/(Lxy*p.Mag());
  C(0,1) = d.X()*p.Y()/(Lxy*p.Mag());
  C(1,0) = d.Y()*p.X()/(Lxy*p.Mag());

  EP = VP;
  EP*= ((double)1.0/p.Mag2());
  EL = VL;
  EL*= ((double)1.0/(Lxy*Lxy));

  //Test
  //EL(0,1) = 0.0;
  //EL(1,0) = 0.0;

  //Calculando Sigma Lxy
  // Sigma Lxy^2 = Tr{A*EL + (B + 4*A - 4*C)*EP}
  // NOTA2: en nuestro caso basicamente es Sigma Lxy^2 = Tr{A*EL), dado que no consideramos el momentum P
  
  TMatrixD A1 = A;
  A1*=(double)4.0;
  A1+=B;
  TMatrixD C1 = C;
  C1*=(double)4.0;
  A1-=C1;
  
  TMatrixD A_EL(A,TMatrixD::kMult,EL);
  TMatrixD A1_EP(A1,TMatrixD::kMult,EP);
  TMatrixD SL = A_EL;SL+=A1_EP;
  double sLxy2 = SL(0,0) + SL(1,1); 
  
  return ect = (double) fabs(lf)*sqrt(sLxy2);
}

void Make_Bujk_AODHI_2016(int era, int vq)
{
  TChain * ch = new TChain("ntuple","");
  //ch->Add("DATA/dataCpPb/0000/Rootuple_Bu_JpsiK_2016HI-AOD_534.root/rootuple/ntuple");// test
  
  if (era==1){
    ch->Add("/cms/jhovanny/mytest/BuHI/DATA/dataCpPb/*.root/rootuple/ntuple");
    
  }
  else{
    ch->Add("/cms/jhovanny/mytest/BuHI/DATA/dataCPbp/*.root/rootuple/ntuple");  
  }
  

  // ****************************************************************************************
       
  TTree *tree = (TTree*)ch;
  classntuple t(tree);
  Long64_t nentries = t.fChain->GetEntries();
  cout<<" Entries : "<<nentries<<endl;
  cout<<" Reading data ..."<<nentries<<endl;

  //TString outfileName( "ROOTS_Bujk_2018_AOD_1_best1.root" );
  TString outfileName(Form("ROOTS_Bujk_AOD_HI2016_Era%1i_best%1i.root",era,vq) );
 
  // **********************  **********************************                                   
  TFile* outputFile = TFile::Open( outfileName, "RECREATE" );
  
  //*******************************
  //*     VARIABLES TO SAVE       *
  //*******************************
   
  Double_t massB, massJ;
  Double_t Bdl,BdlIP;
  Double_t BdlE,BdlIPE; 
  Double_t Brapidity, Bphi;
  Double_t Jrapidity, Jphi;
  
  Double_t sigLxyJ, cosalfaJ;
  Double_t sigLxyB, cosalfaB;
  
  Double_t Jpsipt, Bpt;
  Double_t maxPtMu;
  Double_t minPtMu;
  
  Double_t mu1pt,  mu2pt;
  Double_t mu1phi,  mu2phi;
  Double_t mu1eta,  mu2eta;
  Int_t    mu1soft, mu2soft;
  Int_t muon1StationTight, muon2StationTight;
  
  Double_t pi1pt, pi1eta, pi1phi;
  Int_t    pi1Thits, pi1Phits, pi1q; 
  
  Double_t Bpro, Jpro;
  
  Double_t dxypi1, dxyEpi1;
  Double_t DRmumu;
  
  UInt_t          nPV;
  Int_t Ntrk, Ntrkvip, NtrkvipQ; 
  Int_t runn, evtn, Lumiblock, Kcharge;
  Int_t Tripal1open, Tripal1, Tripal2, Tripal3;
  UInt_t Trigger;
  
  TTree *treeS =  new TTree("treeS","signal");
  
  treeS->Branch("massB",   &massB,    "massB/D");
  treeS->Branch("Bpt",     &Bpt,      "Bpt/D");
  treeS->Branch("Brapidity",    &Brapidity,     "Brapidity/D");
  treeS->Branch("Bphi",    &Bphi,     "Bphi/D");
  treeS->Branch("Bdl",     &Bdl,      "Bdl/D");
  treeS->Branch("BdlE",     &BdlE,      "BdlE/D"); 
  treeS->Branch("BdlIP",     &BdlIP,      "BdlIP/D");
  treeS->Branch("BdlIPE",     &BdlIPE,      "BdlIPE/D"); 
  
  treeS->Branch("massJ",   &massJ,    "massJ/D");
  treeS->Branch("Jpsipt",  &Jpsipt,   "Jpsipt/D");
  treeS->Branch("Jrapidity",    &Jrapidity,     "Jrapidity/D");
  treeS->Branch("Jphi",    &Jphi,     "Jphi/D");
  treeS->Branch("sigLxyJ", &sigLxyJ,  "sigLxyJ/D");
  treeS->Branch("cosalfaJ", &cosalfaJ,  "cosalfaJ/D");
  treeS->Branch("sigLxyB", &sigLxyB,  "sigLxyB/D");
  treeS->Branch("cosalfaB", &cosalfaB,  "cosalfaB/D");
  
  treeS->Branch("Kcharge", &Kcharge, "Kcharge/I");   
  treeS->Branch("pi1pt",    &pi1pt,     "pi1pt/D");
  treeS->Branch("pi1eta",    &pi1eta,      "pi1eta/D");
  treeS->Branch("pi1phi",    &pi1phi,       "pi1phi/D");
  treeS->Branch("pi1Thits", &pi1Thits, "pi1Thits/I");
  treeS->Branch("pi1Phits", &pi1Phits, "pi1Phits/I");
  treeS->Branch("pi1q",    &pi1q,     "pi1q/I");
  
  treeS->Branch("mu1pt",   &mu1pt,     "mu1pt/D");
  treeS->Branch("mu2pt",   &mu2pt,     "mu2pt/D");  
  treeS->Branch("mu1phi",  &mu1phi,    "mu1phi/D");
  treeS->Branch("mu2phi",  &mu2phi,    "mu2phi/D");
  treeS->Branch("mu1eta",  &mu1eta,    "mu1eta/D");
  treeS->Branch("mu2eta",  &mu2eta,    "mu2eta/D");
  treeS->Branch("mu1soft",   &mu1soft,     "mu1soft/I");
  treeS->Branch("mu2soft",   &mu2soft,     "mu2soft/I");
  treeS->Branch("muon1StationTight",   &muon1StationTight,     "muon1StationTight/I");
  treeS->Branch("muon2StationTight",   &muon2StationTight,     "muon2StationTight/I");

  treeS->Branch("maxPtMu", &maxPtMu,   "maxPtMu/D");    
  treeS->Branch("minPtMu", &minPtMu,   "minPtMu/D");
  
  treeS->Branch("Bpro",   &Bpro,    "Bpro/D");
  treeS->Branch("Jpro",   &Jpro,    "Jpro/D");
  
  treeS->Branch("dxypi1", &dxypi1,  "dxypi1/D");
  treeS->Branch("dxyEpi1", &dxyEpi1,  "dxyEpi1/D");
  treeS->Branch("DRmumu", &DRmumu,   "DRmumu/D");

  treeS->Branch("Ntrk", &Ntrk, "Ntrk/I");
  treeS->Branch("Ntrkvip", &Ntrkvip, "Ntrkvip/I");
  treeS->Branch("NtrkvipQ", &NtrkvipQ, "NtrkvipQ/I");
  treeS->Branch("nPV", &nPV, "nPV/i");
  
  treeS->Branch("runn", &runn, "runn/I");
  treeS->Branch("evtn", &evtn, "evtn/I");
  treeS->Branch("Lumiblock", &Lumiblock, "Lumiblock/I");

  treeS->Branch("Tripal1open", &Tripal1open, "Tripal1open/I");
  treeS->Branch("Tripal1", &Tripal1, "Tripal1/I");
  treeS->Branch("Tripal2", &Tripal2, "Tripal2/I");
  treeS->Branch("Tripal3", &Tripal3, "Tripal3/I");
  treeS->Branch("Trigger", &Trigger, "Trigger/i");

  
  TVector3 pT,pv,sv,pvbs,pvrf,bsv,pvip,pvipBSc;
  Double_t ct,ect,M,ctbs,ectbs,ctrf,ectrf,ctip,ectip,ctipBSc,ectipBSc,ctBS,ctBSE;
  TMatrix ESV(3,3);
  TMatrix EPV(3,3);
  TMatrix EPVBS(3,3);
  TMatrix EPVRf(3,3);
  TMatrix EBSV(3,3);
  TMatrix EPVip(3,3);
  TMatrix EPVipBSc(3,3);
  
  Int_t pt_ip_igual=0;
  Int_t pt_ip_diferente=0;
  
  Int_t nTen = nentries/10;

  Int_t nbytes = 0, nb = 0;
  for(Long64_t jentry=0; jentry<nentries;jentry++) 
    {
      Long64_t ientry = t.LoadTree(jentry);
      if (ientry < 0) break;
      nb = t.fChain->GetEntry(jentry);   nbytes += nb;
      if(jentry%nTen==0) cout<<10*(jentry/nTen)<<"%-"<<flush;
      if(jentry==nentries-1) cout<<endl;
      
      pv.SetXYZ(t.priVtxX,t.priVtxY,t.priVtxZ);
      EPV(0,0) = t.priVtxXE;
      EPV(1,1) = t.priVtxYE;
      EPV(2,2) = t.priVtxZE;
      EPV(0,1) = t.priVtxXYE;
      EPV(0,2) = t.priVtxXZE;
      EPV(1,2) = t.priVtxYZE;
      
      
      // para escoger el mejor Chi2 prob vertex     vq=1            
      Double_t proval=-10.0;                                                                                                                                                           
      Int_t idbest=-1;                                                                                                                                                                     
      for(unsigned int j=0;j<t.B_J_mass->size();j++)                                                                                                                                   
	{                                                                                                                                                                                
	  if(t.B_Prob->at(j)>proval)                                                                                                                                                       
	    {                                                                                                                                                                            
	      idbest = j;                                                                                                                                                            
	      proval = t.B_Prob->at(j);                                                                                                                                              
	    }                                                                                                                                                                            
	}
      
      
      
      // para escoger el mejor Pt candidate    vq=2
      /*
      Double_t ptval=-10.0;                                                                                                                                                           
      Int_t idbest=-1;                                                                                                                                                                     
      for(unsigned int j=0;j<t.B_J_mass->size();j++)                                                                                                                                   
	{               
	  
	  TVector3 Bcpt(t.B_px->at(j),t.B_py->at(j),t.B_pz->at(j));
	  
	  if(Bcpt.Pt()>ptval)                                                                                                                                                       
	    {                                                                                                                                                                            
	      idbest = j;                                                                                                                                                            
	      ptval = Bcpt.Pt();                                                                                                                                              
	    }                                                                                                                                                                            
	} 
      */
      
      for(unsigned int i=0;i<t.B_J_mass->size();i++)
	{
	  
	  if(i!=idbest) continue; // this is to chose the best canidate per even                                                                                                       
	  
	  TVector3 B(t.B_px->at(i),t.B_py->at(i),t.B_pz->at(i));	   
	  TVector3 Jpsi( t.B_J_px->at(i), t.B_J_py->at(i), t.B_J_pz->at(i) );
	  
	  TVector3 pi( t.B_k_px_track->at(i), t.B_k_py_track->at(i), t.B_k_pz_track->at(i) );
	  TVector3 mu1( t.B_J_px1->at(i), t.B_J_py1->at(i), t.B_J_pz1->at(i) );
	  TVector3 mu2( t.B_J_px2->at(i), t.B_J_py2->at(i), t.B_J_pz2->at(i) );	  	     
	  
	  pT.SetXYZ(t.B_px->at(i),t.B_py->at(i),0.0);
	  
	  sv.SetXYZ(t.B_DecayVtxX->at(i),t.B_DecayVtxY->at(i),t.B_DecayVtxZ->at(i));
	  
	  ESV(0,0) = t.B_DecayVtxXE->at(i);
	  ESV(1,1) = t.B_DecayVtxYE->at(i);
	  ESV(2,2) = t.B_DecayVtxZE->at(i);
	  ESV(0,1) = t.B_DecayVtxXYE->at(i);
	  ESV(0,2) = t.B_DecayVtxXZE->at(i);
	  ESV(1,2) = t.B_DecayVtxYZE->at(i);

	  pvip.SetXYZ(t.pVtxIPX->at(i),t.pVtxIPY->at(i),t.pVtxIPZ->at(i));
	  EPVip(0,0) = t.pVtxIPXE->at(i);
	  EPVip(1,1) = t.pVtxIPYE->at(i);
	  EPVip(2,2) = t.pVtxIPZE->at(i);
	  EPVip(0,1) = t.pVtxIPXYE->at(i);
	  EPVip(0,2) = t.pVtxIPXZE->at(i);
	  EPVip(1,2) = t.pVtxIPYZE->at(i);

	  V0_Lifetime(pv,sv,EPV,ESV, 5.27932, pT, ct, ect);
	  V0_Lifetime(pvip,sv,EPVip,ESV, 5.27932, pT, ctip, ectip);

	  //*******************
	  // Jpsi(mumu) Cuts  *
	  //*******************
	  
	  Double_t dxJ = t.B_J_DecayVtxX->at(i) - t.pVtxIPX->at(i);
	  Double_t dyJ = t.B_J_DecayVtxY->at(i) - t.pVtxIPY->at(i);
	  
	  Double_t dxJE = t.B_J_DecayVtxXE->at(i);// + t.priVtxXE;
	  Double_t dyJE = t.B_J_DecayVtxYE->at(i);// + t.priVtxYE;
	  
	  Double_t sigLxyJtmp =(dxJ*dxJ +dyJ*dyJ)/sqrt( dxJ*dxJ*dxJE*dxJE + dyJ*dyJ*dyJE*dyJE );
	  
	  if( t.B_J_mass->at(i)<2.9 || t.B_J_mass->at(i)>3.3 ) continue;
	  if(fabs(mu1.Eta())>2.4 || fabs(mu2.Eta())>2.4) continue;
	  if( mu1.Pt()<1.5 || mu2.Pt()<1.5)continue;
	  //if( Jpsi.Pt()<2.0) continue;
	  Double_t cosAlphaXY = ( t.B_J_px->at(i)*dxJ + t.B_J_py->at(i)*dyJ )/( sqrt(dxJ*dxJ+dyJ*dyJ)*Jpsi.Pt() );
	  //if(cosAlphaXY<0.9) continue;
	  //if( sigLxyJtmp < 3.0 )continue;
	  
	  //*  Lxy and cosalpha for B  *
	  Double_t dxB = t.B_DecayVtxX->at(i) - t.priRfVtxX->at(i);
	  Double_t dyB = t.B_DecayVtxY->at(i) - t.priRfVtxY->at(i);
	  Double_t dxBE = t.B_DecayVtxXE->at(i);
	  Double_t dyBE = t.B_DecayVtxYE->at(i);
	  
	  Double_t sigLxyBtmp =(dxB*dxB +dyB*dyB)/sqrt( dxB*dxB*dxBE*dxBE + dyB*dyB*dyBE*dyBE );
	  Double_t cosAlphaBXY = ( t.B_px->at(i)*dxB + t.B_py->at(i)*dyB )/( sqrt(dxB*dxB+dyB*dyB)*B.Pt() );
	  
	  //*******************
	  //*   Final Cuts    *
	  //*******************
	  if((t.B_mass->at(i)>6.0) || (t.B_mass->at(i)<5.0))  continue;
	  if( B.Pt()<=1.0) continue;
	  //if(ctip<0.0) continue;
	  //if( pi.Pt()<0.5) continue;
	  if(fabs(pi.Eta())>2.4 || fabs(pi.Eta())>2.4) continue;
	  
	  //*******************
	  //*   fill tree    *
	  //*******************
	  
	  // for Rapidity variables
	  TLorentzVector J_p4, B_p4, B2_p4;
	  J_p4.SetPtEtaPhiM(Jpsi.Pt(), Jpsi.Eta(), Jpsi.Phi(), t.B_J_mass->at(i));
	  B_p4.SetPtEtaPhiM(B.Pt(), B.Eta(), B.Phi(), t.B_mass->at(i) );
	  
	  runn= t.run;
	  evtn =t.event;
	  Lumiblock = t.lumiblock;

	  //Tripal1open = t.tri_PAL1DoubleMu0_open->at(i);
	  //Tripal1 = t.tri_PAL1DoubleMu0->at(i);
	  Tripal2 = t.tri_PAL2DoubleMu0->at(i);
	  //Tripal3 = t.tri_PAL3DoubleMu0->at(i);
	  Trigger = t.trigger ;
	  
	  massB   = t.B_mass->at(i);
	  massJ   = t.B_J_mass->at(i);
	  sigLxyJ = sigLxyJtmp;
	  cosalfaJ= cosAlphaXY;
	  sigLxyB = sigLxyBtmp;
	  cosalfaB= cosAlphaBXY;
	  
	  Bdl   = ct;
	  BdlIP = ctip;
	  BdlE  =  ect;
	  BdlIPE = ectip;
	  
	  //Beta    = B.Eta();
	  Brapidity    = B_p4.Rapidity();
	  Bphi    = B.Phi();
	  Bpt     = B.Pt();
	  Jpsipt  = Jpsi.Pt();
	  //Jeta    = Jpsi.Eta();
	  Jrapidity    = J_p4.Rapidity();
	  Jphi    = Jpsi.Phi();
	  
	  pi1pt    = pi.Pt();
	  pi1phi   = pi.Phi();
	  pi1eta   = pi.Eta();
	  pi1Thits = t.pi1_trackerhits->at(i);
	  pi1Phits = t.pi1_pixelhits->at(i);
	  pi1q     = t.pi1Q->at(i);
	  
	  mu1pt   = mu1.Pt();
	  mu2pt   = mu2.Pt();	 
	  mu1phi  = mu1.Phi();
	  mu2phi  = mu2.Phi();         
	  mu1eta  = mu1.Eta();
	  mu2eta  = mu2.Eta();
	  mu1soft = t.muon1Soft->at(i);
	  mu2soft = t.muon2Soft->at(i);
	  muon1StationTight = t.muon1StationTight->at(i);
	  muon2StationTight = t.muon2StationTight->at(i);

	  Bpro  = t.B_Prob->at(i);
	  Jpro  = t.B_J_Prob->at(i);
	  Kcharge = t.B_k_charge1->at(i);
	  
	  dxypi1 = t.pi1dxy->at(i);
	  dxyEpi1 = t.pi1dxyE->at(i);
	  DRmumu   = mu1.DeltaR(mu2);
	  nPV = t.nVtx;
	  Ntrk = t.nTrk;
	  Ntrkvip = t.nTrk_VtxIP->at(i);
	  NtrkvipQ = t.nTrk_VtxIP_Q->at(i);

	  if( mu1.Pt() > mu2.Pt() ){
	    maxPtMu = mu1.Pt();
	    minPtMu = mu2.Pt();
	  }
	  else{
	    maxPtMu = mu2.Pt();
	    minPtMu = mu1.Pt();
	  }
	  
	  
	  treeS->Fill();
	  
	  
	}
      
       
    }
  
  
  outputFile->Write("",TObject::kOverwrite);

}
