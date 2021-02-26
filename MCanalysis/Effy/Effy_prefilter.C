#include "classmcrootuple.C"

#include <TROOT.h>
#include "TMath.h"
#include <iostream>
#include <fstream>

#include "RooPolyVar.h"
#include "RooConstVar.h"
#include "RooChi2Var.h"
#include "RooMinuit.h"

#include "RooFit.h"
#include "RooFitResult.h"
#include "RooGaussian.h"
#include "RooExponential.h"
#include "RooPolynomial.h"
#include "RooChebychev.h"
#include "RooAddPdf.h"

#include "TPaveText.h"
#include "TLegend.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TGaxis.h"
#include "TAxis.h"
#include "TLatex.h"
#include "RooPlot.h"
#include "RooHist.h"
#include "TF1.h"
#include "RooNumIntConfig.h"

#include "RooDataSet.h"
#include "RooRealVar.h"

using namespace std;
using namespace RooFit ;

void Effy_prefilter() {

 gStyle->SetOptTitle(0);
 gStyle->SetOptFit(0);
 gStyle->SetOptStat(0); 
 
 TChain * ch = new TChain("ntuple","");
 ch->Add("../dataMC/Rootuple_Bu_JpsiK_2016HI-GenLevelbacht1.root/rootuple/ntuple");

 TTree *tree = (TTree*)ch;
 classmcrootuple t(tree);
 Long64_t nentries = t.fChain->GetEntries();
 cout<<" Entries : "<<nentries<<endl;
 
 cout<<"xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"<<endl;
 cout<<" Reading data ..."<<nentries<<endl;

 //define los arreglos to get the prefilter Efficiency
 const Int_t bin = 4;
 Double_t Effp0[bin];
 Double_t Effp1[bin];
 Int_t EpncT = 0; Double_t EpwcT = 0;
 
 Int_t nTen = nentries/10;
 Int_t nbytes = 0, nb = 0;
 for(Long64_t jentry=0; jentry<nentries;jentry++)
 //for(Long64_t jentry=0; jentry<20000;jentry++)   
   {
     Long64_t ientry = t.LoadTree(jentry);
     if (ientry < 0) break;
     nb = t.fChain->GetEntry(jentry);   nbytes += nb;
     if(jentry%nTen==0) cout<<10*(jentry/nTen)<<"%-"<<flush;
     if(jentry==nentries-1) cout<<endl;
     
     if( (TMath::Abs(t.gen_bc_p4->Rapidity()) > 2.4) )continue;
     if( (t.gen_bc_p4->Pt()<7.0) || (t.gen_bc_p4->Pt()>=50.0) )continue;
     EpncT = EpncT + 1;
     
     //"prefilter Effyciencia  no cuts"
     //if(t.gen_bc_p4->Pt()>=12.0 && t.gen_bc_p4->Pt()<13.0){Effp0[0]++;}
     if(t.gen_bc_p4->Pt()>=7.0 && t.gen_bc_p4->Pt()<10.0){Effp0[0]++;}
     if(t.gen_bc_p4->Pt()>=10.0 && t.gen_bc_p4->Pt()<15.0){Effp0[1]++;}
     if(t.gen_bc_p4->Pt()>=15.0 && t.gen_bc_p4->Pt()<20.0){Effp0[2]++;}
     if(t.gen_bc_p4->Pt()>=20.0 && t.gen_bc_p4->Pt()<50.0){Effp0[3]++;}
 
     //"prefilter Effyciencia  with cuts"
     if(t.gen_muon1_p4->Pt()>1.5 && t.gen_muon2_p4->Pt()>1.5 &&
	(TMath::Abs(t.gen_muon1_p4->Eta()) < 2.5) && (TMath::Abs(t.gen_muon2_p4->Eta()) < 2.5) && t.gen_pion3_p4->Pt()>0.4 && (TMath::Abs(t.gen_pion3_p4->Eta()) < 2.5)){
       
       if(t.gen_bc_p4->Pt()>=7.0 && t.gen_bc_p4->Pt()<10.0){Effp1[0]++;}
       if(t.gen_bc_p4->Pt()>=10.0 && t.gen_bc_p4->Pt()<15.0){Effp1[1]++;}
       if(t.gen_bc_p4->Pt()>=15.0 && t.gen_bc_p4->Pt()<20.0){Effp1[2]++;}
       if(t.gen_bc_p4->Pt()>=20.0 && t.gen_bc_p4->Pt()<50.0){Effp1[3]++;}
             
       if(t.gen_bc_p4->Pt()>=7.0 && t.gen_bc_p4->Pt()<50.0)
	  {
	    EpwcT++;
	  }
	
     }
     
   }
 //cout << Epwc0 << endl;
 cout<< "entries for pt>7 y pt<10. No cuts: "<< Effp0[0] << " with cuts:  " << Effp1[0] << endl; // =
 cout<< "entries for pt>7 y pt<50. No cuts: "<< EpncT << " with cuts:  " << EpwcT << endl; // = 
 cout<< "prefilter effy pt>7 y pt<50: "<< EpwcT/EpncT << " +/-:  " << (EpwcT/EpncT)*sqrt( ((EpwcT+1.0)*(EpncT-EpwcT+1.0))/((EpncT+3.0)*(EpncT+2.0)*(EpncT+2.0)) )  << endl;
 cout << endl;
 ofstream salidaT ("plotseffyBu/output_Buprefiltereffy_pt_Integrada.txt");
 salidaT.is_open();
 salidaT << EpwcT/EpncT << " " << (EpwcT/EpncT)*sqrt( ((EpwcT+1.0)*(EpncT-EpwcT+1.0))/((EpncT+3.0)*(EpncT+2.0)*(EpncT+2.0)) ) ;
 salidaT.close();
 //return;

 //define los arreglos para el bineado
 Double_t Bpt[bin] = {8.5,12.5,17.5,35.0};
 Double_t Bpte[bin] = {1.5,2.5,2.5,15.0};

 //define los arreglos de prefilter para usar en la eficencia total 
 Double_t Bpre[bin] ;
 Double_t BpreE[bin];
 
 RooRealVar x_pre("x_pre","x_pre",7.0,50) ;
 RooRealVar y_pre("y_pre","y_pre",0,200) ;
 RooDataSet dxy_pre("dxy_pre","dxy_pre",RooArgSet(x_pre,y_pre),StoreError(RooArgSet(x_pre,y_pre))) ;     
 for(int i=0; i<bin; i++){
   x_pre = Bpt[i];
   x_pre.setError(Bpte[i]) ;
   
   y_pre=(Effp1[i]/Effp0[i]);
   //Double_t dive_pre = (Effp1[i]/Effp0[i])*sqrt( (1/Effp1[i]) + (1/Effp0[i]) ) ;// usando el error como la raiz cuadrada del numero
   Double_t dive_pre =  sqrt( ((Effp1[i]+1.0)*(Effp0[i]-Effp1[i]+1.0))/((Effp0[i]+3.0)*(Effp0[i]+2.0)*(Effp0[i]+2.0)) ); // formula by d'Agostini
   //Double_t dive_pre =  sqrt( Effp1[i]*(1-(Effp1[i]/Effp0[i])) )/Effp0[i] ;//Binomial model
   
   y_pre.setError(dive_pre);  
   dxy_pre.add(RooArgSet(x_pre,y_pre)) ;
   
   cout<<"prefilter efficiency for pt bins:  "<<Effp1[i]/Effp0[i]<<"  +/- "<<dive_pre<<endl;  
   Bpre[i] = (Effp1[i]/Effp0[i]);
   BpreE[i] = dive_pre;     
 }
 dxy_pre.Print("v");

 ofstream salida ("plotseffyBu/output_Buprefiltereffy_pt.txt");
 salida.is_open();
 //salida << Bpre[0] << " " << BpreE[0] << " " << Bpre[1] << " " <<  BpreE[1] << " " << Bpre[2] << " " <<  BpreE[2] ;
 for(int ii=0; ii<bin; ii++){
   salida << Bpre[ii] << " " << BpreE[ii] << " ";
 }
 salida.close();

 int H = 600;
 int W = 800;
 
 TCanvas *c1 = new TCanvas("c1","c1",50,50,W,H);
 //c1->Divide(1,2);
 c1->SetLeftMargin(0.125);
 c1->SetRightMargin(0.01);
 c1->SetTopMargin(0.09);
 c1->SetBottomMargin(0.14); 
 
 RooPlot* frame_pre = x_pre.frame(5.0,55.0,bin);
 dxy_pre.plotOnXY(frame_pre,YVar(y_pre)) ;

 frame_pre->SetYTitle("pre-filter efficiency");
 //frame_pre->SetYTitle("efficiency"); 
 frame_pre->SetXTitle("p_{T}(J/#psi K^{+}) (GeV)");                                                                                                                    
 frame_pre->SetLabelSize(0.05,"XY");
 frame_pre->SetTitleSize(0.05,"XY");
 frame_pre->GetYaxis()->CenterTitle();   
 frame_pre->GetXaxis()->CenterTitle();
 frame_pre->GetYaxis()->SetNdivisions(505,1);
 //frame_pre->GetXaxis()->SetDecimals(2); 
 frame_pre->GetXaxis()->SetNdivisions(504,1); 
 frame_pre->SetTitleOffset(1.0,"X");
 frame_pre->SetTitleOffset(0.9,"Y");
 frame_pre->SetTitleSize(0.06,"XY");
 //frame_pre->SetMinimum(0.2);
 //frame_pre->SetMaximum(2200.0);
 frame_pre->Draw();

 TLatex *   tex1 = new TLatex(0.94,0.926,"MC simulation");
 tex1->SetNDC();
 tex1->SetTextAlign(31);
 tex1->SetTextFont(42);
 tex1->SetTextSize(0.04); 
 tex1->SetLineWidth(2);
 
 TLatex *tex2 = new TLatex(0.17,0.926,"CMS");
 tex2->SetNDC();
 tex2->SetTextFont(61);
 tex2->SetTextSize(0.04); 
 tex2->SetLineWidth(2);

 TLatex *tex3 = new TLatex(0.24,0.926,"Preliminary");
 tex3->SetNDC();
 tex3->SetTextFont(52);
 tex3->SetTextSize(0.04); 
 tex3->SetLineWidth(2);

 tex1->Draw();  
 tex2->Draw();
 tex3->Draw();

 c1->Modified();
 
 c1->Print("plotseffyBu/Bueffy_prefilter.png");
 c1->Print("plotseffyBu/Bueffy_prefilter.pdf");
 
   
}
