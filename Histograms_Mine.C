#include "classreduce.C"

#include <TROOT.h>
#include "TMath.h"
#include <iostream>
#include <fstream>

#include "RooFit.h"
#include "RooFitResult.h"
#include "RooGlobalFunc.h"
#include "RooArgusBG.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooCBShape.h"
#include "RooExponential.h"
#include "RooGenericPdf.h"
#include "RooFFTConvPdf.h"
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
#include "RooNumIntConfig.h"

using namespace RooFit;
using namespace std;

int Nhisto;

TCanvas* CreateCanvas1D(TString cname, TH1F *h1, TString yhname, TString xhname, TString Lum) 
{

 int H = 600;
 int W = 800;

 TCanvas* canv = new TCanvas(cname,cname,50,50,W,H);
 canv->cd();
 canv->SetLeftMargin(0.12);
 canv->SetRightMargin(0.04);
 canv->SetTopMargin(0.07);
 canv->SetBottomMargin(0.14);
 canv->SetGrid();
 //gPad->SetLogy();

 h1->Draw("e1");
 h1->SetMarkerStyle(24);
 h1->SetMarkerSize(1.0);
 h1->SetMarkerColor(4);
 h1->SetLineColor(4);
 h1->SetLineWidth(2);
 h1->GetYaxis()->CenterTitle(true);
 //h1->SetYTitle("Pt(J/#psi) [GeV]");
 h1->SetYTitle(yhname);
 h1->GetXaxis()->CenterTitle(true); 
 //h1->SetXTitle("|#eta(J/#psi )|"); //#eta(#pi_{Bu}
 h1->SetXTitle(xhname); //#eta(#pi_{Bu}
 h1->SetTitleSize(35,"XY"); 
 h1->SetLabelSize(30,"XY");
 h1->SetTitleOffset(1.0,"Y");
 h1->SetTitleOffset(1.0,"X");
 h1->SetLabelFont(43,"XY");  
 h1->SetTitleFont(43,"XY");
 h1->SetMinimum(1.0);
 //h1->SetMinimum(1e-5);  
 //h1->SetMaximum(0.599);

 //TLatex *   tex1 = new TLatex(0.88,0.926,"36.1 fb^{-1} (13 TeV, 2016)");
 TLatex *   tex1 = new TLatex(0.96,0.94,Lum); 
 tex1->SetNDC();
 tex1->SetTextAlign(31);
 tex1->SetTextFont(42);
 tex1->SetTextSize(0.05); 
 tex1->SetLineWidth(2);
 
 TLatex *tex2 = new TLatex(0.115,0.94,"CMS");
 tex2->SetNDC();
 tex2->SetTextFont(61);
 tex2->SetTextSize(0.05); 
 tex2->SetLineWidth(2);

 TLatex *tex3 = new TLatex(0.21,0.94,"Preliminary");
 tex3->SetNDC();
 tex3->SetTextFont(52);
 tex3->SetTextSize(0.05); 
 tex3->SetLineWidth(2); 

 tex1->Draw();
 tex2->Draw();
 tex3->Draw();

 if(Nhisto==1) // Ntrk
   {
     TLine *l1 =new TLine(35.0,0.0,35.0,1416.0);
     l1->SetLineColor(kRed);
     l1->Draw();
     
     TLine *l2 =new TLine(75.0,0.0,75.0,1416.0);
     l2->SetLineColor(kRed);
     l2->Draw();

     TLine *l3 =new TLine(185.0,0.0,185.0,1416.0);
     l3->SetLineColor(kRed);
     l3->Draw();

     TLine *l4 =new TLine(250.0,0.0,250.0,1416.0);
     l4->SetLineColor(kRed);
     l4->Draw();
   }

  if(Nhisto==2) // Bupt
   {
     TLine *l1 =new TLine(10.0,0.0,10.0,2255.0);
     l1->SetLineColor(kRed);
     l1->Draw();
     
     TLine *l2 =new TLine(15.0,0.0,15.0,2255.0);
     l2->SetLineColor(kRed);
     l2->Draw();

     TLine *l3 =new TLine(20.0,0.0,20.0,2255.0);
     l3->SetLineColor(kRed);
     l3->Draw();

     TLine *l4 =new TLine(50.0,0.0,50.0,2255.0);
     l4->SetLineColor(kRed);
     l4->Draw();
   }

   if(Nhisto==3) // Brapidity
   {
     TLine *l1 =new TLine(-0.8,0.0,-0.8,247.0);
     l1->SetLineColor(kRed);
     l1->Draw();
     
     TLine *l2 =new TLine(0.8,0.0,0.8,247.0);
     l2->SetLineColor(kRed);
     l2->Draw();

     TLine *l3 =new TLine(-1.3,0.0,-1.3,247.0);
     l3->SetLineColor(kRed);
     l3->Draw();

     TLine *l4 =new TLine(1.3,0.0,1.3,247.0);
     l4->SetLineColor(kRed);
     l4->Draw();

     TLine *l5 =new TLine(-1.7,0.0,-1.7,247.0);
     l5->SetLineColor(kRed);
     l5->Draw();
     
     TLine *l6 =new TLine(1.7,0.0,1.7,247.0);
     l6->SetLineColor(kRed);
     l6->Draw();
   }
 
 return canv;
}



void Histograms_Mine(int era = 3, int vq = 1)
{
 gStyle->SetOptTitle(0);
 gStyle->SetOptFit(0);
 gStyle->SetOptStat(0);
 gStyle->SetErrorX(0);
 
 Double_t Mmin = 5.0; 
 Double_t Mmax = 5.6;
 TString Lumi = "";  
 
 TChain *ch = new TChain("butree","");
 if(era==1){ ch->Add(Form("reducetree_Bujk_AOD_HI2016_Era1_best%1i.root/butree",vq));  Lumi = "64.4 nb^{-1} (#it{p}Pb 8.16 TeV)";}  // 64.41
 else if(era==2){ch->Add(Form("reducetree_Bujk_AOD_HI2016_Era2_best%1i.root/butree",vq)); Lumi = "115.3 nb^{-1} (Pb#it{p} 8.16 TeV)";} //115.28
 else
   {
     ch->Add(Form("reducetree_Bujk_AOD_HI2016_Era1_best%1i.root/butree",vq));
     ch->Add(Form("reducetree_Bujk_AOD_HI2016_Era2_best%1i.root/butree",vq));
     Lumi = "189.7 nb^{-1} (#it{p}Pb 8.16 TeV)"; // 179.69
   }
 
 TTree *tree = (TTree*) ch;
 classreduce t(tree);
 Long64_t nentries = t.fChain->GetEntries();
 cout<<" Entries : "<<nentries<<endl;
 
 //------------------------------
 RooRealVar M("M"," M(J/#psi K^{+}) (GeV)",Mmin,Mmax);
 RooDataSet data("data","data",RooArgSet(M));


 Double_t Ntrmax = 450.0;
 Double_t Ntrmin = 0.0;
 Double_t Ntrbin = (Ntrmax - Ntrmin)/10.0;
 TH1F *histoNtr= new TH1F("histoNtr","histoNtr",Ntrbin,Ntrmin,Ntrmax);

 Double_t Ntrvipmax = 300.0;
 Double_t Ntrvipmin = 0.0;
 Double_t Ntrvipbin = (Ntrvipmax - Ntrvipmin)/10.0;
 TH1F *histoNtrvip= new TH1F("histoNtrvip","histoNtrvip",Ntrvipbin,Ntrvipmin,Ntrvipmax);
 TH1F *histoNtrvipQ= new TH1F("histoNtrvipQ","histoNtrvipQ",Ntrvipbin,Ntrvipmin,Ntrvipmax);

 Double_t Buptmax = 70.0;
 Double_t Buptmin = 7.0;
 Double_t Buptbin = (Buptmax - Buptmin)/2.0;
 TH1F *histoBupt = new TH1F("histoBupt","histoBupt",Buptbin,Buptmin,Buptmax);

 /*
 Double_t BurapidityRmax = 2.4;
 Double_t BurapidityRmin = 0.0;
 Double_t BurapidityRbin = (BurapidityRmax - BurapidityRmin)/0.05;
 TH1F *histoBurapidityR = new TH1F("histoBurapidityR","histoBurapidityR",BurapidityRbin,BurapidityRmin,BurapidityRmax);
 
 Double_t BurapidityLmax = 0.0;
 Double_t BurapidityLmin = -2.4;
 Double_t BurapidityLbin = (BurapidityLmax - BurapidityLmin)/0.05;
 TH1F *histoBurapidityL = new TH1F("histoBurapidityL","histoBurapidityL",BurapidityLbin,BurapidityLmin,BurapidityLmax);
 */
 
 Double_t Burapiditymax = 2.4;
 Double_t Burapiditymin = -2.4;
 Double_t Burapiditybin = (Burapiditymax)*2/0.05;
 TH1F *histoBurapidity = new TH1F("histoBurapidity","histoBurapidity",Burapiditybin,Burapiditymin,Burapiditymax);

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
     
     //Mass windows cuts
     if(t.B_mass<=Mmin || t.B_mass>=Mmax) continue;	
     if(t.J_mass<=2.9 || t.J_mass>=3.3) continue;
     //if(t.J_mass<=3.0969-0.150 || t.J_mass>=3.0969+0.150) continue;
     
     if(t.mu1pt<2.0 || t.mu2pt<2.0) continue;
     
     if(t.Bupt<7.0 || t.Bupt>=70.0)continue; 
     if(abs(t.rapidityB)>2.4)continue;
     //if((t.pdl/t.pdle)<5.0)continue;
     
     //if(t.pion1pt<0.5)continue;     
     if(t.pion1pt<0.8)continue;   
  
     M=t.B_mass;    
     data.add(RooArgSet(M));
     
     histoNtr->Fill(t.Ntrk);
     histoNtrvip->Fill(t.Ntrkvip);
     histoNtrvipQ->Fill(t.NtrkvipQ);
     histoBupt->Fill(t.Bupt); 
     histoBurapidity->Fill(t.rapidityB);
     
     

     /*
     if( t.mu1pt > t.mu2pt )
       {
	 histoptmul->Fill(t.mu1pt);
	 histoptmusl->Fill(t.mu2pt);
       }
     else
       {
	 histoptmul->Fill(t.mu2pt);
	 histoptmusl->Fill(t.mu1pt);      
       }
     */
   }


 
cout<<"xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"<<endl;
data.Print("v");
//return;


 
TCanvas* canv_histoNtr = CreateCanvas1D("canv_histoNtr", histoNtr, "Events / 10", "Number of charged  particle  tracks", Lumi); 
canv_histoNtr->SaveAs("plotsHistograms/Histo_histoNtr.png");
canv_histoNtr->SaveAs("plotsHistograms/Histo_histoNtr.pdf");


TCanvas* canv_histoNtrvip = CreateCanvas1D("canv_histoNtrvip", histoNtrvip, "Events / 10", "Tracks originate at the primary vertex", Lumi); 
canv_histoNtrvip->SaveAs("plotsHistograms/Histo_histoNtrvip.png");
canv_histoNtrvip->SaveAs("plotsHistograms/Histo_histoNtrvip.pdf");


Nhisto=1; 
TCanvas* canv_histoNtrvipQ = CreateCanvas1D("canv_histoNtrvipQ", histoNtrvipQ, "Events / 10", "Tracks originate at the primary vertex", Lumi); 
canv_histoNtrvipQ->SaveAs("plotsHistograms/Histo_histoNtrvipQ.png");
canv_histoNtrvipQ->SaveAs("plotsHistograms/Histo_histoNtrvipQ.pdf"); 

Nhisto=2;
TCanvas* canv_histoBupt = CreateCanvas1D("canv_histoBupt", histoBupt, "Events / 10", "B^{+} transverse momentum", Lumi); 
canv_histoBupt->SaveAs("plotsHistograms/Histo_histoBupt.png");
canv_histoBupt->SaveAs("plotsHistograms/Histo_histoBupt.pdf");

/*
TCanvas* canv_histoBurapidityR = CreateCanvas1D("canv_histoBurapidityR", histoBurapidityR, "Events / 10", "B^{+} rapidity (right)", Lumi); 
canv_histoBurapidityR->SaveAs("plotsHistograms/Histo_histoBurapidityR.png");
canv_histoBurapidityR->SaveAs("plotsHistograms/Histo_histoBurapidityR.pdf");
 
TCanvas* canv_histoBurapidityL = CreateCanvas1D("canv_histoBurapidityL", histoBurapidityL, "Events / 10", "B^{+} rapidity (left)", Lumi); 
canv_histoBurapidityL->SaveAs("plotsHistograms/Histo_histoBurapidityL.png");
canv_histoBurapidityL->SaveAs("plotsHistograms/Histo_histoBurapidityL.pdf");
*/

Nhisto=3;
TCanvas* canv_histoBurapidity = CreateCanvas1D("canv_histoBurapidity", histoBurapidity, "Events / 10", "B^{+} rapidity", Lumi); 
canv_histoBurapidity->SaveAs("plotsHistograms/Histo_histoBurapidity.png");
canv_histoBurapidity->SaveAs("plotsHistograms/Histo_histoBurapidity.pdf");

return;
RooAbsReal::defaultIntegratorConfig()->setEpsAbs(1e-10) ;
RooAbsReal::defaultIntegratorConfig()->setEpsRel(1e-10) ;

//---- Mass model -------
 Double_t supM = Mmax;
 Double_t infM = Mmin;
 
 // ****define nominal background ****
 RooRealVar a0("a0","a0",0.23,-10.0,10.0);
 RooRealVar a1("a1","a1",-0.04,-10.0,10.0);
 RooRealVar a2("a2","a2",0.2,-10.0,10.0);
 RooRealVar a3("a3","a3",0.2,-10.0,10.0);
 //RooChebychev bkg1("bkg1","Background",M,RooArgList(a0)); 
 //RooChebychev bkg1("bkg1","Background",M,RooArgList(a0,a1));
 //RooChebychev bkg1("bkg1","Background",M,RooArgList(a0,a1,a2));

 RooRealVar c("c","c",-2.0,-10.0,10.0);
 RooExponential bkg1("bkg1","Exp. Background",M,c);

 //**** define nominal signal ****
 //gaussians                                                                                                                                               
 //RooRealVar mean("mean"," Mass mean",3.686097);
 RooRealVar mean("mean"," Mass mean",5.279,5.250,5.3,"GeV");
 RooRealVar width("width"," Mass width",0.015,0.005,0.020,"GeV");
 RooGaussian Sig("Sig"," Signal PDF",M,mean,width);

 //RooRealVar meanB("#muB"," Mass mean B",5.279,5.250,5.3,"GeV");                                                                                                                       
 RooRealVar width2("width2"," Mass width2 ",0.025,0.020,0.07,"GeV");
 RooGaussian Sig2("Sig2"," Signal PDF B",M,mean,width2);

 RooRealVar d0("d0","d0",0.035,0.01,0.085);// another parameters related to "shoulder" shape on B+ mass BG
 RooRealVar d1("d1","d1",5.14,5.05,5.2);// describes the "shoulder" steepness   
 RooGenericPdf genpdf("genpdf","genpdf","(TMath::Erf((-M + d1)/d0)+1)",RooArgSet(d0,d1,M));

 //********final PDF ********
 RooRealVar Ns("Ns","Ns",0.,100000);
 RooRealVar Nb("Nb","Nb",0.,100000);   
 RooRealVar fs("fs","fs",0.5,0.,1.);
 RooRealVar fb("fb","fb",0.2,0.,1.);

 RooAddPdf sumgau("sumgau","sumgau",RooArgList(Sig,Sig2),RooArgList(fs));
 RooAddPdf sumbkg("sumbkg","add Erf and expB",RooArgList(genpdf,bkg1),RooArgList(fb)); 

 //RooAddPdf MassModel("MassModel","MassModel",RooArgList(Sig,bkg1),RooArgList(Ns,Nb));
 RooAddPdf MassModel("MassModel","MassModel",RooArgList(sumgau,bkg1),RooArgList(Ns,Nb));
 //RooAddPdf MassModel("MassModel","MassModel",RooArgList(sumgau,sumbkg),RooArgList(Ns,Nb));
//------------ Fit procedure -------------------
 Ns.setVal(1000.0);
 Nb.setVal(2000.0);
 //Ns.setConstant(kTRUE);

 RooFitResult* fitres = MassModel.fitTo(data,Extended(),Minos(kFALSE),Save(kTRUE), NumCPU(4));
 data.Print("v"); 
 fitres->Print("v");

 
 
}//End analysis
