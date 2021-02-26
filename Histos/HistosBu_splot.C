#include <TROOT.h>
#include "TMath.h"
#include <iostream>
#include "RooFit.h"
#include "RooFitResult.h"
#include "RooGlobalFunc.h"
#include "RooArgusBG.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooExponential.h"
#include "RooGenericPdf.h"
#include "RooPolynomial.h"
#include "RooChebychev.h"
#include "TPaveText.h"
#include "RooAddPdf.h"
#include "RooWorkspace.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "TAxis.h"
#include "RooHist.h"
#include "RooProdPdf.h"
#include "RooFitResult.h"
#include "TStyle.h"
#include "Roo1DTable.h"
#include "RooExponential.h"
#include "RooAddPdf.h"
#include "RooDecay.h"
#include "RooProdPdf.h"
#include "TRandom3.h"
#include "RooDataHist.h"
#include "RooConstVar.h"
#include "RooGlobalFunc.h"
#include "RooGaussian.h"
#include "RooBinning.h"
#include "RooEffProd.h"
#include "RooGamma.h"
#include "RooWorkspace.h"
#include "RooGenericPdf.h"
#include "RooArgSet.h"
#include "RooEfficiency.h"
#include "RooAbsReal.h"
#include "RooMsgService.h"
#include "RooHist.h"
#include "RooFormulaVar.h"
#include "TH1F.h"
#include "RooChi2Var.h"
#include "RooEffProd.h"
#include "RooPlot.h"
#include <vector>
#include <cmath>
#include "TCut.h"
#include "TEfficiency.h"
#include  "RooRealVar.h"
#include "RooBinning.h"
#include "TFile.h"
#include "TLegend.h"
#include "RooChebychev.h"
#include "RooAbsPdf.h"
#include "RooCustomizer.h"
#include <iostream>
#include "RooDataHist.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TGaxis.h"
#include "TAxis.h"
#include "TLatex.h"
#include "TF1.h"

#include "RooFFTConvPdf.h"

using namespace RooFit;
using namespace std;

TCanvas* CreateCanvas(TString cname, TH1 *h1, TH1 *h2, TString yhname, TString xhname, TString lavel1, TString lavel2, Double_t Bin, Double_t Min, Double_t Max ) 
{

 h1->Scale(1.0/h1->Integral());
 h2->Scale(1.0/h2->Integral());
 //h2->Scale(h1->GetSumOfWeights()/h2->GetSumOfWeights());

 int binmax1 = h1->GetMaximumBin();
 Double_t Yscale1 = h1->GetBinContent(binmax1);
 //cout << "GetMaximu1:   "<< Yscale1 << endl;

 int binmax2 = h2->GetMaximumBin();
 Double_t Yscale2 = h2->GetBinContent(binmax2);
 //cout << "GetMaximu2:   "<< Yscale2 << endl;

 Double_t Ymax = 1.0;
 //Double_t Yplus = 0.2;
 //Double_t Yplus = 0.05;
 Double_t Yplus = 0.02;

 if(Yscale1>Yscale2){
   Ymax = Yscale1+Yplus;
 }
 else{
   Ymax = Yscale2+Yplus;
 }

 int H = 600;
 int W = 800;
 TCanvas* canv = new TCanvas(cname,cname,50,50,W,H);
 //canv->Divide(1,2);
 canv->cd();
 canv->Draw();
 canv->SetLeftMargin(0.01);
 //canv->SetRightMargin(0.01);
 canv->SetTopMargin(0.07);
 canv->SetBottomMargin(0.001);

 TPad *pad1 = new TPad("pad1", "padi",0.04962312,0.3151786,0.9503769, 0.92 );
 pad1->SetBottomMargin(0.015);
 pad1->SetRightMargin(0.03);
 pad1->SetLeftMargin(0.08); 
 pad1->SetTopMargin(0.02);

 TPad *pad2 = new TPad("pad2", "pad2",0.04962312,0.015,0.9503769,0.31);
 
 pad2->SetTopMargin(0.04);
 pad2->SetBottomMargin(0.3);
 pad2->SetRightMargin(0.03);
 pad2->SetLeftMargin(0.08); 
 pad2->SetTickx(0);
 pad2->SetFillColor(0);
 pad2->SetGridx(0);
 pad2->SetGridy(0);

 pad1->Draw();
 pad2->Draw();
 pad1->cd();

 //h1->Draw("COLZ");
 h1->Draw("");
 h1->SetMarkerStyle(25);
 //h1->SetMarkerSize(1.5);
 h1->SetMarkerColor(2);
 h1->SetLineColor(2);
 h1->SetLineWidth(1);
 h1->GetYaxis()->CenterTitle(true);
 //h1->SetYTitle("Pt(J/#psi) [GeV]");
 h1->SetYTitle(yhname);
 //h1->SetYTitle("Distribution normalized."+yhname);
 h1->GetXaxis()->CenterTitle(true); 
 h1->SetXTitle(""); 
 //h1->SetXTitle(xhname); //#eta(#pi_{Bc} 
 h1->SetTitleSize(20,"Y"); 
 h1->SetLabelSize(20,"Y");
 h1->SetTitleOffset(1.2,"Y");
 h1->SetLabelFont(33,"Y");  
 h1->SetTitleFont(43,"Y");
 h1->SetLabelFont(3,"X");  
 h1->SetTitleFont(43,"X");
 //h1->SetMinimum(0.0);
 h1->SetMinimum(1e-5);
 h1->SetMaximum(Ymax); //para los piones de bajo pt
 
 //h2->Draw("same hits");
 h2->Draw("same"); 
 h2->SetMarkerStyle(20);
 h2->SetMarkerColor(4);
 h2->SetMarkerSize(1.0); 
 
 h2->SetLineColor(4);
 h2->SetLineWidth(1);
 //h2->SetMinimum(1.0);

 TLegend* histo_massleg = new TLegend(0.7,0.8,0.9,0.9); // aca ubicamos el tamaño de la legend
 histo_massleg->SetFillColor(0);
 histo_massleg->SetBorderSize(0);
 //histo_massleg->AddEntry(histo_mass,"Normal ", "lep");// si pone este ("lep") sale una cruz en ves de una linea
 histo_massleg->AddEntry(h1,lavel1, "ep");
 histo_massleg->AddEntry(h2,lavel2, "ep");
 histo_massleg->Draw("same");

 pad2->cd();
 
 TH1* h3 = new TH1F("h3","h3",Bin,Min,Max);
 h3->Sumw2();    
 h3->Divide(h1,h2,1,1);
 h3->Draw("");
 h3->SetMarkerStyle(2);
 h3->GetYaxis()->CenterTitle(true);
 h3->SetYTitle("Data/MC");
 h3->GetXaxis()->CenterTitle(true); 
 h3->SetXTitle(xhname);
 //h3->GetYaxis()->SetNdivisions(11006);10404
 h3->GetYaxis()->SetNdivisions(10004);// 
 h3->SetTitleSize(20,"XY"); 
 h3->SetLabelSize(20,"XY");
 h3->SetTitleOffset(1.2,"Y");
 h3->SetTitleOffset(3.6,"X");
 h3->SetLabelFont(33,"XY");  
 h3->SetTitleFont(43,"XY");
 h3->SetMinimum(0.0);
 h3->SetMaximum(2.0); 
 
 TLine *line = new TLine(Min,1,Max,1);
 line->SetLineWidth(2);
 line->SetLineColor(3);
 line->SetLineStyle(6);
 line->Draw("same");

 canv->cd();
 
 TLatex *   tex1 = new TLatex(0.88,0.926,"Data VS MC");
 tex1->SetNDC();
 tex1->SetTextAlign(31);
 tex1->SetTextFont(42);
 tex1->SetTextSize(0.04); 
 tex1->SetLineWidth(2);
 
 TLatex *tex2 = new TLatex(0.2,0.926,"CMS");
 tex2->SetNDC();
 tex2->SetTextFont(61);
 tex2->SetTextSize(0.04); 
 tex2->SetLineWidth(2);

 TLatex *tex3 = new TLatex(0.28,0.926,"Preliminary");
 tex3->SetNDC();
 tex3->SetTextFont(52);
 tex3->SetTextSize(0.04); 
 tex3->SetLineWidth(2); 

 //tex1->Draw();
 tex2->Draw();
 tex3->Draw();

 return canv;
}


void HistosBu_splot(){

 gStyle->SetOptTitle(0);
 gStyle->SetOptFit(0);
 gStyle->SetOptStat(0);
 gStyle->SetErrorX(0);
 //setTDRStyle();

 TString BMass = "M";
 TString yAxisTitlePDL = "Events / 10 MeV";
 TString xAxisTitleMass = "M(J/#psi K^{+}) (GeV)";
 TString signalFitTitle = "M(J/#psi K^{+}) (GeV)"; 
 
 // Open the rooworkspace that contains data, pdfs and observables 
 TFile *g = TFile::Open("MyBuws.root"); 
 RooWorkspace* wbc = (RooWorkspace*)g->Get("MyBu");

 TFile *g2MC = TFile::Open("MyBuMCws.root"); 
 RooWorkspace* wbc2MC = (RooWorkspace*)g2MC->Get("MyBuMC");
 
 RooRealVar* M = (RooRealVar*) wbc->var( BMass );
 M->SetTitle(xAxisTitleMass);
 M->setUnit("");
 //Double_t nbin = ( ( M->getMax() - M->getMin() )/0.010 );

 RooDataSet* Data = (RooDataSet*) wbc->data("dataWithSWeightsBu");
 RooDataSet* Data2MC = (RooDataSet*) wbc2MC->data("dataWithSWeightsBuMC");

 Double_t etaBmin = 0.0;
 //Double_t etaBmin = -2.4;  
 Double_t etaBmax = 2.4; 
 Double_t etaBbin = ( ( etaBmax- etaBmin )/0.2 );

 TH1* h1_etamul = Data->createHistogram("etamul",etaBbin,etaBmin,etaBmax) ;
 TH1* h1MC_etamul = Data2MC->createHistogram("etamul",etaBbin,etaBmin,etaBmax) ;

 TCanvas* canv_histoetamul = CreateCanvas("canv_histoetamul",  h1_etamul, h1MC_etamul, "Normalized distribution per event", "|y|(#mu_{1})", "Data B^{+}", "MC B^{+}", etaBbin, etaBmin, etaBmax); 
 canv_histoetamul->SaveAs("plots/Histo_DatavsMC_etamul.png");
 canv_histoetamul->SaveAs("plots/Histo_DatavsMC_etamul.pdf");

 TH1* h1_etamusl = Data->createHistogram("etamusl",etaBbin,etaBmin,etaBmax) ;
 TH1* h1MC_etamusl = Data2MC->createHistogram("etamusl",etaBbin,etaBmin,etaBmax) ;

 TCanvas* canv_histoetamusl = CreateCanvas("canv_histoetamusl",  h1_etamusl, h1MC_etamusl, "Normalized distribution per event", "|y|(#mu_{2})", "Data B^{+}", "MC B^{+}", etaBbin, etaBmin, etaBmax); 
 canv_histoetamusl->SaveAs("plots/Histo_DatavsMC_etamusl.png");
 canv_histoetamusl->SaveAs("plots/Histo_DatavsMC_etamusl.pdf");

 TH1* h1_etatkl = Data->createHistogram("etatkl",etaBbin,etaBmin,etaBmax) ;
 TH1* h1MC_etatkl = Data2MC->createHistogram("etatkl",etaBbin,etaBmin,etaBmax) ;

 TCanvas* canv_histoetatkl = CreateCanvas("canv_histoetatkl",  h1_etatkl, h1MC_etatkl, "Normalized distribution per event", "|y|(Tk_{1})", "Data B^{+}", "MC B^{+}", etaBbin, etaBmin, etaBmax); 
 canv_histoetatkl->SaveAs("plots/Histo_DatavsMC_etatkl.png");
 canv_histoetatkl->SaveAs("plots/Histo_DatavsMC_etatkl.pdf");

 //return;
 Double_t pttkmin = 0.8; 
 Double_t pttkmax = 4.8; 
 Double_t pttkbin = ( ( pttkmax- pttkmin )/0.5 );
 TH1* h1_pttkl = Data->createHistogram("pttkl",pttkbin,pttkmin,pttkmax);
 TH1* h1MC_pttkl = Data2MC->createHistogram("pttkl",pttkbin,pttkmin,pttkmax) ;

 TCanvas* canv_histopttkl = CreateCanvas("canv_histopttkl", h1_pttkl, h1MC_pttkl, "Normalized distribution per event", "p_{T}(Tk_{1}) (GeV)", "Data B^{+}", "MC B^{+}", pttkbin,pttkmin,pttkmax);
 canv_histopttkl->SaveAs("plots/Histo_DatavsMC_pttkl.png");
 canv_histopttkl->SaveAs("plots/Histo_DatavsMC_pttkl.pdf");
 
 //return;
 Double_t ptmuslmin = 1.5; 
 Double_t ptmuslmax = 5.0; 
 Double_t ptmuslbin = ( ( ptmuslmax - ptmuslmin )/0.5 );
 TH1* h1_ptmusl = Data->createHistogram("ptmusl",ptmuslbin,ptmuslmin,ptmuslmax) ;
 TH1* h1MC_ptmusl = Data2MC->createHistogram("ptmusl",ptmuslbin,ptmuslmin,ptmuslmax) ;
 
 TCanvas* canv_histoptmusl = CreateCanvas("canv_histoptmusl", h1_ptmusl, h1MC_ptmusl, "Normalized distribution per event", "p_{T}(#mu_{2}) (GeV)", "Data B^{+}", "MC B^{+}", ptmuslbin,ptmuslmin,ptmuslmax);
 canv_histoptmusl->SaveAs("plots/Histo_DatavsMC_ptmusl.png");
 canv_histoptmusl->SaveAs("plots/Histo_DatavsMC_ptmusl.pdf");

 //return;
 Double_t ptmumin = 1.5; 
 Double_t ptmumax = 15.5; 
 Double_t ptmubin = ( ( ptmumax- ptmumin )/0.5 );
 TH1* h1_ptmul = Data->createHistogram("ptmul",ptmubin,ptmumin,ptmumax);
 TH1* h1MC_ptmul = Data2MC->createHistogram("ptmul",ptmubin,ptmumin,ptmumax) ;

 TCanvas* canv_histoptmul = CreateCanvas("canv_histoptmul", h1_ptmul, h1MC_ptmul, "Normalized distribution per event", "p_{T}(#mu_{1}) (GeV)", "Data B^{+}", "MC B^{+}", ptmubin,ptmumin,ptmumax);
 canv_histoptmul->SaveAs("plots/Histo_DatavsMC_ptmul.png");
 canv_histoptmul->SaveAs("plots/Histo_DatavsMC_ptmul.pdf");

 
 //return;
 Double_t ptBmin = 7.0; 
 Double_t ptBmax = 49.0; 
 Double_t ptBbin = ( ( ptBmax- ptBmin )/2.0 );
 TH1* h1_ptB = Data->createHistogram("ptBu",ptBbin,ptBmin,ptBmax) ;
 TH1* h1MC_ptB = Data2MC->createHistogram("ptBu",ptBbin,ptBmin,ptBmax) ;

 TCanvas* canv_histoptB = CreateCanvas("canv_histoptB", h1_ptB, h1MC_ptB, "Normalized distribution per event", "p_{T}(J/#psi K^{+}) (GeV)", "Data B^{+}", "MC B^{+}", ptBbin, ptBmin, ptBmax);
 canv_histoptB->SaveAs("plots/Histo_DatavsMC_ptB.png");
 canv_histoptB->SaveAs("plots/Histo_DatavsMC_ptB.pdf");

 //return;
 TH1* h1_etaB = Data->createHistogram("etaBu",etaBbin,etaBmin,etaBmax) ;
 TH1* h1MC_etaB = Data2MC->createHistogram("etaBu",etaBbin,etaBmin,etaBmax) ;

 TCanvas* canv_histoetaB = CreateCanvas("canv_histoetaB",  h1_etaB, h1MC_etaB, "Normalized distribution per event", "|y|(J/#psi K^{+})", "Data B^{+}", "MC B^{+}", etaBbin, etaBmin, etaBmax); 
 canv_histoetaB->SaveAs("plots/Histo_DatavsMC_etaB.png");
 canv_histoetaB->SaveAs("plots/Histo_DatavsMC_etaB.pdf");

 
 
}




