#include "classreduceMC.C"

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

//#include "../../analisis/tdrstyle.C"

using namespace RooFit;
using namespace std;

void save_parameters(ofstream& salida, RooRealVar m,  RooRealVar g1, RooRealVar g2, RooRealVar fg)
{
  //salida <<  m.getVal() << " " << m.getError() << " " << g1.getVal() << " " << g1.getError()
  //	 << " " << g2.getVal() << " " << g2.getError() << " " << fg.getVal() << " " << fg.getError();
  salida <<  m.getVal() << " "  << g1.getVal() << " " << g2.getVal() << " " << fg.getVal();
  cout << " el archivo de parametros se escribio bien" << endl; 
  return;
}

void save_result(ofstream& salida, RooDataSet data,RooRealVar Ns, RooRealVar Nb)
{
  //salida <<  Ns.getVal() << " " << Ns.getError() << " " << Nb.getVal() << " " << Nb.getError();
  salida <<  (double)data.numEntries() << " " << sqrt(data.numEntries()) ;
  cout << " el archivo se escribio bien" << endl; 
  return;
}

TCanvas* CreateCanvasNomPull(TString cname, RooFitResult* result, RooDataSet data, RooRealVar M, Double_t supM, Double_t infM,  RooAbsPdf *MassModel, RooAbsPdf *sumgau, RooAbsPdf *bkg1, RooRealVar Ns, RooRealVar Nb, RooRealVar width, RooRealVar width2, RooRealVar fs, RooRealVar mean, Double_t ptl, Double_t pth, Int_t valsT)  
{
 //Double_t nbin = ((supM-infM)/0.010) + 1;
 Double_t nbin = ((supM-infM)/0.001);

 int H = 600;
 int W = 800;
 TCanvas *c1 = new TCanvas(cname,cname,50,50,W,H);
 //TCanvas *c1 = new TCanvas(cname,cname,W,H);
 //c1->Divide(1,2);
 //c1->cd(1) ;
 c1->cd() ;  
 c1->SetLeftMargin(0.005);
 c1->SetRightMargin(0.01);
 c1->SetTopMargin(0.09);
 c1->SetBottomMargin(0.1);
 //gPad->SetLogy();
 
 
 TPad *pad1 = new TPad("pad1", "padi",0.01,0.411,0.9903769, 0.99 );
 pad1->SetLeftMargin(0.09);   
 pad1->SetRightMargin(0.019);
 pad1->SetTopMargin(0.09);
 pad1->SetBottomMargin(0.0);
 
 
 TPad *pad2 = new TPad("pad2", "pad2",0.01,0.01,0.9903769,0.41);
 pad2->SetLeftMargin(0.09);
 pad2->SetRightMargin(0.019);  
 pad2->SetTopMargin(0.0);
 pad2->SetBottomMargin(0.25);
 //pad2->SetTickx(0);
 pad2->SetFillColor(0);
 pad2->SetGridx(0);
 pad2->SetGridy(0);
 
 pad1->Draw();
 pad2->Draw();
 pad1->cd();
 
 RooPlot* Mframe = M.frame(infM,supM,nbin);
 data.plotOn(Mframe,DataError(RooAbsData::SumW2),MarkerSize(1.0));
 //MassModel->plotOn(Mframe);
 MassModel->plotOn(Mframe,DrawOption("F"),FillColor(0),LineWidth(2),Name("fittotal"));

 RooHist* hpullm2 = Mframe->pullHist() ;
 
 //MassModel->plotOn(Mframe,Components(*sumgau),LineColor(kRed),LineWidth(2),Name("Signal")); 
 //MassModel->plotOn(Mframe,Components(*bkg1),LineColor(kGreen),LineWidth(2),Name("bkg"));
 //MassModel->plotOn(Mframe,Components(genpdf),LineColor(kOrange),LineWidth(2),Name("jpsikx"));  
 //MassModel->plotOn(Mframe,DrawOption("F"),FillColor(kBlue-9),LineColor(1),LineWidth(1),Name("fittotal"));
 //MassModel->plotOn(Mframe,Components(sumgau),DrawOption("F"),FillColor(kBlue-9),LineColor(kBlue-9),LineWidth(2),Name("Signal")); 
 //MassModel->plotOn(Mframe,Components(bkg1),DrawOption("F"),FillColor(kCyan),LineColor(kCyan),LineWidth(2),Name("bkg"));
 data.plotOn(Mframe,DataError(RooAbsData::SumW2),MarkerSize(1.0),XErrorSize(0),Name("Data"));
 MassModel->plotOn(Mframe);
 
 Mframe->SetYTitle("Events / 1 MeV");                                                                                                                   
 Mframe->SetLabelSize(0.07,"XY");
 Mframe->SetTitleSize(0.08,"XY");
 Mframe->GetYaxis()->CenterTitle();   
 Mframe->GetXaxis()->CenterTitle();
 Mframe->GetYaxis()->SetNdivisions(505,1);
 Mframe->GetXaxis()->SetNdivisions(505,1);   
 Mframe->GetXaxis()->SetDecimals(1); 
 Mframe->SetTitleOffset(0.8,"X");
 Mframe->SetTitleOffset(0.6,"Y");
 Mframe->SetMinimum(1.0); 
 Mframe->Draw();
 
 //TLegend *leg = new TLegend(0.15,0.68,0.35,0.88);
 TLegend *leg = new TLegend(0.18,0.6,0.38,0.88); 
 leg->SetTextSize(0.06);
 leg->SetFillColor(0);
 leg->SetBorderSize(0);
 leg->SetFillStyle(0);
 leg->AddEntry(Mframe->findObject("Data")," Data","ep"); 
 leg->AddEntry(Mframe->findObject("fittotal")," Fit result","l");
 //leg->AddEntry(Mframe->findObject("Signal")," Signal","l");
 //leg->AddEntry(Mframe->findObject("jpsikx"),"J/#psi kX","l"); 
 //leg->AddEntry(Mframe->findObject("bkg"),"Comb. backg.","l"); 
 //leg->AddEntry(Mframe->findObject("Signal")," Signal","f");
 //leg->AddEntry(Mframe->findObject("bkg"),"Comb. backg.","f");
 leg->Draw();
 
 Double_t Mpsi = mean.getVal()*1000.0;
 Double_t MpsiE = mean.getError()*1000.0;
 Double_t G ;
 Double_t GE;
 if(valsT==0){
    G = sqrt( fs.getVal()*width.getVal()*width.getVal() + (1-fs.getVal())*width2.getVal()*width2.getVal() )*1000.0;
    GE = (1/G)*sqrt( (fs.getVal()*fs.getVal())*(width.getVal()*width.getVal())*(width.getError()*width.getError()) + ((1-fs.getVal())*(1-fs.getVal()))*(width2.getVal()*width2.getVal())*(width2.getError()*width2.getError()) )*1000.0*1000.0;
  }
  else{
    G = width.getVal()*1000.0;
    GE = width.getError()*1000.0;
  }

 TLegend *legpar = new TLegend(0.6,0.6,0.8,0.88);
 legpar->SetTextSize(0.06);
 legpar->SetFillColor(0);
 legpar->SetBorderSize(0);
 legpar->SetFillStyle(0);
 legpar->AddEntry("",Form("M(B^{+}) = %1.2f #pm %1.2f MeV",Mpsi,MpsiE),"");
 legpar->AddEntry("",Form("#sigma = %1.2f #pm %1.2f MeV",G,GE),"");
 //legpar->AddEntry("",Form("N_{B^{+}} = %1.0f #pm %1.0f",Ns.getVal(),Ns.getError()),"");
 legpar->AddEntry("",Form("N_{B^{+}} = %1i #pm %1.0f",data.numEntries(),sqrt(data.numEntries())),"");
 //legpar->AddEntry(Form("%1.1f #leq p_{T}(B^{+}) < %1.1f GeV ",ptl,pth)); 
 legpar->Draw();

 TLegend *legMass = new TLegend(0.64,0.45,0.83,0.50);
 legMass->SetTextFont(43); 
 legMass->SetTextSize(20);  
 legMass->SetFillColor(0); 
 legMass->SetBorderSize(0);
 legMass->SetFillStyle(0); 
 legMass->SetHeader(Form("%1.1f #leq p_{T}(B^{+}) < %1.1f GeV ",ptl,pth));
 legMass->Draw();

 pad2->cd();

 // Create a new frame to draw the pull distribution 
 RooPlot* framem2 = M.frame(infM,supM,nbin) ;
 //framem2->addPlotable(hpullm2,"P") ;
 framem2->addPlotable(hpullm2,"XP") ;// setting Y error to cero
 
 framem2->SetYTitle(" (Data-Fit)/#sigma");
 framem2->SetLabelSize(0.1,"XY");
 framem2->SetTitleSize(0.13,"X");
 framem2->SetTitleSize(0.11,"Y");  
 framem2->GetYaxis()->CenterTitle();   
 framem2->GetXaxis()->CenterTitle();
 framem2->GetYaxis()->SetNdivisions(505,1);
 framem2->GetXaxis()->SetNdivisions(505,1);
 framem2->GetXaxis()->SetTickLength(0.07);  
 framem2->SetTitleOffset(0.9,"X");
 framem2->SetTitleOffset(0.4,"Y");
 framem2->SetMaximum(4.9);
 framem2->SetMinimum(-4.9);
 
 framem2->Draw();
 
 TLine *line1 = new TLine(infM,0.0,supM,0.0);
 line1->SetLineColor(1);
 line1->SetLineWidth(1);
 line1->Draw();
 TLine *line2 = new TLine(infM,4.0,supM,4.0);
 line2->SetLineColor(2);
 line2->SetLineWidth(2);
 line2->SetLineStyle(2);
 line2->Draw();
 TLine *line3 = new TLine(infM,-4.0,supM,-4.0);
 line3->SetLineColor(2);
 line3->SetLineWidth(2);
 line3->SetLineStyle(2);
 line3->Draw();
 
 c1->cd();

 //TLatex *   tex1 = new TLatex(0.92,0.926,"61.2 fb^{-1} (13 TeV, 2018)");
 TLatex *   tex1 = new TLatex(0.98,0.95,"MC simulation");

 tex1->SetNDC();
 tex1->SetTextAlign(31);
 tex1->SetTextFont(42);
 tex1->SetTextSize(0.05); 
 tex1->SetLineWidth(2);
 
 TLatex *tex2 = new TLatex(0.1,0.95,"CMS");
 tex2->SetNDC();
 tex2->SetTextFont(61);
 tex2->SetTextSize(0.05); 
 tex2->SetLineWidth(2);

 TLatex *tex3 = new TLatex(0.19,0.95,"Preliminary");
 tex3->SetNDC();
 tex3->SetTextFont(52);
 tex3->SetTextSize(0.05); 
 tex3->SetLineWidth(2);

 tex1->Draw();  
 tex2->Draw();
 tex3->Draw();

 c1->Modified();
 gPad->Update();
 gPad->RedrawAxis();
 TLine l;
 l.DrawLine(gPad->GetUxmax(), gPad->GetUymax(), gPad->GetUxmax(), gPad->GetUymin());
 return c1; 
 
}


TCanvas* CreateCanvasNom(TString cname, RooFitResult* result, RooDataSet data, RooRealVar M, Double_t supM, Double_t infM,  RooAbsPdf *MassModel, RooAbsPdf *sumgau, RooAbsPdf *bkg1, RooRealVar Ns, RooRealVar Nb, RooRealVar width, RooRealVar width2, RooRealVar fs, RooRealVar mean, Double_t ptl, Double_t pth, Int_t valsT)  
{

 //Double_t nbin = ((supM-infM)/0.010) + 1;
 Double_t nbin = ((supM-infM)/0.010);

 int H = 600;
 int W = 800;
 TCanvas *c1 = new TCanvas("c1","",50,50,W,H);
 c1->cd();
 c1->SetLeftMargin(0.09);
 c1->SetRightMargin(0.02);
 c1->SetTopMargin(0.09);
 c1->SetBottomMargin(0.13); 
 
 RooPlot* Mframe = M.frame(infM,supM,nbin);
 data.plotOn(Mframe,DataError(RooAbsData::SumW2),MarkerSize(1.0));
 //MassModel->plotOn(Mframe);
 MassModel->plotOn(Mframe,DrawOption("F"),FillColor(0),LineWidth(2),Name("fittotal"));

 Double_t nfloatpars = result->floatParsFinal().getSize();
 Double_t ndof = nbin - nfloatpars;
 Double_t chi2tmp = Mframe->chiSquare()*nbin;
 Double_t probChi2 = TMath::Prob(chi2tmp, ndof)*100.0;
 cout<<" Chi2/ndof : "<<Mframe->chiSquare()<< endl;
 cout<<" Ndof : "<<ndof<< endl;
 cout<<" Chi2 : "<<chi2tmp<< endl; 

 
 //MassModel->plotOn(Mframe,Components(*sumgau),LineColor(kRed),LineWidth(2),Name("Signal")); 
 //MassModel->plotOn(Mframe,Components(*bkg1),LineColor(kGreen),LineWidth(2),Name("bkg"));
 //MassModel->plotOn(Mframe,Components(genpdf),LineColor(kOrange),LineWidth(2),Name("jpsikx"));  
 //MassModel->plotOn(Mframe,DrawOption("F"),FillColor(kBlue-9),LineColor(1),LineWidth(1),Name("fittotal"));
 //MassModel->plotOn(Mframe,Components(sumgau),DrawOption("F"),FillColor(kBlue-9),LineColor(kBlue-9),LineWidth(2),Name("Signal")); 
 //MassModel->plotOn(Mframe,Components(bkg1),DrawOption("F"),FillColor(kCyan),LineColor(kCyan),LineWidth(2),Name("bkg"));
 data.plotOn(Mframe,DataError(RooAbsData::SumW2),MarkerSize(1.0),XErrorSize(0),Name("Data"));
 MassModel->plotOn(Mframe);
 
 Mframe->SetYTitle("Events / 10 MeV");                                                                                                                   
 Mframe->SetLabelSize(0.04,"XY");
 Mframe->SetTitleSize(0.05,"XY");
 Mframe->GetYaxis()->CenterTitle();   
 Mframe->GetXaxis()->CenterTitle();
 Mframe->GetYaxis()->SetNdivisions(505,1);
 Mframe->GetXaxis()->SetNdivisions(505,1);   
 Mframe->GetXaxis()->SetDecimals(1); 
 Mframe->SetTitleOffset(0.9,"X");
 Mframe->SetTitleOffset(0.7,"Y");
 Mframe->SetTitleSize(0.06,"XY");
 Mframe->SetMinimum(1.0); 
 Mframe->Draw();
 
 //TLegend *leg = new TLegend(0.15,0.68,0.35,0.88);
 TLegend *leg = new TLegend(0.18,0.68,0.38,0.88); 
 leg->SetTextSize(0.035);
 leg->SetFillColor(0);
 leg->SetBorderSize(0);
 leg->SetFillStyle(0);
 leg->AddEntry(Mframe->findObject("Data")," Data","ep"); 
 leg->AddEntry(Mframe->findObject("fittotal")," Fit result","l");
 //leg->AddEntry(Mframe->findObject("Signal")," Signal","l");
 //leg->AddEntry(Mframe->findObject("jpsikx"),"J/#psi kX","l"); 
 //leg->AddEntry(Mframe->findObject("bkg"),"Comb. backg.","l"); 
 //leg->AddEntry(Mframe->findObject("Signal")," Signal","f");
 //leg->AddEntry(Mframe->findObject("bkg"),"Comb. backg.","f");
 leg->Draw();
 
 Double_t Mpsi = mean.getVal()*1000.0;
 Double_t MpsiE = mean.getError()*1000.0;
 Double_t G ;
 Double_t GE;
 if(valsT==0){
    G = sqrt( fs.getVal()*width.getVal()*width.getVal() + (1-fs.getVal())*width2.getVal()*width2.getVal() )*1000.0;
    GE = (1/G)*sqrt( (fs.getVal()*fs.getVal())*(width.getVal()*width.getVal())*(width.getError()*width.getError()) + ((1-fs.getVal())*(1-fs.getVal()))*(width2.getVal()*width2.getVal())*(width2.getError()*width2.getError()) )*1000.0*1000.0;
  }
  else{
    G = width.getVal()*1000.0;
    GE = width.getError()*1000.0;
  }

 TLegend *legpar = new TLegend(0.6,0.68,0.8,0.88);
 legpar->SetTextSize(0.035);
 legpar->SetFillColor(0);
 legpar->SetBorderSize(0);
 legpar->SetFillStyle(0);
 legpar->AddEntry("",Form("M(B^{+}) = %1.2f #pm %1.2f MeV",Mpsi,MpsiE),"");
 legpar->AddEntry("",Form("#sigma = %1.2f #pm %1.2f MeV",G,GE),"");
 //legpar->AddEntry("",Form("N_{B^{+}} = %1.0f #pm %1.0f",Ns.getVal(),Ns.getError()),"");
 legpar->AddEntry("",Form("N_{B^{+}} = %1i #pm %1.0f",data.numEntries(),sqrt(data.numEntries())),"");
 //legpar->AddEntry(Form("%1.1f #leq p_{T}(B^{+}) < %1.1f GeV ",ptl,pth)); 
 legpar->Draw();

 TLegend *legMass = new TLegend(0.64,0.57,0.83,0.65);
 legMass->SetTextFont(43); 
 legMass->SetTextSize(20);  
 legMass->SetFillColor(0); 
 legMass->SetBorderSize(0);
 legMass->SetFillStyle(0); 
 legMass->SetHeader(Form("%1.1f #leq p_{T}(B^{+}) < %1.1f GeV ",ptl,pth));
 legMass->Draw(); 

 //TLatex *   tex1 = new TLatex(0.92,0.926,"61.2 fb^{-1} (13 TeV, 2018)");
 TLatex *   tex1 = new TLatex(0.92,0.926,"MC simulation");

 tex1->SetNDC();
 tex1->SetTextAlign(31);
 tex1->SetTextFont(42);
 tex1->SetTextSize(0.05); 
 tex1->SetLineWidth(2);
 
 TLatex *tex2 = new TLatex(0.2,0.926,"CMS");
 tex2->SetNDC();
 tex2->SetTextFont(61);
 tex2->SetTextSize(0.05); 
 tex2->SetLineWidth(2);

 TLatex *tex3 = new TLatex(0.29,0.926,"Preliminary");
 tex3->SetNDC();
 tex3->SetTextFont(52);
 tex3->SetTextSize(0.05); 
 tex3->SetLineWidth(2);

 tex1->Draw();  
 tex2->Draw();
 tex3->Draw();

 c1->Modified();
 gPad->Update();
 gPad->RedrawAxis();
 TLine l;
 l.DrawLine(gPad->GetUxmax(), gPad->GetUymax(), gPad->GetUxmax(), gPad->GetUymin());
 return c1; 
 
}

void BuFitMC_ptbins(int vq=1, Double_t ptl=7.0, Double_t pth=50.0)
{
gStyle->SetOptTitle(0);
gStyle->SetOptFit(0);
gStyle->SetOptStat(0);
//gStyle->SetErrorX(0);
 
Double_t Mmin = 5.0; 
Double_t Mmax = 5.6;
//Double_t Mmin = 5.0; 
//Double_t Mmax = 5.5; 

TString Lumi = ""; 
 
TChain *ch = new TChain("butree","");
ch->Add(Form("reducetree_Bujk_AOD_HI2016_MC_best%1i.root/butree",vq));

TTree *tree = (TTree*) ch;
classreduceMC t(tree);
Long64_t nentries = t.fChain->GetEntries();
cout<<" Entries : "<<nentries<<endl;
  
//------------------------------
RooRealVar M("M"," M(J/#psi K^{+}) (GeV)",Mmin,Mmax);
RooDataSet data("data","data",RooArgSet(M));
 
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

      //Muon cuts
      //if(t.mu1pt<2.0 || t.mu2pt<2.0) continue;
      if( abs(t.etamu1)<1.1 && t.mu1pt<3.3)continue;
      if( ( abs(t.etamu1)>1.1 && abs(t.etamu1)<2.1 ) && t.mu1pt<(5.5 - 2.0 * abs(t.etamu1)) )continue;
      if( abs(t.etamu1)>2.1 && t.mu1pt<1.5)continue;
      
      if( abs(t.etamu2)<1.1 && t.mu2pt<3.3)continue;
      if( ( abs(t.etamu2)>1.1 && abs(t.etamu2)<2.1 ) && t.mu2pt<(5.5 - 2.0 * abs(t.etamu2)) )continue;
      if( abs(t.etamu2)>2.1 && t.mu2pt<1.5)continue;
      
      //B fiducial region cuts
      //if(t.Bupt<7.0 || t.Bupt>=70.0)continue;
      if(t.Bupt<ptl || t.Bupt>=pth)continue; 
      if(abs(t.rapidityB)>2.4)continue;

      // optimized cuts
      if(t.pion1pt<1.0)continue;
      if(t.Bpro<0.05)continue;
      if((t.pdl/t.pdle)<4)continue;
            
      M=t.B_mass;    
      data.add(RooArgSet(M));
    }

 cout<<"xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"<<endl;
 data.Print("v");
 //return; 
 RooAbsReal::defaultIntegratorConfig()->setEpsAbs(1e-10) ;
 RooAbsReal::defaultIntegratorConfig()->setEpsRel(1e-10) ;
 
 //---- Mass model -------
 Double_t supM = Mmax;
 Double_t infM = Mmin;
 
 // ****define nominal background ****
 //RooRealVar a0("a0","a0",0.23,-10.0,10.0);
 //RooChebychev bkg1("bkg1","Background",M,RooArgList(a0)); 

 RooRealVar c("c","c",-2.0,-10.0,10.0);
 //RooExponential bkg1("bkg1","Exp. Background",M,c);
 RooAbsPdf *bkg1 = 0;
 bkg1 = new RooExponential("bkg1","Exp. Background",M,c);

 //**** define nominal signal ****
 //RooRealVar mean("mean"," Mass mean",3.686097);
 RooRealVar mean("mean"," Mass mean",5.279,5.250,5.3,"GeV");
 RooRealVar width("width"," Mass width",0.015,0.005,0.030,"GeV");
 RooGaussian Sig("Sig"," Signal PDF",M,mean,width);

 //RooRealVar meanB("#muB"," Mass mean B",5.279,5.250,5.3,"GeV"); 
 RooRealVar width2("width2"," Mass width2 ",0.040,0.020,0.07,"GeV");
 RooGaussian Sig2("Sig2"," Signal PDF B",M,mean,width2);

 //******** Final PDF ********
 RooRealVar Ns("Ns","Ns",0.,1000000);
 RooRealVar Nb("Nb","Nb",0.,100000);   
 RooRealVar fs("fs","fs",0.5,0.,1.);
 RooRealVar fs2("fs2","fs2",0.3,0.,1.); 
 RooRealVar fb("fb","fb",0.2,0.,1.);

 //RooAddPdf sumgau("sumgau","sumgau",RooArgList(Sig,Sig2),RooArgList(fs));
 RooAbsPdf *sumgau = 0;
 sumgau = new RooAddPdf("sumgau","sumgau",RooArgList(Sig,Sig2),RooArgList(fs));

 //RooAddPdf MassModel("MassModel","MassModel",RooArgList(Sig,bkg1),RooArgList(Ns,Nb));
 //RooAddPdf MassModel("MassModel","MassModel",RooArgList(sumgau,bkg1),RooArgList(Ns,Nb));

 RooAbsPdf *MassModel = 0;
 MassModel = new RooAddPdf("MassModel","MassModel",RooArgList(*sumgau,*bkg1),RooArgList(Ns,Nb));
 //MassModel = new RooAddPdf("MassModel","MassModel",RooArgList(Sig,Sig2),RooArgList(fs));
 //MassModel = new RooAddPdf("MassModel","MassModel",RooArgList(Sig,Sig2,Sig3),RooArgList(fs,fs2));

//------------ Fit procedure -------------------
 Ns.setVal(10000.0);
 Nb.setVal(200.0);
 //Ns.setConstant(kTRUE);

 RooFitResult* fitres = MassModel->fitTo(data,Extended(),Minos(kFALSE),Save(kTRUE), NumCPU(4));
 //RooFitResult* fitres = MassModel->fitTo(data,Minos(kFALSE),Save(kTRUE), NumCPU(4)); 
 data.Print("v"); 
 fitres->Print("v");

 //Double_t Gt = sqrt( fs.getVal()*width.getVal()*width.getVal() + (1-fs.getVal())*width2.getVal()*width2.getVal() )*1000.0;
 //Double_t GtE = (1/Gt)*sqrt( (fs.getVal()*fs.getVal())*(width.getVal()*width.getVal())*(width.getError()*width.getError()) + ((1-fs.getVal())*(1-fs.getVal()))*(width2.getVal()*width2.getVal())*(width2.getError()*width2.getError()) )*1000.0*1000.0;

 //parameters output file
 ofstream salida_parameters(Form("plots_ptbinsMC/parameters_BuFit_ptbins_%1.0f_%1.0f.txt",ptl,pth));
 salida_parameters.is_open();
 save_parameters(salida_parameters, mean, width, width2, fs);
 
 //yields output file
 ofstream salida_nominal(Form("plots_ptbinsMC/output_BuFit_%1i_ptbins_%1.0f_%1.0f.txt",vq,ptl,pth));
 salida_nominal.is_open();
 save_result(salida_nominal, data, Ns, Nb);
 
 //made canvas
 TCanvas* canv_nominal = CreateCanvasNom("canv_nominal", fitres, data, M, supM, infM, MassModel, sumgau, bkg1, Ns, Nb, width, width2, fs, mean, ptl, pth, 0); 
 canv_nominal->Print(Form("plots_ptbinsMC/mass_BuFitMC_best%1i_ptbins_%1.0f_%1.0f.png",vq,ptl,pth));
 canv_nominal->Print(Form("plots_ptbinsMC/mass_BuFitMC_best%1i_ptbins_%1.0f_%1.0f.pdf",vq,ptl,pth));
 //made canvas with pull
 //TCanvas* canv_nominalPull = CreateCanvasNomPull("canv_nominalPull", fitres, data, M, supM, infM, MassModel, sumgau, bkg1, Ns, Nb, width, width2, fs, mean, ptl, pth, 0); 
 //canv_nominalPull->Print(Form("plots_ptbinsMC/mass_BuFitMC_best%1i_ptbins_%1.0f_%1.0f_pull.png",vq,ptl,pth));
 //canv_nominalPull->Print(Form("plots_ptbinsMC/mass_BuFitMC_best%1i_ptbins_%1.0f_%1.0f_pull.pdf",vq,ptl,pth));
 
}//End analysis
