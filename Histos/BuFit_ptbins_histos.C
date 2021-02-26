#include "../classreduce.C"

#include <TROOT.h>
#include "TMath.h"
#include <iostream>
#include <fstream>

#include "RooStats/SPlot.h"
#include "RooWorkspace.h"

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

void DoFit(RooDataSet data,  RooRealVar M, Double_t maxM, Double_t minM, TString Lumi, Double_t ptl, Double_t pth);

TCanvas* CreateCanvasNom(TString cname, RooFitResult* result, RooDataSet data, RooRealVar M, Double_t supM, Double_t infM,  RooAbsPdf* MassModel, RooAbsPdf* sumgau, RooAbsPdf *bkg1, RooAbsPdf* genpdf, RooRealVar Ns, RooRealVar Nb, RooRealVar width, RooRealVar width2, RooRealVar fs, RooRealVar mean, Double_t ptl, Double_t pth, TString Lum, Int_t valsT)
{

  Double_t nbin = ((supM-infM)/0.010) + 1;
  //Double_t nbin = ((supM-infM)/0.010);
  //Double_t nbin = ((supM-infM)/0.020);

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
  data.plotOn(Mframe,DataError(RooAbsData::SumW2),MarkerSize(1.0),XErrorSize(0));
  //MassModel->plotOn(Mframe);
  MassModel->plotOn(Mframe,DrawOption("F"),FillColor(0),LineWidth(2),Name("fittotal"));
  
  RooHist* hpullm2 = Mframe->pullHist() ;
  
  MassModel->plotOn(Mframe,Components(*sumgau),LineColor(kRed),LineWidth(2),Name("Signal"));
  MassModel->plotOn(Mframe,Components(*genpdf),DrawOption("F"),FillColor(kOrange),LineColor(kOrange),LineWidth(2),Name("jpsikx"));
  MassModel->plotOn(Mframe,Components(*bkg1),LineStyle(2),LineColor(kGreen),LineWidth(2),Name("bkg"));  
  //MassModel->plotOn(Mframe,Components(Gaussjk),LineStyle(6),LineColor(kViolet),LineWidth(2),Name("jpsipi")); 
  data.plotOn(Mframe,DataError(RooAbsData::SumW2),MarkerSize(1.0),XErrorSize(0),Name("Data"));
  MassModel->plotOn(Mframe);
  
  Mframe->SetYTitle("Events / 10 MeV");                                                                                                                   
  Mframe->SetLabelSize(0.07,"XY");
  Mframe->SetTitleSize(0.08,"XY");
  Mframe->GetYaxis()->CenterTitle();   
  Mframe->GetXaxis()->CenterTitle();
  Mframe->GetYaxis()->SetNdivisions(505,1);
  Mframe->GetXaxis()->SetNdivisions(505,1);   
  Mframe->GetXaxis()->SetDecimals(1); 
  Mframe->SetTitleOffset(0.9,"X");
  Mframe->SetTitleOffset(0.6,"Y");
  Mframe->SetMinimum(1.0);
  //Mframe->SetMaximum(800.0);   
  Mframe->Draw();
  
  //TLegend *leg = new TLegend(0.15,0.68,0.35,0.88);
  TLegend *leg = new TLegend(0.18,0.48,0.38,0.88); 
  leg->SetTextSize(0.06);
  leg->SetTextFont(42);  
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(Mframe->findObject("Data")," Data","ep"); 
  leg->AddEntry(Mframe->findObject("fittotal")," Fit result","l");
  leg->AddEntry(Mframe->findObject("Signal")," Signal","l");
  leg->AddEntry(Mframe->findObject("bkg"),"Comb. backg.","l");
  leg->AddEntry(Mframe->findObject("jpsikx"),"J/#psi K^{+} + X","f");
  //leg->AddEntry(Mframe->findObject("jpsipi"),"J/#psi #pi^{+}","l");    
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
  
  TLegend *legpar = new TLegend(0.6,0.58,0.8,0.88);
  legpar->SetTextSize(0.06);
  legpar->SetTextFont(42);
  legpar->SetFillColor(0);
  legpar->SetBorderSize(0);
  legpar->SetFillStyle(0);
  legpar->AddEntry("",Form("M(B^{+}) = %1.2f #pm %1.2f MeV",Mpsi,MpsiE),"");
  //legpar->AddEntry("",Form("#sigma = %1.2f #pm %1.2f MeV",G,GE),"");
  legpar->AddEntry("",Form("N_{B^{+}} = %1.0f #pm %1.0f",Ns.getVal(),Ns.getError()),"");
  legpar->AddEntry("",Form("N_{bkg} = %1.0f #pm %1.0f",Nb.getVal(),Nb.getError()),"");
  legpar->Draw();

  TLegend *legMass = new TLegend(0.64,0.35,0.83,0.4);
  legMass->SetTextFont(42); 
  legMass->SetTextSize(0.06);  
  legMass->SetFillColor(0); 
  legMass->SetBorderSize(0);
  legMass->SetFillStyle(0); 
  legMass->SetHeader(Form("%1.1f #leq p_{T}(B^{+}) < %1.1f GeV ",ptl,pth));
  legMass->Draw();
  
  pad2->cd();
  
  // Create a new frame to draw the pull distribution 
  RooPlot* framem2 = M.frame(infM,supM,nbin) ;
  framem2->addPlotable(hpullm2,"P") ;

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
  framem2->SetMaximum(3.9);
  framem2->SetMinimum(-3.9);
  
  framem2->Draw();

  TLine *line1 = new TLine(infM,0.0,supM,0.0);
  line1->SetLineColor(1);
  //line1->SetLineStyle(2);
  line1->SetLineWidth(1);
  line1->Draw();
  
  c1->cd();

 TLatex *   tex1 = new TLatex(0.98,0.95,Lum); 
 tex1->SetNDC();
 tex1->SetTextAlign(31);
 tex1->SetTextFont(42);
 tex1->SetTextSize(0.05); 
 tex1->SetLineWidth(2);
 
 TLatex *tex2 = new TLatex(0.10,0.95,"CMS");
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


void BuFit_ptbins_histos(Double_t ptl=7.0, Double_t pth=50.0)
{
 gStyle->SetOptTitle(0);
 gStyle->SetOptFit(0);
 gStyle->SetOptStat(0);
 //gStyle->SetErrorX(0);
 
 Double_t Mmin = 5.0; 
 Double_t Mmax = 5.6;
 TString Lumi = "";  
 Lumi = "189.7 nb^{-1} (#it{p}Pb 8.16 TeV)";

 TChain *ch = new TChain("butree","");
 ch->Add("../reducetree_Bujk_AOD_HI2016_Era1_best1.root/butree");
 ch->Add("../reducetree_Bujk_AOD_HI2016_Era2_best1.root/butree");

 TTree *tree = (TTree*) ch;
 classreduce t(tree);
 Long64_t nentries = t.fChain->GetEntries();
 cout<<" Entries : "<<nentries<<endl;
 
 //------------------------------
 RooRealVar M("M"," M(J/#psi K^{+}) (GeV)",Mmin,Mmax);
 //RooDataSet data("data","data",RooArgSet(M));

 RooRealVar etaBu("etaBu","etaBu",0.0,2.4);  
 RooRealVar ptBu("ptBu","ptBu",7.0,50.0);
 RooRealVar ptmul("ptmul","ptmul",1.5,20.5);
 RooRealVar ptmusl("ptmusl","ptmusl",1.5,20.5);
 RooRealVar etamul("etamul","etamul",0.0,2.4);
 RooRealVar etamusl("etamusl","etamusl",0.0,2.4);
 RooRealVar pttkl("pttkl","pttkl",0.8,10.8);
 RooRealVar etatkl("etatkl","etatkl",0.0,2.4);
 
 RooDataSet data("data","data",RooArgSet(M,etaBu,ptBu,ptmul,ptmusl,pttkl,etamul,etamusl,etatkl)); 
 
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

     //Muon cuts
     //if(t.mu1pt<2.0 || t.mu2pt<2.0) continue;
     if( abs(t.etamu1)<1.1 && t.mu1pt<3.3)continue;
     if( ( abs(t.etamu1)>1.1 && abs(t.etamu1)<2.1 ) && t.mu1pt<(5.5 - 2.0 * abs(t.etamu1)) )continue;
     if( abs(t.etamu1)>2.1 && t.mu1pt<1.5)continue;

     if( abs(t.etamu2)<1.1 && t.mu2pt<3.3)continue;
     if( ( abs(t.etamu2)>1.1 && abs(t.etamu2)<2.1 ) && t.mu2pt<(5.5 - 2.0 * abs(t.etamu2)) )continue;
     if( abs(t.etamu2)>2.1 && t.mu2pt<1.5)continue;

     //B fiducial region cuts
     if(t.Bupt<7.0 || t.Bupt>=50.0)continue;
     //if(t.Bupt<ptl || t.Bupt>=pth)continue; 
     if(abs(t.rapidityB)>2.4)continue;

     // optimized cuts
     if(t.pion1pt<0.8)continue;
     if(t.Bpro<0.01)continue;
     if((t.pdl/t.pdle)<6)continue;
  
     M=t.B_mass;    
     //data.add(RooArgSet(M));
     etaBu = fabs(t.rapidityB);
     ptBu = t.Bupt;
     ptmul = t.mu1pt;
     ptmusl = t.mu2pt;
     etamul = fabs(t.etamu1);
     etamusl = fabs(t.etamu2);
     pttkl = t.pion1pt;
     etatkl = fabs(t.etapi1);
     data.add(RooArgSet(M,etaBu,ptBu,ptmul,ptmusl,pttkl,etamul,etamusl,etatkl));

   }

cout<<"xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"<<endl;
data.Print("v");
//return;

DoFit(data, M, Mmax, Mmin, Lumi, ptl, pth);


}//End "main" function 

void DoFit(RooDataSet data,  RooRealVar M, Double_t maxM, Double_t minM, TString Lumi, Double_t ptl, Double_t pth)
{
 RooAbsReal::defaultIntegratorConfig()->setEpsAbs(1e-10) ;
 RooAbsReal::defaultIntegratorConfig()->setEpsRel(1e-10) ;

 //---- Mass model -------
 Double_t supM = maxM;
 Double_t infM = minM;
 //Int_t valsT = 0;
 
 //***********************************
 //   define combinatorial background
 //**********************************
 RooRealVar c("c","c",-2.0,-10.0,10.0);
 //RooExponential bkg1("bkg1","Exp. Background",M,c);
 RooAbsPdf *bkg1 = 0;
 bkg1 = new RooExponential("bkg1","Exp. Background",M,c);

 //**************************
 //   Define nominal signal
 //**************************
 // **** gaussians ****
 RooRealVar mean("mean"," Mass mean",5.279,5.250,5.3,"GeV");
 RooRealVar width("width"," Mass width",0.015,0.005,0.030,"GeV");
 RooGaussian Sig("Sig"," Signal PDF",M,mean,width);

 RooRealVar width2("width2"," Mass width2 ",0.040,0.020,0.07,"GeV");
 RooGaussian Sig2("Sig2"," Signal PDF B",M,mean,width2);

 // ****define J/psi K + X background ****
 RooAbsPdf *genpdf = 0;
 RooRealVar d0("d0","d0",0.035,0.001,0.085);// another parameters related to "shoulder" shape on B+ mass BG
 RooRealVar d1("d1","d1",5.14,5.05,5.2);// describes the "shoulder" steepness   
 //RooGenericPdf genpdf("genpdf","genpdf","(TMath::Erf((-M + d1)/d0)+1)",RooArgSet(d0,d1,M));
 //RooGenericPdf genpdf("genpdf","genpdf","(TMath::Erf((-@0 + @1)/@2)+1)",RooArgSet(M,d1,d0));
 //genpdf = new RooGenericPdf ("genpdf","genpdf","(TMath::Erf((-@0 + @1)/@2)+1)",RooArgSet(M,d1,d0));

 RooRealVar mbd("mbd","mbd",5.14,5.0,5.2);
 RooRealVar wbd("wbd","wbd",0.04,0.020,0.1,"GeV");
 genpdf = new RooGaussian("genpdf","genpdf",M,mbd,wbd);

 // ******** final PDF ********
 RooRealVar Ns("Ns","Ns",0.,10000);
 RooRealVar Nb("Nb","Nb",0.,10000);   
 RooRealVar fs("fs","fs",0.7,0.,1.);
 RooRealVar fb("fb","fb",0.2,0.,1.);

 RooAbsPdf *sumgau = 0;
 sumgau = new RooAddPdf("sumgau","sumgau",RooArgList(Sig,Sig2),RooArgList(fs));
 RooAbsPdf *sumbkg =0;
 sumbkg = new RooAddPdf("sumbkg","add Erf and bkg1",RooArgList(*genpdf,*bkg1),RooArgList(fb)); 

 RooAbsPdf *MassModel = 0;
 //RooAddPdf MassModel("MassModel","MassModel",RooArgList(Sig,bkg1),RooArgList(Ns,Nb));
 //RooAddPdf MassModel("MassModel","MassModel",RooArgList(sumgau,bkg1),RooArgList(Ns,Nb));
 //RooAddPdf MassModel("MassModel","MassModel",RooArgList(sumgau,*sumbkg),RooArgList(Ns,Nb));
 MassModel = new RooAddPdf("MassModel","MassModel",RooArgList(*sumgau,*sumbkg),RooArgList(Ns,Nb));

//------------ Fit procedure -------------------
 Ns.setVal(1500.0);
 Nb.setVal(1000.0);

 RooFitResult* fitres = MassModel->fitTo(data,Extended(),Minos(kFALSE),Save(kTRUE), NumCPU(4));
 data.Print("v"); 
 fitres->Print("v");

 //made canvas
 //TCanvas* canv_nominal = CreateCanvasNom("canv_nominal", fitres, data, M, supM, infM, MassModel, sumgau, bkg1, genpdf, Ns, Nb, width, width2, fs, mean, ptl, pth, Lumi, 0); 
 //canv_nominal->Print(Form("plots_ptbins/mass_BuFit_ptbins_%1.0f_%1.0f.png",ptl,pth));
 //return;

 //https://root.cern/doc/v610/rs301__splot_8C_source.html
 //https://root.cern.ch/doc/master/rs301__splot_8C.html
 //https://root.cern.ch/root/html/tutorials/roostats/rs301_splot.C.html
 // The sPlot technique requires that we fix the parameters
 // of the model that are not yields after doing the fit.

 c.setConstant(kTRUE);  
 mbd.setConstant(kTRUE); 
 wbd.setConstant(kTRUE);
 fs.setConstant(kTRUE);
 fb.setConstant(kTRUE);
 mean.setConstant(kTRUE);
 width.setConstant(kTRUE);
 width2.setConstant(kTRUE);

 RooMsgService::instance().setSilentMode(true);

 // Now we use the SPlot class to add SWeights to our data set
 // based on our model and our yield variables
 RooStats::SPlot* sData = new RooStats::SPlot("sData","sData", data, MassModel, RooArgList(Ns,Nb) );

 std::cout << "\n\nThe dataset after creating sWeights:\n";
 data.Print();
 //return;

 // Check that our weights have the desired properties
 std::cout << "Check SWeights:" << std::endl;

 std::cout <<std::endl<<"Yield of peak1 is "<< Ns.getVal() << ".  From sWeights it is " << sData->GetYieldFromSWeight("Ns") << std::endl;
 std::cout <<std::endl<<"Yield of bkg is "<< Nb.getVal() << ".  From sWeights it is " << sData->GetYieldFromSWeight("Nb") << std::endl;
 
 RooWorkspace* wspace = new RooWorkspace("myWSBu");
 wspace->import(data, Rename("dataWithSWeightsBu"));
 //RooDataSet dataWS2(dataWS.GetName(),dataWS.GetTitle(),&dataWS,*dataWS.get(),0,w->GetName()) ;
 RooDataSet* datasp = (RooDataSet*) wspace->data("dataWithSWeightsBu");
 cout << "datasp name: " << datasp->GetName() << endl;

 // create weightfed data set for peak signal
 RooDataSet * dataw_s1 = new RooDataSet(datasp->GetName(),datasp->GetTitle(),datasp,*datasp->get(),0,"Ns_sw") ;

 RooWorkspace myWSBu("MyBu");
 // Import data into the workspace
 myWSBu.import(*dataw_s1) ;
 cout << "Workspace contains " << endl;
 myWSBu.Print();
 // save your workspace in a file
 myWSBu.writeToFile("MyBuws.root");
 
}


