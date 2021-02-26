#define newTree_cxx
#include "newTree.h"

#include <TStyle.h>
#include "TMath.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include "TPaveStats.h"
#include "TCut.h"
#include "TH1.h"
#include "TF1.h"
#include "TLegend.h"



#include "RooFit.h"
#include "RooRealVar.h"
#include "RooCategory.h"
#include "RooAbsPdf.h"
#include "RooResolutionModel.h"
#include "RooFitResult.h"
#include "RooDataSet.h"
#include "RooGlobalFunc.h"
#include "RooProdPdf.h"
#include "RooGaussian.h"
#include "RooAddPdf.h"
#include "RooPlot.h"
#include "RooGaussModel.h"
#include "RooTruthModel.h"
#include "RooLandau.h"
#include "RooProdPdf.h"
#include "RooChebychev.h"
#include "RooExponential.h"
#include "RooDecay.h"
#include "RooPolynomial.h"
#include "RooGenericPdf.h"
#include "RooHistPdf.h"
#include "RooHist.h"
#include "RooBreitWigner.h"
#include "RooVoigtian.h"
#include "RooWorkspace.h"

#include "RooDataHist.h"
#include "RooKeysPdf.h"
#include "RooNDKeysPdf.h"

#include <vector>

using namespace RooFit;

Bool_t dominos = kTRUE;
Int_t ncpu = 6;
Int_t binfit=100;

TH1F* bkgsubs(vector<TTree*> tree,vector<string> hname, vector<TCut> cut, Double_t frac, string fillfunc, string descr, string savename, Double_t inibin, Double_t finbin, Int_t bin=100, Int_t normalize=0,Bool_t doweight=0,string xdescr="",Bool_t log=0, string addstr="bkgsubs");

void  dobkgsub(vector<TTree*> tree,vector<string> hname, vector<TCut> cut, Double_t frac,string addstr="bkgsubs");
void  weight(vector<TTree*> tree,vector<string> hname, vector<TCut> cut, Double_t frac,string addstr="weight");

Double_t turn(Double_t *x, Double_t* p) { 
  Double_t xx = x[0];
  return (2.0*p[0]-p[3])*TMath::Erf(p[1]*(xx-p[2]))/2.0 + p[3]/2.0; 
}

string thedataset = "";
Double_t min_Masa = 0;
Double_t max_Masa = 0;
Double_t min_J_lifetime_sig_data = 0;
Double_t max_J_lifetime_sig_data = 0;
Double_t min_J_lifetime_sig_mc = 0;
Double_t max_J_lifetime_sig_mc = 0;
Double_t min_abs_J_eta1 = 0;
Double_t max_abs_J_eta1 = 0;
Double_t min_muon_thresh = 0;
string thecutset = "";

void  reweight(string datasample="run2a", string cutset = "psi2s_central", string fitset="psi2s", int ngauss=1 ,int bkgdeg=1, bool compareonly = false) {
  TChain *tree = new TChain("Tree");
  thedataset = datasample;

  ifstream ifstextfile( datasample.c_str() );
  string templine;
  while( getline( ifstextfile, templine ) ) tree->Add(templine.c_str());

  newTree *Tree = new newTree(tree); 

  RooRealVar* Masa = new RooRealVar("Masa","Invariant Mass (GeV/c^{2})",0);
  RooRealVar* J_lifetime_sig_data = new RooRealVar("J_lifetime_sig_data","J/psi Lifetime significance for data",0);
  RooRealVar* J_lifetime_sig_mc = new RooRealVar("J_lifetime_sig_mc","J/psi Lifetime significance for mc",0);
  RooRealVar* abs_J_eta1 = new RooRealVar("abs_J_eta1","Absolute value of leading muon psedorapidity",0);
  RooRealVar* muon_thresh = new RooRealVar("muon_thresh","Muon threshold",0);
  RooArgSet* variables = new RooArgSet(*Masa,*J_lifetime_sig_data, *J_lifetime_sig_mc, *abs_J_eta1,*muon_thresh);
  variables->readFromFile("variables_default.txt","READ",cutset.c_str()) ;
  variables->Print("v") ;

  min_Masa = Masa->getMin();
  max_Masa = Masa->getMax();
  min_J_lifetime_sig_data = J_lifetime_sig_data->getMin();
  max_J_lifetime_sig_data = J_lifetime_sig_data->getMax();
  min_J_lifetime_sig_mc = J_lifetime_sig_mc->getMin();
  max_J_lifetime_sig_mc = J_lifetime_sig_mc->getMax();
  min_abs_J_eta1 = abs_J_eta1->getMin();
  max_abs_J_eta1 = abs_J_eta1->getMax();
  min_muon_thresh = muon_thresh->getMin();
  thecutset = cutset;

  Long64_t nevent = Tree->fChain->GetEntries(); //Numero de eventos
  cout << "Sample " << datasample << " has been chosen" << endl;
  cout << "Loaded correct sign: " << nevent << " candidates" << endl;

  //Load wrong sign:
  TChain *treews = new TChain("Tree");
  string textfilews = Form("%s_ws",datasample.c_str());
  ifstream ifstextfilews( textfilews.c_str() ); 
  while( getline( ifstextfilews, templine ) ) treews->Add(templine.c_str()); 
  newTree *Treews = new newTree(treews); 
  Long64_t neventws = Treews->fChain->GetEntries(); //Numero de eventos
  cout << "Loaded wrong sign: " << neventws << " candidates" << endl;

  RooRealVar *therun = new RooRealVar("therun", "run",0);
  RooRealVar *theevent = new RooRealVar("theevent", "event",0);
  RooDataSet* datawstmp = new RooDataSet("datawstmp","datawstmp",RooArgSet(*Masa,*therun,*theevent));
  RooDataSet* dataws = new RooDataSet("dataws","dataws",RooArgSet(*Masa));
  map<string,int> run_event;
  string thisrunevent;

  for (Long64_t i=0; i< neventws;i++) {
    Treews->GetEntry(i);

    if ( min_Masa > Treews->P_massC || max_Masa < Treews->P_massC ) continue;
    if ( Treews->J_lifetime/Treews->J_lifetimeError < min_J_lifetime_sig_data) continue;
    if ( Treews->J_lifetime/Treews->J_lifetimeError > max_J_lifetime_sig_data) continue;
    if ( min_abs_J_eta1 > TMath::Abs(Treews->J_eta1) || max_abs_J_eta1 < TMath::Abs(Treews->J_eta1) ) continue;
    if ( min_muon_thresh > Treews->J_pt1 || min_muon_thresh > Treews->J_pt2 ) continue;

    //if (i>5000) break;

    *Masa = Treews->P_massC;
    *therun = Treews->run;
    *theevent = Treews->event;
    datawstmp->add(RooArgSet(*Masa,*therun,*theevent));

    thisrunevent = Form("%d_%d",Treews->run,Treews->event);
    run_event[thisrunevent.c_str()]++;   

  }
  cout << "Wrong sign entries after selection " <<  datawstmp->numEntries() <<endl;  

  RooArgSet* set;
  for(Long64_t i = 0 ; i < datawstmp->numEntries(); i++){
    set = (RooArgSet*) datawstmp->get(i);
    *Masa = set->getRealValue(Masa->GetName());
    *therun = set->getRealValue(therun->GetName());
    *theevent = set->getRealValue(theevent->GetName());
    thisrunevent = Form("%d_%d", (int) therun->getVal(), (int) theevent->getVal()); 
    //cout << thisrunevent.c_str() << " " <<run_event[thisrunevent.c_str()] << " ";
    if ( run_event[thisrunevent.c_str()] > 300) { 
      //cout << endl; 
      continue; 
    }
    //cout << run_event[thisrunevent.c_str()] << endl;
    dataws->add(RooArgSet(*Masa));
  }

  delete datawstmp;
  cout << "Wrong sign entries after selection / number of candidates cut " <<  dataws->numEntries() <<endl;  

  run_event.clear();

  RooDataSet* datatmp = new RooDataSet("datatmp","datatmp",RooArgSet(*Masa,*therun,*theevent));
  RooDataSet* data = new RooDataSet("data","data",RooArgSet(*Masa));

  for (Long64_t i=0;i< nevent;i++) {
    Tree->GetEntry(i);

    if ( min_Masa > Tree->P_massC || max_Masa < Tree->P_massC ) continue;
    if ( Tree->J_lifetime/Tree->J_lifetimeError < min_J_lifetime_sig_data) continue;
    if ( Tree->J_lifetime/Tree->J_lifetimeError > max_J_lifetime_sig_data) continue;
    if ( min_abs_J_eta1 > TMath::Abs(Tree->J_eta1) || max_abs_J_eta1 < TMath::Abs(Tree->J_eta1) ) continue;
    if ( min_muon_thresh > Tree->J_pt1 || min_muon_thresh > Tree->J_pt2 ) continue;

    *Masa = Tree->P_massC;
    *therun = Tree->run;
    *theevent = Tree->event;
    datatmp->add(RooArgSet(*Masa,*therun,*theevent));

    thisrunevent = Form("%d_%d",Tree->run,Tree->event);
    run_event[thisrunevent.c_str()]++;   

  }

  cout << "Correct sign entries after selection " <<  datatmp->numEntries() <<endl;  

  for(Long64_t i = 0 ; i < datatmp->numEntries(); i++){
    set = (RooArgSet*) datatmp->get(i);
    *Masa = set->getRealValue(Masa->GetName());
    *therun = set->getRealValue(therun->GetName());
    *theevent = set->getRealValue(theevent->GetName());
    thisrunevent = Form("%d_%d", (int) therun->getVal(), (int) theevent->getVal()); 
    //cout << thisrunevent.c_str() << " " <<run_event[thisrunevent.c_str()] << " ";
    if ( run_event[thisrunevent.c_str()] > 300) { 
      //cout << endl; 
      continue; 
    }
    //cout << run_event[thisrunevent.c_str()] << endl;
    data->add(RooArgSet(*Masa));
  }

  delete datatmp;
  cout << "Correct sign entries after selection / number of candidates cut " <<  data->numEntries() <<endl;  

  Masa->setBins(100) ;
  RooDataHist* hist_bkg_mass = dataws->binnedClone() ;
  RooHistPdf* histpdf_bkg_mass = new RooHistPdf("histpdf_bkg_mass", "histpdf_bkg_mass", *Masa, *hist_bkg_mass, 0) ; //0 is for not interpolation, 2 for second order interpolation.
  RooHistPdf* histpdf_bkg_mass_pol1 = new RooHistPdf("histpdf_bkg_mass_pol1", "histpdf_bkg_mass_pol1", *Masa, *hist_bkg_mass, 1) ; //0 is for not interpolation, 2 for second order interpolation.
  RooHistPdf* histpdf_bkg_mass_pol2 = new RooHistPdf("histpdf_bkg_mass_pol2", "histpdf_bkg_mass_pol2", *Masa, *hist_bkg_mass, 2) ; //0 is for not interpolation, 2 for second order interpolation.

//  RooKeysPdf* keyspdf_bkg_mass = new RooKeysPdf("keyspdf_bkg_mass","keyspdf_bkg_mass",*Masa,*dataws,RooKeysPdf::NoMirror) ;
//  RooKeysPdf* keyspdf_bkg_mass_mirr = new RooKeysPdf("keyspdf_bkg_mass_mirr","keyspdf_bkg_mass_mirr",*Masa,*dataws,RooKeysPdf::MirrorBoth) ;
//  RooKeysPdf* keyspdf_bkg_mass_mirr2 = new RooKeysPdf("keyspdf_bkg_mass_mirr2","keyspdf_bkg_mass_mirr2",*Masa,*dataws,RooKeysPdf::MirrorBoth,2) ;
//  RooKeysPdf* keys2pdf_bkg_mass = new RooKeysPdf("keys2pdf_bkg_mass","keys2pdf_bkg_mass",*Masa,*dataws,RooKeysPdf::NoMirror,2) ;//increased bandwidth scale factor (promotes smoothness over detail preservation)  //considera demasiado los extremos

  // Plot unbinned data and histogram pdf overlaid
  RooPlot* framews = Masa->frame(Title("Wrong sign mass"),Bins(100),Range(min_Masa, max_Masa)) ;
//  data->plotOn(frame1,MarkerColor(kBlue));
  dataws->plotOn(framews) ;
  histpdf_bkg_mass->plotOn(framews,LineColor(kBlue)) ;
  histpdf_bkg_mass_pol1->plotOn(framews,LineColor(kRed)) ;
  histpdf_bkg_mass_pol2->plotOn(framews,LineColor(kGreen)) ;    
//  keyspdf_bkg_mass->plotOn(framews,LineColor(kRed)) ;
//  keyspdf_bkg_mass_mirr->plotOn(framews,LineColor(kMagenta)) ;
//  keyspdf_bkg_mass_mirr2->plotOn(framews,LineColor(kOrange)) ;
//  keys2pdf_bkg_mass->plotOn(framews,LineColor(kGreen)) ;        
  TCanvas* cws = new TCanvas("cws","cws",800,700);
  cws->SetFillColor(0);
  cws->cd(1) ; gPad->SetLeftMargin(0.15) ; 
  framews->Draw() ;
  cws->Print(Form("massws_%sR%s.eps",datasample.c_str(),thecutset.c_str() ));


//  delete Tree;
//  delete tree;

  Int_t dentries = data->numEntries();


 //Fine tunned:
  Float_t normsig = 2.*dentries*0.50 ; 
  Float_t normbkg = 2.*dentries*0.50 ;  

  RooRealVar* sidebandlow = new RooRealVar("sidebandlow", "sideband low bound",0);
  RooRealVar* sidebandhigh = new RooRealVar("sidebandhigh", "sideband high bound",0);

  RooRealVar* peak_yield = new RooRealVar("peak_yield", "yield signal peak", 0, normsig);
  RooRealVar* bkgd_yield = new RooRealVar("bkgd_yield", "yield of background", 0, normbkg);

  RooRealVar* peak_gaussian_mean = new RooRealVar("peak_gaussian_mean", "mean of gaussian for signal peak",0);
  RooRealVar* peak_gaussian_sigma = new RooRealVar("peak_gaussian_sigma", "sigma of gaussian for signalpeak",0);
  RooRealVar* peak_gaussian_sigma2 = new RooRealVar("peak_gaussian_sigma2", "sigma of gaussian for signalpeak2",0);
  RooRealVar* peak_gaussian_width = new RooRealVar("peak_gaussian_width", "width of gaussian for signalpeak",0);
  RooRealVar* f_gauss = new RooRealVar("f_gauss", "fraction btw gauss",0);
  RooRealVar* bkgd_poly_c1 = new RooRealVar("bkgd_poly_c1", "slope of background",0); 
  RooRealVar* bkgd_poly_c2 = new RooRealVar("bkgd_poly_c2", "curvture of background",0);
  RooRealVar* bkgd_poly_c3 = new RooRealVar("bkgd_poly_c3", "3rd background",0);

  RooArgSet* massparams = new RooArgSet(*peak_gaussian_mean,*peak_gaussian_sigma,*peak_gaussian_sigma2,*peak_gaussian_width,*f_gauss);
  massparams->add(RooArgSet(*bkgd_poly_c1,*bkgd_poly_c2,*bkgd_poly_c3,*sidebandlow,*sidebandhigh));

  //RooPolynomial* bkgd_poly = new RooPolynomial("bkgd_poly", " function for background",*Masa,RooArgList(*bkgd_poly_c1,*bkgd_poly_c2));
  //RooChebychev* bkgd_poly;
  RooAbsPdf* bkgd_poly = 0;

  RooRealVar *A =  new RooRealVar("A", "A",Masa->getMin()); 
  RooRealVar *B =  new RooRealVar("B", "B",Masa->getMax()); 
  RooFormulaVar* bkgd_poly_c0 = new RooFormulaVar("bkgd_poly_c0", "bkgd_poly_c0", 
                     "1./(@1-@0) - @2*(@1+@0)/2. - @3*(@1*@1+@1*@0+@0*@0)/3. + @4*(@1+@0)*(@1*@1+@0*@0)/4.",
                      RooArgList(*A,*B,*bkgd_poly_c1,*bkgd_poly_c2,*bkgd_poly_c3)); 
  if (bkgdeg==1) 
    bkgd_poly = new RooPolynomial("bkgd_poly", " function for background",*Masa,RooArgList(*bkgd_poly_c0,*bkgd_poly_c1),0);
  if (bkgdeg==2)
    bkgd_poly = new RooPolynomial("bkgd_poly", " function for background",*Masa,RooArgList(*bkgd_poly_c0,*bkgd_poly_c1,*bkgd_poly_c2),0);
  if (bkgdeg==3)
    bkgd_poly = new RooPolynomial("bkgd_poly", " function for background",*Masa,RooArgList(*bkgd_poly_c0,*bkgd_poly_c1,*bkgd_poly_c2,*bkgd_poly_c3),0);

  RooGaussian* peak_gaussian = 0;
  RooGaussian* peak_gaussian1 = 0;
  RooGaussian* peak_gaussian2 = 0;
  RooAddPdf* doublepeak_gaussian = 0;
  if(ngauss==1)  
    peak_gaussian = new RooGaussian("peak_gaussian", "gaussian for signal peak", *Masa, *peak_gaussian_mean, *peak_gaussian_sigma);
  if(ngauss==2){
    peak_gaussian1 = new RooGaussian("peak_gaussian1", "gaussian for signal peak1", *Masa, *peak_gaussian_mean, *peak_gaussian_sigma);
    peak_gaussian2 = new RooGaussian("peak_gaussian2", "gaussian for signal peak2", *Masa, *peak_gaussian_mean, *peak_gaussian_sigma2);
    doublepeak_gaussian = new RooAddPdf("doublepeak_gaussian", "suma de gausianas", RooArgList(*peak_gaussian1,*peak_gaussian2),*f_gauss);
  }

//  RooBreitWigner* peak_gaussian = new RooBreitWigner("peak_gaussian", "gaussian for signal peak", *Masa, *peak_gaussian_mean, *peak_gaussian_sigma);
//  RooVoigtian* peak_gaussian = new RooVoigtian("peak_gaussian", "gaussian for signal peak", *Masa, *peak_gaussian_mean, *peak_gaussian_width, *peak_gaussian_sigma );

  RooAddPdf* totalPdf = 0;
  if(ngauss==1) totalPdf = new RooAddPdf("totalPdf", "suma S+B", RooArgList(*peak_gaussian,*bkgd_poly),RooArgList(*peak_yield,*bkgd_yield));
  if(ngauss==2) totalPdf = new RooAddPdf("totalPdf", "suma S+B", RooArgList(*doublepeak_gaussian,*bkgd_poly),RooArgList(*peak_yield,*bkgd_yield));

  //Initialize mass parameters:
  massparams->readFromFile("massparams_default.txt","READ",fitset.c_str() );
  massparams->Print("v") ;

  if (bkgdeg==1) { bkgd_poly_c2->setVal(0.0); bkgd_poly_c3->setVal(0.0); bkgd_poly_c2->setConstant(kTRUE); bkgd_poly_c3->setConstant(kTRUE); }
  if (bkgdeg==2) { bkgd_poly_c3->setVal(0.0); bkgd_poly_c3->setConstant(kTRUE);}
  if (ngauss==1) {
    f_gauss->setVal(1.0);
    peak_gaussian_sigma2->setVal(0.0); //will go to lowest allowed value
  }
 
  Masa->setRange("R1",min_Masa,sidebandlow->getVal()) ; //mean_mass - 4.0*width_mass; //5.5;
  Masa->setRange("R2",sidebandhigh->getVal(),max_Masa) ; //mean_mass + 4.0*width_mass; //5.75;
  RooDataSet* sidebands = (RooDataSet*) data->reduce(Cut( Form("Masa< %f || Masa> %f",sidebandlow->getVal(),sidebandhigh->getVal())));

  // FIT BACKGROUND

  TCanvas *c1 = new TCanvas("c1",Form("Analisis para el Background Fit"),100,10,900,600);
  c1->SetFillColor(0);
  c1->Divide(1);

  RooFitResult* fitres = bkgd_poly->fitTo(*dataws,Minos(kTRUE),Save(kTRUE),NumCPU(ncpu));//,Range("R1,R2")
  fitres->Print("v");

  int itb=0;
  while (fitres->status()!=0 || fitres->covQual()!=3 
          || bkgd_poly_c1->getAsymErrorLo() == 0. || bkgd_poly_c1->getAsymErrorHi() == 0. 
          || (bkgdeg > 1 && bkgd_poly_c2->getAsymErrorLo() == 0.) || (bkgdeg >1 && bkgd_poly_c2->getAsymErrorHi() == 0.)  ){
   if (itb%2){
    fitres = bkgd_poly->fitTo(*dataws,Strategy(1),Minos(kTRUE),Save(kTRUE),NumCPU(ncpu));
    fitres->Print("v");
   }
   else{
    fitres = bkgd_poly->fitTo(*dataws,Strategy(2),Minos(kTRUE),Save(kTRUE),NumCPU(ncpu));//NumCPU(4)
    fitres->Print("v");
   }
   itb++;
   if (itb>2) break;
  }

  c1->cd();
  TString titulo_bkg = Form("#psi(2S) wrong sign");
  RooPlot* Mframe = Masa->frame(Range(min_Masa, max_Masa),Bins(binfit),Title(titulo_bkg));
  dataws->plotOn(Mframe);
  bkgd_poly->plotOn(Mframe); 
  Double_t chi2_bkg = Mframe->chiSquare(fitres->floatParsFinal().getSize()+1); //(dof)+norm
  cout << "The chi2/ndf of bkg is: " << chi2_bkg << endl;
  //RooHist* hresid = Mframetot->pullHist(); //->residHist() ;
  //bkgd_poly->plotOn(Mframe,Range(min_Masa, max_Masa),NormRange("R1,R2"),LineStyle(kDashed));
  bkgd_poly_c1->setPlotLabel("c1") ;
  bkgd_poly_c2->setPlotLabel("c2") ;
  bkgd_poly_c3->setPlotLabel("c3") ;
  bkgd_poly->paramOn(Mframe,Parameters(RooArgSet(*bkgd_poly_c1,*bkgd_poly_c2,*bkgd_poly_c3)), Layout(0.57,0.99,0.3) );
  Mframe->Draw(); 
  c1->Print(Form("massbkg_%sR%s.eps",datasample.c_str(),thecutset.c_str() ) );


  //Fixing to wrong sign fit!
  bkgd_poly_c1->setConstant(kTRUE);
  bkgd_poly_c2->setConstant(kTRUE);
  bkgd_poly_c3->setConstant(kTRUE);

  TCanvas *c2 = new TCanvas("c2","Signal Analysis",200,10,700,600);
  c2->SetFillColor(0);
  c2->Divide(1);
  c2->cd();


  RooFitResult* fittot = totalPdf->fitTo(*data,Strategy(1),Extended(),Minos(dominos),Save(kTRUE),NumCPU(ncpu));
  fittot->Print("v");


  int it=0;
  while (fittot->status()!=0 || fittot->covQual()!=3
         || !(peak_yield->hasAsymError()) 
         || !(bkgd_yield->hasAsymError()) 
         || !(peak_gaussian_sigma->hasAsymError()) ){
   if (it%2){
    fittot = totalPdf->fitTo(*data,Strategy(1),Extended(),Minos(dominos),Save(kTRUE),NumCPU(ncpu));
    fittot->Print("v");
   }
   else{
    fittot = totalPdf->fitTo(*data,Strategy(2),Extended(),Minos(dominos),Save(kTRUE),NumCPU(ncpu));//NumCPU(4)
    fittot->Print("v");
   }
   it++;
   if (it>2) break;
  }

  TString titulo = Form("#psi(2S)");
  RooPlot* Mframetot = Masa->frame(Range(min_Masa, max_Masa),Bins(binfit),Title(titulo));
  data->plotOn(Mframetot);
  totalPdf->plotOn(Mframetot);
  Double_t chi2_m = Mframetot->chiSquare(fittot->floatParsFinal().getSize()+1); //(dof)+norm
  cout << "The chi2/ndf is: " << chi2_m << endl;
  RooHist* hresid = Mframetot->pullHist(); //->residHist() ;
  if (ngauss==1) totalPdf->plotOn(Mframetot,Components(*peak_gaussian),LineStyle(7),LineColor(kRed));
  if (ngauss==2) totalPdf->plotOn(Mframetot,Components(*doublepeak_gaussian),LineStyle(7),LineColor(kRed));
  totalPdf->plotOn(Mframetot,Components(*bkgd_poly),LineStyle(kDotted),LineColor(kGray+2));

  //Caculando signficancia en ventana de masa:
  Double_t sigmaw2= f_gauss->getVal()*peak_gaussian_sigma->getVal()*peak_gaussian_sigma->getVal()
                   + (1.-f_gauss->getVal())*peak_gaussian_sigma2->getVal()*peak_gaussian_sigma2->getVal();
  Double_t sigmaw=TMath::Sqrt(sigmaw2);
  Double_t mn3sig = peak_gaussian_mean->getVal() - 3.0*sigmaw; //sidebandlow->getVal();//2.85; //peak_gaussian_mean.getVal() - 3*peak_gaussian_sigma.getVal();
  Double_t pl3sig = peak_gaussian_mean->getVal() + 3.0*sigmaw;//sidebandhigh->getVal();//3.3;  //peak_gaussian_mean.getVal() + 3*peak_gaussian_sigma.getVal();
  Masa->setRange("window",mn3sig,pl3sig) ;
  RooAbsReal* fracSigRange = 0;
  if (ngauss==1) fracSigRange = peak_gaussian->createIntegral(*Masa,*Masa,"window");
  if (ngauss==2) fracSigRange = doublepeak_gaussian->createIntegral(*Masa,*Masa,"window");
  RooAbsReal* fracBkgRange = bkgd_poly->createIntegral(*Masa,*Masa,"window");
  Double_t nsigWindow = peak_yield->getVal() *  fracSigRange->getVal();
  Double_t nbkgWindow = bkgd_yield->getVal() *  fracBkgRange->getVal();
  cout << "Eventos de señal en ventana de masa: " << nsigWindow << endl;
  cout << "Eventos de background en ventana de masa: " << nbkgWindow << endl;
  Double_t signif = nsigWindow/TMath::Sqrt(nsigWindow + nbkgWindow);

  TPaveText *tbox = new TPaveText(0.65,0.15,0.995,0.495,"BRNDC");
  tbox->AddText(Form("Events = %d",dentries));
  tbox->AddText(Form("N_{sig} = %5.1f #pm %5.1f",peak_yield->getVal(),peak_yield->getError() ) );
  tbox->AddText(Form("N_{bkg} = %5.1f #pm %5.1f",bkgd_yield->getVal(),bkgd_yield->getError() ) );
  tbox->AddText(Form("#mu = %5.5f #pm %5.5f",peak_gaussian_mean->getVal(),peak_gaussian_mean->getError() ) );
//  tbox->AddText(Form("#sigma = %5.5f #pm %5.5f",peak_gaussian_sigma->getVal(),peak_gaussian_sigma->getError() ) );
//  tbox->AddText(Form("w = %5.5f #pm %5.5f",peak_gaussian_width->getVal(),peak_gaussian_width->getError() ) );
  tbox->AddText(Form("#sigma_{1} = %5.5f #pm %5.5f",peak_gaussian_sigma->getVal(),peak_gaussian_sigma->getError() ) );
  if (ngauss==2) tbox->AddText(Form("#sigma_{2} = %5.5f #pm %5.5f",peak_gaussian_sigma2->getVal(),peak_gaussian_sigma2->getError() ) );
  if (ngauss==2) tbox->AddText(Form("f_{1} = %5.2f #pm %5.2f",f_gauss->getVal(),f_gauss->getError() ) );
  tbox->AddText(Form("S / #sqrt{S+B} = %5.2f",signif ) );
//  tbox->AddText(Form("#chi^{2} = %5.2f",chi2_m ) );
  tbox->SetFillColor(0);
  tbox->SetTextColor(1);
  tbox->SetTextSize(0.04);
  tbox->SetBorderSize(1);
  Mframetot->addObject(tbox);

  Mframetot->Draw(); 
  c2->Print(Form("massfit_%sR%s.eps",datasample.c_str(),thecutset.c_str() ));

  TCanvas *c3 = new TCanvas("c3",Form("Normalized Residuals"),100,10,900,300);
  c3->SetFillColor(0);
  c3->Divide(1);
  c3->cd();
  RooPlot* frameres = Masa->frame(Range(min_Masa,max_Masa) ) ;
  frameres->addPlotable(hresid,"PZ") ;
  frameres->Draw();
  c3->Print(Form("masspulls_%sR%s.eps",datasample.c_str(),thecutset.c_str() ));


  string ffitresults = Form("fitresults_%sR%s.dat",datasample.c_str(),thecutset.c_str());
  massparams->writeToFile(ffitresults.c_str());
  massparams->Print("v");
  
  Double_t signalpercent = peak_yield->getVal()/(bkgd_yield->getVal()+ peak_yield->getVal())*100.0;

  string fsalida2 = Form("results_%sR%s.dat",datasample.c_str(),thecutset.c_str());
//  ofstream salida2(fsalida2.c_str(),ios::app);
  ofstream salida2(fsalida2.c_str());

  salida2 << "Entries "  << 
             "Signal " << "SignalE " << 
             "Bkgnd " << "BkgndE " <<
             "Signal% " << "Signifcn " << 
             "Liklhood " << "FitStat " <<  "MtrxE " << 
             "Chi2 " <<
             "mean " << 
             "sigma " << 
//             "width " << endl;
             "sigma2 " << 
             "sigmaw " << 
             "f_1 " << endl;

  salida2 << dentries << " "  << 
             peak_yield->getVal() << " " << peak_yield->getError() << " " << 
             bkgd_yield->getVal() <<  " " << bkgd_yield->getError() << " " <<
             signalpercent << " " << signif << " " << 
             fittot->minNll() << " " << fittot->status() << " " <<  fittot->covQual() << " " << 
             chi2_m  <<  " " <<
             peak_gaussian_mean->getVal() << " " << 
             peak_gaussian_sigma->getVal() <<  " " << 
//             peak_gaussian_width->getVal() << endl;
             peak_gaussian_sigma2->getVal() << " " <<
             sigmaw << " " <<
             f_gauss->getVal() << endl;

  salida2.close();

  string catsalida = Form("cat %s",fsalida2.c_str()); 
  int rts = system(catsalida.c_str());
  if (rts) cout << "System command problem" << endl;

  RooWorkspace *w = new RooWorkspace("w","workspace") ;
  w->import(*totalPdf) ; //,Silence()
  w->import(*data) ; 
//  w->Print() ;
  w->writeToFile(Form("workspace_%sR%s.root",datasample.c_str(),thecutset.c_str())) ;


  /////////////////  Background subtraction //////////////////
  Double_t winlowcut = mn3sig; //3.6439; //
  Double_t winhighcut = pl3sig; // 3.72618; //
  Double_t winnsig = nsigWindow; // 10316.2; //
  Double_t winnbkg = nbkgWindow; //42986.2; //
  Double_t fbkgwin = winnbkg/(winnbkg + winnsig); 

  //Now pull MC:
  TChain *tree_mc = new TChain("Tree");
  tree_mc->Add(Form("/higgs/data/iheredia/.data/conversion/b-psi2s-p17p20.root"));
  TCut mcnocut   = Form("(P_massC>%f && P_massC<%f && TMath::Abs(J_eta1)>%f && TMath::Abs(J_eta1)<%f && J_pt1>%f && J_pt2>%f && P_mcmatch>0)", min_Masa, max_Masa, min_abs_J_eta1, max_abs_J_eta1,min_muon_thresh,min_muon_thresh);
  TCut datanocut = Form("(P_massC>%f && P_massC<%f && TMath::Abs(J_eta1)>%f && TMath::Abs(J_eta1)<%f && J_pt1>%f && J_pt2>%f && J_lifetime/J_lifetimeError >=%f )", min_Masa, max_Masa, min_abs_J_eta1, max_abs_J_eta1, min_muon_thresh, min_muon_thresh, min_J_lifetime_sig_data); 
  TCut signalwincut = datanocut && Form("(P_massC>%f && P_massC<%f)",winlowcut,winhighcut);
  TCut sidebandscut = datanocut && Form("(P_massC<%f || P_massC>%f)",winlowcut,winhighcut);
  TCut mccut = mcnocut ;
  //If needed, can apply a previous weight in the reweight by checking if weight branch exists and multiplying this weight in mccut.
  //static TString invalidBranch("w_J_pt");
  //TBranch* br = (TBranch*) tree_mc->GetListOfBranches()->FindObject(invalidBranch);
  //TCut mccut = br ? mcnocut : mcnocut * "w_J_pt";
  //delete br;

  //Do bkg substraction
  vector<TTree*> vtree;
  vtree.push_back(tree); //push back All data
  vtree.push_back(tree); //treews //tree //push back once All data for sidebands or push back Wrong-Sign
  vtree.push_back(tree_mc); //"compare to" MC or sidebands or whatever

  vector<string> hname;
  hname.push_back("Data_bkgsub"); //must label bkg subs.
  hname.push_back("MC"); //"compare to" name (sidebands or MC or whatever)

  vector<TCut> vcut;
  vcut.push_back(signalwincut); //data
  vcut.push_back(sidebandscut); //sidebandscut //signalwincut //sidebands or Wrong-sign (can take just the signal window in this case for better modeling)
  vcut.push_back(mccut); //"compare to" cut.
//  dobkgsub(vtree,hname, vcut,fbkgwin,Form("nobkgsub_%sR%s",datasample.c_str(),thecutset.c_str())); //compare using the selected data samples and MC sample without weights.
  if (!compareonly) weight(vtree,hname, vcut,fbkgwin,  Form("nobkgsub_%sR%s",datasample.c_str(),thecutset.c_str())); //weight using the selected data samples and MC sample.

  TCut mccutww = mcnocut * "w_J_pt"; 

  //Loop on all MC samples to compare:
  ifstream ifsin( "MCstoCompare" );
  string thesample;
  while( getline( ifsin, thesample ) ) {
    //if ( thesample.find("_ws") == string::npos ) continue;
    TChain *tree_mc_w = new TChain("Tree");
    tree_mc_w->Add(Form("%s.root",thesample.c_str())); //_w_J_pt
    //update trees:
    vtree.pop_back(); //might have been MC, but maybe not the one which I want to compare to w/ and w/o weights. 
    vtree.push_back(tree_mc_w); //"compare to" MC w/o weights
    vtree.push_back(tree_mc_w); //"compare to" MC weighted
    //update names:
    hname.pop_back(); //remove name for MC used to weight.
    hname.push_back("MC"); //"compare to" MC w/o weights
    hname.push_back("MC_weighted"); //"compare to" MC weighted
    //update cuts:
    vcut.pop_back(); //remove cut for MC used to weight. (additional step convienient in loop).
    vcut.push_back(mccut); //"compare to" cut w/o weights.
    vcut.push_back(mccutww); //"compare to" cut with weights.
    dobkgsub(vtree,hname,vcut,fbkgwin,Form("bkgsub_%s",thesample.c_str() ) );
    //recover state before loop.
    delete tree_mc_w;
    vtree.pop_back();
    hname.pop_back(); 
    vcut.pop_back();
  }

}



TH1F* bkgsubs(vector<TTree*> tree,vector<string> hname, vector<TCut> cut, Double_t frac, string fillfunc, string descr, string savename, Double_t inibin, Double_t finbin, Int_t bin, Int_t normalize,Bool_t doweight,string xdescr,Bool_t log, string addstr) {

  //uses first tree (cut[0] == signalwin) and second tree (cut[1] == sidebands) for subtraction.
  //then plot following trees (if exist) using cut[i>1].

  int tsize= tree.size();
  int nadd = 1;

  TCanvas *cbkgsub = new TCanvas("cbkgsub",Form("%s Distribution",descr.c_str()),0,10,700,500);
  gStyle->SetOptStat(1101);

  TH1F *h1,*h2,*h3,*h4,*h5;
  
  h1 = new TH1F(hname[0].c_str(),Form("%s",descr.c_str()),bin,inibin,finbin);
  if (tsize>1+nadd) h2 = new TH1F(hname[1].c_str(),Form("%s",descr.c_str()),bin,inibin,finbin);
  if (tsize>2+nadd) h3 = new TH1F(hname[2].c_str(),Form("%s",descr.c_str()),bin,inibin,finbin);
  if (tsize>3+nadd) h4 = new TH1F(hname[3].c_str(),Form("%s",descr.c_str()),bin,inibin,finbin);
  if (tsize>4+nadd) h5 = new TH1F(hname[4].c_str(),Form("%s",descr.c_str()),bin,inibin,finbin);

  vector<TH1F*> h;
  h.push_back(h1);
  h.push_back(h2);
  h.push_back(h3);
  h.push_back(h4);
  h.push_back(h5);

  Color_t colors[5]={kBlue,kRed,kBlack,kMagenta,kGreen};

  //Will build h[0] = h1 (bkg subtracted) by hand
  TH1F *g1  = new TH1F("g1","g1",bin,inibin,finbin);
  TH1F *g2  = new TH1F("g2","g2",bin,inibin,finbin);
  g1->Sumw2(); g2->Sumw2();
  tree[0]->Project("g1",fillfunc.c_str(),cut[0]);
  tree[1]->Project("g2",fillfunc.c_str(),cut[1]);
  Double_t ns = g1->GetEntries();
  Double_t nb = g2->GetEntries();
  Double_t Factor = -1*frac*ns/nb;
  h[0]->Add(g1,g2,1.,Factor);
  h[0]->SetLineColor(colors[0]);
//  h[0]->SetLineWidth(2);

  for (int i = 1; i < tsize - nadd; i++){
    tree[i+nadd]->Project(hname[i].c_str(),fillfunc.c_str(),cut[i+nadd]);
    h[i]->SetLineColor(colors[i]);
    h[i]->SetLineWidth(2);
    h[i]->SetLineStyle(i); 
    h[i]->Sumw2();
  }

  h1->GetXaxis()->SetTitle(Form("%s",xdescr.c_str()));

  h1->SetMarkerStyle(20);
  h1->SetMarkerColor(colors[0]);
  h1->Draw("p"); //p show marker
  for (int i = 1; i < tsize-nadd; i++) {
    h[i]->Draw("samesHISTE"); //sames to paint stats, HIST to force paint histogram, E to force paint errors.
  }
  
  gPad->Update();

  TPaveStats *st1,*st2,*st3,*st4,*st5;
  st1 = (TPaveStats*)h1->FindObject("stats");
  if (tsize>1+nadd) st2 = (TPaveStats*)h2->FindObject("stats");
  if (tsize>2+nadd) st3 = (TPaveStats*)h3->FindObject("stats");
  if (tsize>3+nadd) st4 = (TPaveStats*)h4->FindObject("stats");
  if (tsize>4+nadd) st5 = (TPaveStats*)h5->FindObject("stats");

  vector<TPaveStats*> st;
  st.push_back(st1);
  st.push_back(st2);
  st.push_back(st3);
  st.push_back(st4);
  st.push_back(st5);

  int j=0;
  for (int i = 0; i < tsize-nadd; i++){
    h[i]->GetListOfFunctions()->Remove(st[i]);
    h[i]->SetStats(0);
    st[i]->SetTextColor(colors[i]);

    st[i]->SetX1NDC(0.85);
    st[i]->SetX2NDC(0.99);

    st[i]->SetY1NDC(0.86 - j*0.15);
    st[i]->SetY2NDC(0.98 - j*0.15);
    j++;
  }

  for (int i = 0; i < tsize-nadd; i++) st[i]->Draw();

  Double_t themax = h[0]->GetMaximum();
  for (int i = 1; i < tsize-nadd; i++) {
    Double_t maxi = h[i]->GetMaximum();
    if (maxi>themax) themax = maxi;
  }
  Double_t themaxint = h[0]->Integral();
  //Let's normalize to data only for convinience.
/*
  for (int i = 1; i < tsize-nadd; i++) {
    Double_t maxiint = h[i]->Integral();
    if (maxiint>themaxint) themaxint = maxiint;
  }
*/
  
  if (normalize==1){ //normalize by maximum
    for (int i = 0; i < tsize-nadd; i++) {
      h[i]->Scale(themax/h[i]->GetMaximum());
    }
  }
  if (normalize==2){ //normalize by integral
    for (int i = 0; i < tsize-nadd; i++) {
      h[i]->Scale(themaxint/h[i]->Integral());
      //h[i]->Scale(1.0/h[i]->Integral());
    }
  }

  Double_t themaxforint = h[0]->GetMaximum();
  for (int i = 1; i < tsize-nadd; i++) {
    Double_t maxi = h[i]->GetMaximum();
    if (maxi>themaxforint) themaxforint = maxi;
  }
  Double_t binerror = h[0]->GetBinError(h[0]->GetMaximumBin());
  //	if (themax!=h[1]->GetMaximum()) binerror = binerror/2.0; //??
  for (int i = 0; i < tsize-nadd; i++) {
    if (normalize<2) h[i]->SetMaximum(themax*1.1+1.*binerror); 
    if (normalize==2) h[i]->SetMaximum(themaxforint*1.1); 
  }

  if (log) gPad->SetLogy();
  //else for (int i = 0; i < tsize-nadd; i++) h[i]->SetMinimum(0);
  
  cout << Form("\n %s:",descr.c_str()) << endl;
  if (normalize==1) cout << "Normalized by maximum" <<endl;
  if (normalize==2) cout << "Normalized by integral" <<endl;

  //Report kolmogorov-smirnov
  TPaveText *ptt = new TPaveText(0.55,0.4,0.79,0.5,"brNDC");
  ptt->SetFillColor(4000);//ptt->SetFillColor(29);
  ptt->SetTextAlign(12);
  ptt->SetTextFont(72);
  ptt->SetTextSize(0.03);
  ptt->SetTextColor(1);

  if (1){ //Report KS controlled by boolean variable.
    for (int i = 0; i < tsize-nadd -1 ; i++) {
      for (int l = i+1; l < tsize-nadd; l++) {

        Double_t kolmo = (h[i]->KolmogorovTest(h[l]));
        cout << "Kolmogorov " << i+1 << "-" << l+1 <<": " << kolmo <<endl;
        //cout << "Chi2_ndf Data-MCnoncorr = " << (hjpsipt1->Chi2Test(hjpsipt3,"CHI2/NDF")) << endl;
        ptt->AddText(Form("KS(%d - %d)[%%]: %2.2f",i+1,l+1,kolmo*100.));//hname[i].c_str(),hname[l].c_str() //%5.2e

      }//for l
    }//for i

    ptt->Draw();
  }

  TH1F *ww = 0;
  if (doweight){ //Now it's only working for J/psi pt distribution.
    const int binvariable = 14;
    Double_t cotas[binvariable+1] = {inibin,4.,5.,6.,7.,8.,9.,10.,11.,12.,14.,16.,20.,26.,finbin};
    //ww = new TH1F("ww","Weight",bin,inibin,finbin);
    ww = new TH1F("ww","Weight",binvariable,cotas);
    TH1F* h1new = (TH1F*)h[0]->Rebin(binvariable,"h1new",cotas);
    TH1F* h2new = (TH1F*)h[1]->Rebin(binvariable,"h2new",cotas);
    ww->Sumw2();
    //ww->Divide(h[0],h[1],1.,1.,"");
    ww->Divide(h1new,h2new,1.,1.,"");
  }
  else {
   delete ww;
  }

  cbkgsub->Print(Form("%s/%s.eps",addstr.c_str(),savename.c_str()));

  delete st1;
  delete h1;
  if (tsize>1+nadd) { delete st2; delete h2; }
  if (tsize>2+nadd) { delete st3; delete h3; }
  if (tsize>3+nadd) { delete st4; delete h4; }
  if (tsize>4+nadd) { delete st5; delete h5; }
  delete g1; delete g2;
  delete cbkgsub;
  delete ptt;

  if (doweight) return ww;
  else return 0;
}


void weight(vector<TTree*> tree,vector<string> hname, vector<TCut> cut, Double_t frac,string addstr){

  string sys1 = Form("if [ ! -d %s ] ; then mkdir %s ; fi",addstr.c_str(),addstr.c_str());
  int rts = system (sys1.c_str());
  //string sysclean = Form("if [ -d %s ] ; then \"echo cleaning\" ; rm %s/*eps ; fi",addstr.c_str(),addstr.c_str());
  //system (sysclean.c_str());
  if (rts) cout << "System command problem" << endl;

  //normalization: 1 (by max) , 2 (by area).

  int bns =50;
  int n = 2;
/*
  vector<TCut> cutforward;
  vector<TCut> cutcentral;
  TCut cutc[5];
  TCut cutf[5];

  int tsize= tree.size();
  for(int i=0; i< tsize ; i++) {
    cutc[i] = cut[i] && "TMath::Abs(J_eta1)<1.";
    cutf[i] = cut[i] && "TMath::Abs(J_eta1)>1.";
    cutcentral.push_back(cutc[i]);
    cutforward.push_back(cutf[i]);
  } 
  TH1F *w2 = new TH1F("w1","Weight",bns,2.5,30);
*/

  //repeso:
  TH1F *w1 = new TH1F("w1","Weight",bns,2.5,30);
  w1 = bkgsubs(tree,hname,cut,frac,"J_pt","J/#psi p_{T}", "J_pt_reweight" ,2.5,35., bns,n,true,"GeV/c",false,addstr);
  TCanvas *cw1nn = new TCanvas("cw1nn","Weight1nn",10,10,700,500);
  TCanvas *cw1 = new TCanvas("cw1","Weight1",10,10,700,500);
  cw1nn->cd();
  w1->Draw();

  //Ajustando:
  //TF1 *fw = new TF1("fw","pol5",ini,fin); //fw->SetNpx(500);
  //TF1 *fw = new TF1("fw","pol0(0)+expo(1)",inibin,finbin);
  const int nparams = 4;
  TF1 *fw1 = new TF1("fw1",turn,2.5,35.,nparams);
  fw1->SetLineWidth(4);
  fw1->SetLineColor(kBlue);    
  fw1->SetParameters(0.,-0.2,10.,3.5);
  //w->Fit("fw","F0");
  w1->Fit("fw1","F0","ep"); //Use minuit, do not plot, 
  w1->Fit("fw1","FE","ep"); //Use minuit, use minos, save.
  cw1nn->Print(Form("w_J_pt_notnorm_%sR%s.eps",thedataset.c_str(),thecutset.c_str()));

  //now will normalize curve to 1 in infinity.
  cw1->cd();
  Double_t params1[nparams];
  fw1->GetParameters(params1); 
  w1->Scale(1./(params1[3]-params1[0])); //used analitic limit.

  w1->Fit("fw1","F0","ep"); //Use minuit, do not plot, 
  w1->Fit("fw1","FE+","ep"); //Use minuit, use minos, save.
  fw1->GetParameters(params1); 
  cw1->Print(Form("w_J_pt_%sR%s.eps",thedataset.c_str(),thecutset.c_str()));
  fw1->SetParameter(0,0.);
  fw1->SetParameter(3,1.);
  ifstream ifsin( "MCstoReweight" );
  string thesample;
  while( getline( ifsin, thesample ) ) {
    //tree->Add(templine.c_str());
    ifstream ifsin2( thesample.c_str() );
    string thefile;
    TChain *treein = new TChain("Tree");
    //string filein="/higgs/data/iheredia/.data/conversion/b-psi2s-p17.root";
    while( getline( ifsin2, thefile ) ) treein->Add(thefile.c_str());
    //string fileout="b-psi2s-p17_w_J_pt.root";
    string fileout = Form("%sR%s.root",thesample.c_str(),thecutset.c_str()); //_w_J_pt
    TFile *fout = new TFile(fileout.c_str(),"recreate");
    TTree *treeout = (TTree*)treein->CloneTree(0); //Cloning original tree

    Double_t pt,eta1;
    treein->SetBranchAddress("J_pt",&pt);
    treein->SetBranchAddress("J_eta1",&eta1);

    Double_t weight;
    static TString invalidBranch("w_J_pt");
    TBranch* br = (TBranch*) treein->GetListOfBranches()->FindObject(invalidBranch);
    if (br) { 
      cout << "Branch " << invalidBranch << " exists..." << endl; 
      treein->SetBranchAddress(invalidBranch,&weight); 
    }
    else    { 
      cout << "Branch " << invalidBranch << " does not exist..." << endl; 
      treeout->Branch(invalidBranch,&weight,"weight/D"); 
    } 

    Double_t ptout;
    treeout->SetBranchAddress("J_pt",&ptout); //for some reason the value is not copied automatically.
    Double_t etaout1;
    treeout->SetBranchAddress("J_eta1",&etaout1); //for some reason the value is not copied automatically.
    Double_t newweight;
    if (br) {
      treeout->SetBranchAddress(invalidBranch,&newweight); //for some reason the value is not copied automatically.
    }

    Long64_t nentries = (Long64_t)treein->GetEntries();
    cout << "Entradas de MC a repesar: " << nentries << endl;

    bool first = true;
    for (Long64_t i=0;i<nentries;i++) {
      treein->GetEntry(i);
      if (! (i%20000) ) cout << i << " entries processed" << endl;
      Double_t ppt[1] = {pt};

      if ( thecutset.find("central") != string::npos) {
        if (first) cout << "Assigning weights for CENTRAL region..." <<endl;
        if ( TMath::Abs(eta1) < 1. ) weight = br ? weight * turn(ppt,params1) : turn(ppt,params1) ;
        else weight = br ? weight * 1. :  1. ;
        first = false;
        //cout << pt << " "<< TMath::Abs(eta1) << " " << turn(ppt,params1) << " " << weight << endl;
      }
      else if ( thecutset.find("forward") != string::npos) {
        //double oldweight = weight;
        if (first) cout << "Assigning weights for FORWARD region..." <<endl;
        if ( TMath::Abs(eta1) > 1. ) weight = br ? weight * turn(ppt,params1) : turn(ppt,params1) ;
        else weight = br ? weight * 1. :  1. ;
        first = false;
        //cout << pt << " "<< TMath::Abs(eta1) << " " << turn(ppt,params1) << " " << weight << " " << oldweight << endl;
      }
      else {
        if (first) cout << "Assigning weights in FULL region..." <<endl;
        weight = br ? weight * turn(ppt,params1) : turn(ppt,params1) ;
        first = false;
      }

      if (weight < 0.0 ) { //Just give a very small weight... not 0 since it generates segviolations.
       //cout << i << "  " << pt << "     " << w->FindBin(pt) << "   " << (w->GetBinCenter(w->FindBin(pt))) << "   " << (w->GetBinContent(w->FindBin(pt))) << endl;
       //weight = 0.1; //0;
       weight=1.0/double(nentries)/10.0;
      }
      ptout = pt;
      etaout1 = eta1;
      newweight = weight;
      treeout->Fill();
    }

    treeout->Write();
    delete br;
    delete treeout;
    delete treein;
    delete fout;
    cout << thesample.c_str() << " " << thecutset.c_str() << " reweight is done!" << endl;

  } //while thesample

  delete fw1;
  delete w1;
  delete cw1;
  delete cw1nn;

}


void  dobkgsub(vector<TTree*> tree,vector<string> hname, vector<TCut> cut, Double_t frac,string addstr){

  string sys1 = Form("if [ ! -d %s ] ; then mkdir %s ; fi",addstr.c_str(),addstr.c_str());
  int rts = system (sys1.c_str());
  //string sysclean = Form("if [ -d %s ] ; then \"echo cleaning\" ; rm %s/*eps ; fi",addstr.c_str(),addstr.c_str());
  //system (sysclean.c_str());
  if (rts) cout << "System command problem" << endl;

  //normalization: 1 (by max) , 2 (by area).

  int bns =50;
  int n = 2;



  //masses

  n=2;
  bkgsubs(tree,hname,cut,frac,"J_mass","J/#psi Invariant Mass", "Jmass" ,2.8,3.35,bns,n,0,"M(#mu^{+}#mu^{-}) (GeV)",false,addstr);
  bkgsubs(tree,hname,cut,frac,"P_massC","J/#psi#pi^{+}#pi^{-} Invariant Mass", "P_massC" ,min_Masa,max_Masa,bns,n,0,"M(#mu^{+}#mu^{-}#pi^{+}#pi^{-}) (GeV)",false,addstr);

  //Vertex chi2
  n=2; //0
  bkgsubs(tree,hname,cut,frac,"J_chi2","#chi^{2}_{vtx}(#mu^{+}#mu^{-})", "J_chi2" ,0,16, bns,n,0,"#chi^{2}_{vtx}(#mu^{+}#mu^{-})",false,addstr);
  bkgsubs(tree,hname,cut,frac,"P_chi2_J1","#chi^{2}_{vtx}(J/#psi#pi^{(lead)})", "P_chi2_J1" ,0,25, bns,n,0,"#chi^{2}_{vtx}(#mu^{+}#mu^{-}#pi^{(lead)})",false,addstr);
  bkgsubs(tree,hname,cut,frac,"P_chi2_J2","#chi^{2}_{vtx}(J/#psi#pi^{(trail)})", "P_chi2_J2" ,0,25, bns,n,0,"#chi^{2}_{vtx}(#mu^{+}#mu^{-}#pi^{(trail)})",false,addstr);


  //pt's
  n=2;
  bkgsubs(tree,hname,cut,frac,"J_pt","p_{T}(J/#psi)", "J_pt" ,2.5,35., bns,n,false,"p_{T}(#mu^{+}#mu^{-}) (GeV)",false,addstr);
  bkgsubs(tree,hname,cut,frac,"P_pt","p_{T}(J/#psi#pi^{+}#pi^{-})", "P_pt" ,3.,40., bns,n,0,"p_{T}(#mu^{+}#mu^{-}#pi^{+}#pi^{-}) (GeV)",false,addstr);
  bkgsubs(tree,hname,cut,frac,"J_pt1","p_{T}(#mu^{(lead)})", "J_pt1" ,1.5,22., bns,n,0,"p_{T}(#mu^{(lead)}) (GeV)",false,addstr);
  bkgsubs(tree,hname,cut,frac,"J_pt2","p_{T}(#mu^{(trail)})", "J_pt2" ,min_muon_thresh,12., bns,n,0,"p_{T}(#mu^{(trail)}) (GeV)",false,addstr);
  bkgsubs(tree,hname,cut,frac,"P_pt1","p_{T}(#pi^{(lead)})", "P_pt1" ,0.5,7, bns,n,0,"p_{T}(#pi^{(lead)}) (GeV)",false,addstr);
  bkgsubs(tree,hname,cut,frac,"P_pt2","p_{T}(#pi^{(trail)})", "P_pt2" ,0.5,3.5, bns,n,0,"p_{T}(#pi^{(trail)}) (GeV)",false,addstr);

  //other
  //"P_massPP"
  bkgsubs(tree,hname,cut,frac,"P_mass - J_mass - P_massPP","Q", "Q" ,0,0.72,bns,n,0,"M(#mu^{+}#mu^{-}#pi^{+}#pi^{-} - #mu^{+}#mu^{-} - #pi^{+}#pi^{-}) (GeV)",false,addstr);
  bkgsubs(tree,hname,cut,frac,"P_massPP","#pi^{+}#pi^{-} Invariant Mass", "P_massPP" ,0.139*2,1.,bns,n,0,"M(#pi^{+}#pi^{-}) (GeV)",false,addstr);
  bkgsubs(tree,hname,cut,frac,"TMath::Max(P_mjp1, P_mjp2)","max(m(J/#psi#pi^{#pm}))", "maxMJpi" ,3.25,3.95,bns,n,0,"max(M(#mu^{+}#mu^{-}#pi^{#pm})) (GeV)",false,addstr);
  bkgsubs(tree,hname,cut,frac,"TMath::Abs(P_eta1 + P_eta2)/2.0","#bar{#eta}(#pi^{#pm})", "etapi" ,0,2.5,bns,n,0,"#bar{#eta}(#pi^{#pm})",false,addstr);
  bkgsubs(tree,hname,cut,frac,"TMath::Max(P_dR1, P_dR2)","max(#DeltaR(#pi^{#pm}))", "maxDR" ,0,1.45,bns,n,0,"max(#DeltaR(#pi^{#pm}))",false,addstr);

/*
 //pv size, trks.
  n=0;
  bkgsubs(tree,hname,cut,frac,"pv_size","Number of tracks in Primary vertex", "pv_size",0,120,120,n,0,"N_{tracks}(PV)",true,addstr);
  bkgsubs(tree,hname,cut,frac,"nvrt","Number Primary vertices", "nvrt",0,20,20,n,0,"N_{PV}",true,addstr);
*/
}


