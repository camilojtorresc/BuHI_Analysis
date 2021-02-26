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

void Effy_Total()
{

 gStyle->SetOptTitle(0);
 gStyle->SetOptFit(0);
 gStyle->SetOptStat(0);

 // Read files
 const Int_t bin = 4;
 Double_t Bcpre[bin];
 Double_t Bcpree[bin];
 Double_t Bcreco[bin];
 Double_t Bcrecoe[bin];
 
 ifstream entrada1 ("plotseffyBu/output_Buprefiltereffy_pt.txt");
 if ( !entrada1 ) 
    {
      cout << "No se pudo abrir el archivo" << endl;
      exit( 1 );
    }
 entrada1 >> Bcpre[0] >> Bcpree[0] >> Bcpre[1] >> Bcpree[1] >> Bcpre[2] >> Bcpree[2] >> Bcpre[3] >> Bcpree[3] ;

 ifstream entrada2 ("plotseffyBu/output_Burecofilterffy_pt.txt");
 if ( !entrada2 ) 
    {
      cout << "No se pudo abrir el archivo" << endl;
      exit( 1 );
    }       
 entrada2 >> Bcreco[0] >> Bcrecoe[0] >> Bcreco[1] >> Bcrecoe[1] >> Bcreco[2] >> Bcrecoe[2] >> Bcreco[3] >> Bcrecoe[3] ;
 
 cout << Bcpre[3] << " " << Bcreco[3] << " " << Bcpre[3]*Bcreco[3] << endl;
 //cout << Bcpre[11] << " " << Bcpree[11] << " " << Bcpre[11]*Bcreco[11] << endl;

 //define los arreglos para el bineado
 Double_t Bpt[bin] = {8.5,12.5,17.5,35.0};
 Double_t Bpte[bin] = {1.5,2.5,2.5,15.0};


 //define los arreglos de Total Effi para usar en la cross-section
 Double_t Bte[bin] ;
 Double_t BteE[bin];
 
 RooRealVar x("x","x",7.0,50.0) ;
 RooRealVar y("y","y",0,200) ;
 RooDataSet dxy("dxy","dxy",RooArgSet(x,y),StoreError(RooArgSet(x,y))) ;  
 
 // Fill an example dataset with X,err(X),Y,err(Y) values
 for (int i=0 ; i<bin ; i++) 
   {
     x = Bpt[i];
     x.setError(Bpte[i]) ;
     
     y= Bcpre[i]*Bcreco[i];
     Double_t multe = (Bcpre[i]*Bcreco[i])*sqrt( (Bcpree[i]/Bcpre[i])*(Bcpree[i]/Bcpre[i]) + (Bcrecoe[i]/Bcreco[i])*(Bcrecoe[i]/Bcreco[i]) ) ;
     y.setError(multe);
     
     dxy.add(RooArgSet(x,y)) ;
     cout<<"Total efficiency for pt bins:  "<<Bcpre[i]*Bcreco[i]<<" +/- "<< multe <<endl;

     Bte[i] = Bcpre[i]*Bcreco[i];
     BteE[i] = multe;
     
   }
 dxy.Print("v");
 //return;
 
 ofstream salida ("plotseffyBu/output_BuTotaleffy_pt.txt");
 salida.is_open();
 for(int ii=0; ii<bin; ii++){
   salida << Bte[ii] << " " << BteE[ii] << " ";
 }
 salida.close(); 

 int H = 600;
 int W = 800;
 
 TCanvas *c3 = new TCanvas("c3","c3",50,50,W,H);
 //c3->Divide(1,2);
 c3->SetLeftMargin(0.125);
 c3->SetRightMargin(0.01);
 c3->SetTopMargin(0.09);
 c3->SetBottomMargin(0.14); 
 
 RooPlot* frame = x.frame(5.0,55,bin);
 dxy.plotOnXY(frame,YVar(y)) ;

 frame->SetYTitle("Total efficiency");
 //frame->SetYTitle("efficiency"); 
 frame->SetXTitle("p_{T}(J/#psi K^{+}) (GeV)");                                                                                                                    
 frame->SetLabelSize(0.05,"XY");
 frame->SetTitleSize(0.05,"XY");
 frame->GetYaxis()->CenterTitle();   
 frame->GetXaxis()->CenterTitle();
 frame->GetYaxis()->SetNdivisions(505,1);
 //frame->GetXaxis()->SetDecimals(2); 
 frame->GetXaxis()->SetNdivisions(504,1); 
 frame->SetTitleOffset(1.0,"X");
 frame->SetTitleOffset(0.9,"Y");
 frame->SetTitleSize(0.06,"XY");
 //frame->SetMinimum(0.2);
 //frame->SetMaximum(2200.0);
 frame->Draw();

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

 c3->Modified();
 
 c3->Print("plotseffyBu/Bueffy_Total.png");
 c3->Print("plotseffyBu/Bueffy_Total.pdf");

}
