import math
import ROOT
from array import array
from math import sqrt

from ROOT import gSystem, gROOT
from ROOT import *
import sys

def getEffy(filepat) :
    
    filename = open(filepat, 'r')
    n1, n1e, n2, n2e, n3, n3e, n4, n4e = filename.readline().split()

    fVAL = [float(n1), float(n2), float(n3), float(n4)]
    fVALE = [float(n1e), float(n2e), float(n3e), float(n4e)]  
    #fVAL.append(float(n1val))
    #fVALE.append(float(n1valE))

    return fVAL, fVALE

def getyieldsBu(nbin,etabin,filepat,sys) :
    
    nbkgVAL = []
    nbkgVALE = []
    n1VAL = []
    n1VALE = []

    for mm in range(nbin):           
        filename = open(str(filepat)+'output_BuFit_ptbins_{}{}_{}.txt'.format(sys,etabin[mm],etabin[mm+1]), 'r')   
        n1val, n1valE, nbkgval, nbkgvalE = filename.readline().split()
        n1VAL.append(float(n1val))
        n1VALE.append(float(n1valE))
        nbkgVAL.append(float(nbkgval))
        nbkgVALE.append(float(nbkgvalE))  
        #print n1val
    
    return n1VAL, n1VALE

def save_CS(list_YieldBu,list_YieldBuE,EffyBu,EffyBuE,nbin,ptbin,Hname,Xname,Yname,legname,L,Le,bini,binf,isptH):

    Canv1 = TCanvas('Canv1', "", 900, 800 )
    Canv1.cd()
    gPad.SetBottomMargin( .1 )
    gPad.SetLeftMargin( .12)
    gPad.SetTopMargin( .05 )
    gPad.SetRightMargin(0.03)

    BRbu  = 0.001026
    BRbue = 0.000031
    BRJ = 0.05961   
    BRJe = 0.00033
    L = 189.7
    #( Nu/(DBpte*Eft) )*( 1/(L*BRbu*BRJ*2.0) );
       
    HistoP = TH1D("HistoP","HistoP",nbin,array('d',ptbin))
    HistoP2 = TH1D("HistoP2","HistoP2",nbin,array('d',ptbin))

    for jj in range(len(list_YieldBu)):
        #Debugging Prints
        #print list_YieldBu[jj]
        #print list_YieldBu/EffyBu[jj]
        #print "bin width:  ", ptbin[jj+1]-ptbin[jj]
        #R=(list_YieldBu[jj]/EffyBu[jj])*(1.0/(ptbin[jj+1]-ptbin[jj]))*( 1.0/(L*BRbu*BRJ*2.0))
        R=(list_YieldBu[jj]/EffyBu[jj])*(1.0/(ptbin[jj+1]-ptbin[jj]))*( 1.0/(L*BRbu*BRJ*2.0))*(1.0/1000.0)  #factor 1000 is to move to microbarns(#mub)
        #print R
        # These are the different contributions to the uncertainty estimation
        statistical_contribution_to_sqrt = (list_YieldBuE[jj]/list_YieldBu[jj])**2 
        efficiency_contribution_to_sqrt =  (EffyBuE[jj]/EffyBu[jj])**2
        
        
        ## RE is the statistical uncertainty
        RE = R * sqrt( statistical_contribution_to_sqrt)

        ## RE2 is the sum of all the uncertainties computed so far
        RE2 = sqrt( RE*RE + R*R*efficiency_contribution_to_sqrt)  

        #Fill histograms               
        HistoP.SetBinContent(jj+1,R)
        HistoP.SetBinError(jj+1,RE)
        HistoP2.SetBinContent(jj+1,R)
        HistoP2.SetBinError(jj+1,RE2)
        print "R: ", R, " +/- ", RE, " +/- ", RE2
        #print " "
        
      
    draw_histogram(HistoP,HistoP2,Hname,Xname,legname,Yname,L,Le,bini,binf,isptH)


def draw_histogram(HistoCS1,HistoCS2,hname,Xname,legname,Yname,L,Le,bini,binf,ispth):

    # pPb 5.0 TeV Results
    STevbin  = [5,10,15,20,25,30,60]
    STevH    = [0,149.54,43.81,11.37,4.25,0.85]
    STevHSTE = [0,10.65,3.62,1.49,0.87,0.14] #jus statistical  uncertainti (AN2013-322, Table 18)
    STevHSE  = [0,22.43,6.26,1.63,0.60,0.13] #jus systematical uncertainti (AN2013-322, Table 18)
    Histo5Tev = TH1D("Histo5Tev","Histo5Tev",6,array('d',STevbin))

    for jj in range(len(STevH)):
        Histo5Tev.SetBinContent(jj+1,STevH[jj])
        #Histo5Tev.SetBinError(jj+1,STevHSTE[jj])
        Histo5Tev.SetBinError(jj+1,sqrt(STevHSTE[jj]*STevHSTE[jj]+STevHSE[jj]*STevHSE[jj]) )

    Canv1 = TCanvas('Canv1', "", 900, 800 )
    Canv1.cd()
    gPad.SetBottomMargin( .11 ) 
    gPad.SetLeftMargin( .12)   
    gPad.SetTopMargin( .05 )   
    gPad.SetRightMargin(0.02)
    gPad.SetLogy();

    Histo5Tev.Draw("e1")
    Histo5Tev.SetStats(0)
    Histo5Tev.SetTitle("")
    Histo5Tev.SetMarkerSize(1)
    Histo5Tev.SetMarkerStyle(25)
    Histo5Tev.SetMarkerColor(2)
    Histo5Tev.SetLineColor(2)
    Histo5Tev.SetLineWidth(1)
    Histo5Tev.GetYaxis().CenterTitle(True)
    Histo5Tev.SetYTitle(Yname)
    Histo5Tev.GetXaxis().CenterTitle(True)
    #Histo5Tev.SetXTitle("p_{T}(J/#psi) (GeV)")
    Histo5Tev.SetXTitle(Xname)    
    Histo5Tev.SetTitleSize(35,"XY")
    Histo5Tev.SetLabelSize(30,"XY")
    Histo5Tev.SetTitleOffset(1.1,"Y")
    Histo5Tev.SetTitleOffset(1.0,"X")
    Histo5Tev.SetLabelFont(43,"XY")
    Histo5Tev.SetTitleFont(43,"XY")
    #Histo5Tev.SetMinimum(0.0)  ## Minimum for RAW is 0.0445 ## The usual is 0.15
    Histo5Tev.SetMaximum(1000.0)  ## Maximum for RAW is 0.15   ## For the rest 0.3

    #if ispth >= 1:
    #print "here"
    HistoCS1.Draw("same e1")
    HistoCS1.SetMarkerSize(1)
    HistoCS1.SetMarkerStyle(20)
    HistoCS1.SetMarkerColor(1)
    HistoCS1.SetLineColor(1)
    HistoCS1.SetLineWidth(1)

    HistoCS2.Draw("same e1")
    HistoCS2.SetMarkerSize(1)
    HistoCS2.SetMarkerStyle(20)
    HistoCS2.SetMarkerColor(1)
    HistoCS2.SetLineColor(1)
    HistoCS2.SetLineWidth(1)
        
    legend = TLegend(.45,.80,.85,.90)
    legend.SetFillStyle(0)
    legend.SetBorderSize(0)
    legend.SetFillColor(0)
    legend.SetTextSize(0.03)
    legend.AddEntry(HistoCS1,legname, "ep")
    #if ispth >= 1:
    legend.AddEntry(Histo5Tev,"pPb 5.02 TeV (PRL 116, 032301)", "ep")
    legend.Draw("same")

    line3 = TLine(bini,L,binf,L) 
    line3.SetLineStyle(2)
    line3.SetLineWidth(3)
    #line3.Draw("same") 

    tex1 = TLatex(0.98,0.96,"pPb 8.16 TeV");
    tex1.SetNDC()
    tex1.SetTextAlign(31)
    tex1.SetTextFont(42)
    tex1.SetTextSize(0.04)
    tex1.SetLineWidth(2)

    tex2 = TLatex(0.92,0.96,"#sqrt{s} = 8.16 TeV");
    tex2.SetNDC()
    tex2.SetTextAlign(31)
    tex2.SetTextFont(42)
    tex2.SetTextSize(0.04)
    tex2.SetLineWidth(2)

    tex3 = TLatex(0.12,0.96,"CMS")
    tex3.SetNDC()
    tex3.SetTextFont(61)
    tex3.SetTextSize(0.04)
    tex3.SetLineWidth(2)

    tex4 = TLatex(0.20,0.96,"Preliminary")
    tex4.SetNDC()
    tex4.SetTextFont(52)
    tex4.SetTextSize(0.03)
    tex4.SetLineWidth(2)

    tex1.Draw()
    #tex2.Draw()
    tex3.Draw()
    tex4.Draw()

    Canv1.Modified()
    gPad.Update()

    Canv1.SaveAs('plots/Histo_'+str(hname)+'.png')
    Canv1.SaveAs('plots/Histo_'+str(hname)+'.pdf')

    del Canv1
    return 0


def ptbins(Rc,RcE) :
    #gROOT.SetBatch(kTRUE);
    #gStyle.SetErrorX(0)

    nbin=4
    ptbin = [7,10,15,20,50]

    ## These are the nominal yields on pt bins
    list_YieldBu, list_YieldBuE = getyieldsBu(nbin,ptbin,"/cms/jhovanny/mytest/BuHI/plots_ptbins/","")
    print "yields Bu: ", list_YieldBu
    print " "

    EffyBu, EffyBuE  = getEffy("/cms/jhovanny/mytest/BuHI/MCanalysis/Effy/plotseffyBu/output_BuTotaleffy_pt.txt")
    print "Effy Bu: ", EffyBu

    save_CS(list_YieldBu,list_YieldBuE,EffyBu,EffyBuE,nbin,ptbin,"Bu_CS_pt","p_{T}(B^{+})","d#sigma/dp_{T}(#mub/GeV)","pPb 8.16 TeV", Rc, RcE, 7.0,50.0, 1)

def main() :
    gROOT.SetBatch(kTRUE);
    #gStyle.SetErrorX(0)

    BRbu  = 0.001026   # BR(B+ to JpsiK)  from PDG ( 1.026 +/- 0.031 )x10^-3
    BRbue = 0.000031
    BRJ   = 0.05961    # BR(J/psi to mumu)  from PDG ( 5.961 +/- 0.033 )%
    BRJe  = 0.00033
    L     = 189.7  #189.7 nb
    Eft   = 0.05 # OJO, PONER EL VALOR REAL
    EftE  = 0.001 # OJO, PONER EL VALOR REAL

    #( Nu/(Bpte*Eft) )*( 1/(L*BRbu*BRJ*2.0) );
    #The factor of 2 accounts for our choice of quoting the Bc cross section times branching fraction for positive charge only while N sig includes both Bp and Bn                                     
    #Ojo, el factor Bpte[i] es el ancho del bin     

    # First Total Value
    Filename = open("/cms/jhovanny/mytest/BuHI/MCanalysis/plots_ptbinsMC/output_BuFit_1_ptbins_7_50.txt", 'r')   
    N1val, N1valE =Filename.readline().split()
    N1val = float(N1val)
    N1valE = float(N1valE)

    #Rc = N1val, N1valE
    #Rc, RcE = Ratio(float(N2val), float(N2valE), float(N1val), float(N1valE))
    Rc = (N1val/Eft)*( 1.0/(L*BRbu*BRJ*2.0) )
    RcE = Rc*sqrt( (N1valE*N1valE/(N1val*N1val)) + (EftE*EftE/(Eft*Eft)) ) 
    #print "RT: ", Rc, " +/- ", RcE
    #print " "
    
    print "Now lest see pt bins"
    ptbins(Rc,RcE)
    print " "
  
    
if __name__ == '__main__':
    main()    
 
