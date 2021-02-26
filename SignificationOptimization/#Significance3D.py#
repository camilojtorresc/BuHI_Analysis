import math
import ROOT
from array import array
from math import sqrt

from ROOT import gSystem
from ROOT import *
import sys

import random

def main() :

    #https://root.cern/doc/master/classTHistPainter.html#HP01d
    # en esta pagina uno puede ver las opciones soprtadas para dibujar histogramas 3D
    kpt = [0.5,0.6,0.7,0.8,0.9,1.0]
    Bpro = [0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1]
    ct = [2,3,4,5,6,7]
    print("began")

    #H2pt = TH2F("H2pt","H2pt",5,3.0,5.0,9,10,18)
    #H2pt = TH3F("H2pt","H2pt",6,5.0,10.0,10,1,10,6,2,7) #poner bien los bines
    H2pt = TH3F("H2pt","H2pt",6,0.5,1.0,10,1,10,6,2,7) 

    for mm in range(6):
        #print mm
        for nn in range(10):
            #print nn
            for zz in range(6):
                #print zz
                PTK=int(kpt[mm]*10.0)
                BPRO=int(Bpro[nn]*100.0)
                PDL=int(ct[zz])
                
                text_file = open('/cms/CamiloTorres/Bu2/plots_signi/output_Signif_ptk{}_Bpro{}_cts{}.txt'.format(PTK,BPRO,PDL), 'r')
                lines = text_file.readlines()
                Sig=float(lines[0])
                
                H2pt.SetBinContent(mm+1,nn+1,zz+1,Sig)
            

    canv = TCanvas('canv', "", 900, 800 )
    #canv.Divide(1)
    canv.cd()
    gPad.SetBottomMargin( .1 )
    gPad.SetLeftMargin( .1)
    gPad.SetTopMargin( .05 )
    gPad.SetRightMargin(0.1)
    #gPad.SetFillColor(10)

    #H2pt.Draw("")
    #H2pt.Draw("ISO")
    #H2pt.Draw("BOX")
    H2pt.Draw("BOX2Z")
    #H2pt.Draw("LEGO")
    H2pt.SetStats(0)
    H2pt.SetTitle("")
    H2pt.GetXaxis().CenterTitle(True)
    H2pt.SetXTitle("pT(k) (GeV)")
    H2pt.GetYaxis().CenterTitle(True)
    H2pt.SetYTitle("Bpro(%)")
    H2pt.GetZaxis().CenterTitle(True)
    H2pt.SetZTitle("ct/#sigma_{ct}")
    H2pt.SetTitleOffset(1.6,"X")
    H2pt.SetTitleOffset(1.8,"Y")
    H2pt.SetTitleOffset(1.3,"Z")
    #H2pt.GetYaxis().SetNdivisions(10004)
    H2pt.GetXaxis().SetNdivisions(505,1)
    H2pt.GetZaxis().SetNdivisions(505,1)
    """
    H2pt.SetMarkerSize(3)
    H2pt.SetMarkerStyle(23)
    H2pt.SetMarkerColor(2)
    H2pt.SetLineColor(2)
    H2pt.SetLineWidth(1)
    H2pt.SetTitleSize(35,"XY")
    H2pt.SetLabelSize(30,"XY")
    H2pt.SetTitleOffset(1.1,"Y")
    H2pt.SetTitleOffset(1.0,"X")
    H2pt.SetLabelFont(43,"XY")
    H2pt.SetTitleFont(43,"XY")"""

    
    canv.SaveAs('plots_optsigni/Histo2D.png')
    canv.SaveAs('plots_optsigni/Histo2D.pdf')
    canv.SaveAs('plots_optsigni/Histo2D.C')


if __name__ == '__main__':
    main()      
