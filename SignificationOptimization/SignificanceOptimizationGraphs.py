#!/usr/bin/env python

#import commands
import ROOT
import math, time
import sys

from ROOT import gSystem
from ROOT import *

from ROOT import TCanvas, TGraph
from ROOT import gROOT


import numpy as np
import matplotlib.pyplot as plt
import math
from math import sin
from array import array

ptKbin = [0.5,0.6,0.7,0.8,0.9,1.0]
Bprobbin = [0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1]
ctbin = [2.0,3.0,4.0,5.0,6.0,7.0]


#graph(inicio ptKbin, final ptKbin,inicio Bprobbin, final Bprobbin,inicio ctbin, final ctbin, nombre del array que va a variar)
def CreateCanvas(n1_i,n1_f,n2_i,n2_f,n3_i,n3_f,signif,array):
    
    W=800
    H=600
    Lum = "189.7 nb^{-1} (#it{p}Pb 8.16 TeV)"

    canv = TCanvas('canv', "", W ,H )
    y =np.array(signif)

    if(array==ptKbin):
        x =np.array(ptKbin)
        gr = TGraph(len(ptKbin),x,y)
        gr.GetXaxis().CenterTitle(True)
        gr.GetXaxis().SetTitle('pT(k) (GeV)')
    elif(array==Bprobbin):
        x =np.array(Bprobbin)
        gr = TGraph(len(Bprobbin),x,y)
        gr.GetXaxis().CenterTitle(True)
        gr.GetXaxis().SetTitle('Bpro(%)')
    elif(array==ctbin):
        x = np.array(ctbin)
        gr = TGraph(len(ctbin),x,y)
        gr.GetXaxis().CenterTitle(True)
        gr.GetXaxis().SetTitle("ct/#sigma_{ct}")
        
    gr.SetMarkerColor( 4 )
    gr.SetMarkerStyle( 21 )
    gr.GetYaxis().CenterTitle(True)
    gr.GetYaxis().SetTitle( 'Significance' )
    gr.SetTitle("");
    gr.Draw("ap")
    
    tex1 = TLatex(0.90,0.92,Lum)
    tex1.SetNDC()
    tex1.SetTextAlign(31)
    tex1.SetTextFont(42)
    tex1.SetTextSize(0.05) 
    tex1.SetLineWidth(2)
    
    tex2 = TLatex(0.10,0.92,"CMS")
    tex2.SetNDC()
    tex2.SetTextFont(61)
    tex2.SetTextSize(0.05) 
    tex2.SetLineWidth(2)
    
    tex3 = TLatex(0.19,0.92,"Preliminary")
    tex3.SetNDC()
    tex3.SetTextFont(52)
    tex3.SetTextSize(0.05)
    tex3.SetLineWidth(2)
    
    tex1.Draw()
    tex2.Draw()
    tex3.Draw()

    canv.SetGrid()

    if(array==ptKbin):
        canv.SaveAs('plots_optsigni/optimizedsignificance_fixedBpro{}_fixedct{}.png'.format(Bprobbin[n2_i],ctbin[n3_i]))
        print('Plot with fixed Bpro = {} and fixed ct = {} DONE'.format(Bprobbin[n2_i],ctbin[n3_i]))
    elif(array==Bprobbin):
        canv.SaveAs('plots_optsigni/optimizedsignificance_fixedPtk{}_fixedct{}.png'.format(ptKbin[n1_i],ctbin[n3_i]))
        print("Plot with fixed pT(K) = {} and fixed ct = {} DONE".format(ptKbin[n1_i],ctbin[n3_i]))
    else:  
        canv.SaveAs('plots_optsigni/optimizedsignificance_fixedPtk{}_fixedBpro{}.png'.format(ptKbin[n1_i],Bprobbin[n2_i]))
        print("Plot with fixed pT(K) = {} and fixed Bprob = {} DONE".format(ptKbin[n1_i],Bprobbin[n2_i]))
    

def Graph(n1_i,n1_f,n2_i,n2_f,n3_i,n3_f,array):

    
    signif = []
    
    for i in range(n1_i,n1_f):
        
        for j in range(n2_i,n2_f):
            
            for k in range(n3_i,n3_f):
                PTK=int(ptKbin[i]*10.0)
                BPRO=int(Bprobbin[j]*100.0)
                CT=int(ctbin[k])
                
                text_file = open('../plots_signi/output_Signif_ptk{}_Bpro{}_cts{}.txt'.format(PTK,BPRO,CT), 'r')
                lines = text_file.readlines()
                sig=float(lines[0])

                signif.append(sig)
                
                text_file.close()

    CreateCanvas(n1_i,n1_f,n2_i,n2_f,n3_i,n3_f,signif,array)


def main():

    #bins

    
    print(" ")
    print('START')
    print(" ")
    
    Graph(0,6,0,1,4,5,ptKbin)

    Graph(3,4,0,10,4,5,Bprobbin)

    Graph(3,4,0,1,0,6,ctbin)


if __name__ == '__main__':
    main()      
