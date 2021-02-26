#!/usr/bin/env python

import commands
import math, time
import sys

import math
from array import array

print
print 'START'
print

n1=6  # number of Kaon pT cuts
n2=10 # number of B+ Probability cuts
n3=6  # number of B+ life-time cuts

#bins
ptKbin = [0.5,0.6,0.7,0.8,0.9,1.0]
Bprobbin = [0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1]
pdlbin = [2,3,4,5,6,7]

BestSig = 0.0

# to test Kaon pT cuts
for i in range(n1):

    #to test B+ probability cuts
    for j in range(n2):

        #to test B+ life-time cuts
        for k in range(n3):
            PTK=int(ptKbin[i]*10.0)
            BPRO=int(Bprobbin[j]*100.0)
            PDL=int(pdlbin[k])
            
            text_file = open('plots_signi/output_Signif_ptk{}_Bpro{}_cts{}.txt'.format(PTK,BPRO,PDL), 'r')
            lines = text_file.readlines()
            sig=float(lines[0])
            
            if(BestSig<sig):
                BestSig=sig
                bestptcut=ptKbin[i]
                bestBprocut=Bprobbin[j]
                bestpdlcut = pdlbin[k]
                filename = 'output_Signif_ptk{}_Bpro{}_cts{}.txt'.format(PTK,BPRO,PDL)                
            text_file.close()

print 'File name with best significance: ',filename
print 'Best significance: ',BestSig
print ''
print "Best K(pT) cut : ",bestptcut
print "Best Bprob cut : ",bestBprocut
print "Best pdl   cut : ",bestpdlcut


