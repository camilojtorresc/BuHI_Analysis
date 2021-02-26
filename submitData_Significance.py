#!/usr/bin/env python

import os, re
#import commands
import math, time
import sys

import math
import ROOT
from array import array

from ROOT import gSystem

print ('')
print ('START')
print ('')

n1=6  # number of Kaon pT cuts
n2=10 # number of B+ Probability cuts66.9505
n3=6  # number of B+ life-time cuts

#bins
ptKbin = [0.5,0.6,0.7,0.8,0.9,1.0]
Bprobbin = [0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1]
pdlbin = [2,3,4,5,6,7]

# to test Kaon pT cuts
for i in range(n1):

    #to test B+ probability cuts
    for j in range(n2):

        #to test B+ life-time cuts
        for k in range(n3):
            os.system( ' root -l -b -q "BuFit_Sinificancia.C({},{},{})" '.format(ptKbin[i],Bprobbin[j],pdlbin[k]) ) # (ptk,prob,cts)



