#!/usr/bin/env python

import os, re
import commands
import math, time
import sys

import math
import ROOT
from array import array

from ROOT import gSystem

print
print 'START'
print

#os.system( 'root -l -b -q "BuFit_ptbins.C(7.0,50.0)" ')
#return

n=4
ptbin = [7,10,15,20,50]

for mm in range(n):           
    os.system( ' root -l -b -q "BuFit_ptbins.C({},{})" '.format(ptbin[mm],ptbin[mm+1]) )

