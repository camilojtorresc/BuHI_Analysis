#!/bin/bash
sourcedir='/higgs/data/iheredia/.data/reweight_psi2s/fit_run2'

runs="run2" #"run2a" #run2b1 run2b2 run2b3 run2b4 run2b run2"
fitset="X3872"
w="Rpsi2s_centralRpsi2s_forward"

for run in $runs ; do
 thedir="compare_${run}_${fitset}" #g${ngauss}_b${nbkg}
 if [ ! -d $thedir ] ; then mkdir $thedir ; fi
 cd $thedir
 cp ../${run} .
 cp ../${run}_ws .
 cp ../variables_default.txt .
 cp ../massparams_default.txt .
 ln -fs ../newTree.h
 cp ../reweight.C .

 #MC to compare in this pass:
 region="X3872full"
 if [ -f MCstoCompare ] ; then rm MCstoCompare ; fi
 hadd -f0 AllXMCR${region}.root ${sourcedir}/b-x3872-p17${w}.root ${sourcedir}/b-x3872-p20${w}.root ${sourcedir}/x3872-p17${w}.root ${sourcedir}/x3872-p20${w}.root
 echo "AllXMCR${region}" >> MCstoCompare

 ngauss=1
 nbkg=2
 fs="${run}R${region}_Fit${fitset}_g${ngauss}_b${nbkg}"
 root -l -b -q "reweight.C+(\"$run\",\"${region}\",\"${fitset}\",$ngauss,$nbkg,1 )" | tee out_${fs} #-b -q

 #Exiting thedir.
 cd ..

done #runs

