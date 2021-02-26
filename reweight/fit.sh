#!/bin/bash

sufix=""
if [ $1 ] ; then sufix="_""$1" ; fi

sourcedir='/higgs/data/iheredia/.data/conversion'

ls ${sourcedir}/run2a.root > run2a
ls ${sourcedir}/run2b1.root > run2b1
ls ${sourcedir}/run2b2.root > run2b2
ls ${sourcedir}/run2b3.root > run2b3
ls ${sourcedir}/run2b4.root > run2b4

ls ${sourcedir}/run2a_ws.root > run2a_ws
ls ${sourcedir}/run2b1_ws.root > run2b1_ws
ls ${sourcedir}/run2b2_ws.root > run2b2_ws
ls ${sourcedir}/run2b3_ws.root > run2b3_ws
ls ${sourcedir}/run2b4_ws.root > run2b4_ws

cat run2b[1-4] > run2b
cat run2a run2b[1-4] > run2

cat run2b[1-4]_ws > run2b_ws
cat run2a_ws run2b[1-4]_ws > run2_ws

runs="run2" #"run2a" #run2b1 run2b2 run2b3 run2b4 run2b run2"

fitset="psi2s"
for run in $runs ; do
 thedir="fit_${run}${sufix}" #g${ngauss}_b${nbkg}
 if [ ! -d $thedir ] ; then mkdir $thedir ; fi
 cd $thedir
 cp ../${run} .
 cp ../${run}_ws .
 cp ../variables_default.txt .
 cp ../massparams_default.txt .
 ln -fs ../newTree.h
 cp ../reweight.C .

 #REWEIGHT CENTRAL REGION:  
 
 #MC to add weights:
 if [ -f MCstoReweight ] ; then rm MCstoReweight ; fi
 for ws in "" "_ws" 
 do
  for x in "qcd-x3872-p17-5pt" "qcd-x3872-p17-7pt" "qcd-x3872-p17-10pt" "b-psi2s-p17p20" "b-psi2s-p17" "b-psi2s-p20" "b-x3872-p17" "b-x3872-p20" "psi2s-p17" "psi2s-p20" "x3872-p17" "x3872-p20"  
  do
    ls ${sourcedir}/${x}${ws}.root > ${x}${ws}
    echo "${x}${ws}" >> MCstoReweight
  done
 done

 region="psi2s_central"
 #MC to compare in this pass:
 if [ -f MCstoCompare ] ; then rm MCstoCompare ; fi
 for x in "b-psi2s-p17p20" "b-psi2s-p17" "b-psi2s-p20"
 do
  echo "${x}R${region}.root" > ${x}R${region}
  echo "${x}R${region}" >> MCstoCompare
 done

 ngauss=2
 nbkg=2
 fs="${run}R${region}_Fit${fitset}_g${ngauss}_b${nbkg}"
 root -l -b -q "reweight.C+(\"$run\",\"${region}\",\"${fitset}\",$ngauss,$nbkg )" | tee out_${fs} #-b -q

 #REWEIGHT FORWARD REGION:  

 #MC to add weights:
 if [ -f MCstoReweight ] ; then rm MCstoReweight ; fi
 for ws in "" "_ws" 
 do
  for x in "qcd-x3872-p17-5pt" "qcd-x3872-p17-7pt" "qcd-x3872-p17-10pt" "b-psi2s-p17p20" "b-psi2s-p17" "b-psi2s-p20" "b-x3872-p17" "b-x3872-p20" "psi2s-p17" "psi2s-p20" "x3872-p17" "x3872-p20"  
  do
    ls ./${x}${ws}R${region}.root > ${x}${ws}R${region}
    echo "${x}${ws}R${region}" >> MCstoReweight
  done
 done

 region="psi2s_forward"
 #MC to compare in this pass:
 if [ -f MCstoCompare ] ; then rm MCstoCompare ; fi
 for x in "b-psi2s-p17p20" "b-psi2s-p17" "b-psi2s-p20"
 do
  echo "${x}Rpsi2s_centralR${region}.root" > ${x}Rpsi2s_centralR${region}
  echo "${x}Rpsi2s_centralR${region}" >> MCstoCompare
 done

 ngauss=2
 nbkg=2
 fs="${run}Rpsi2s_centralR${region}_Fit${fitset}_g${ngauss}_b${nbkg}"
 root -l -b -q "reweight.C+(\"$run\",\"${region}\",\"${fitset}\",$ngauss,$nbkg )" | tee out_${fs} #-b -q


 #Exiting thedir.
 cd ..

done #runs

