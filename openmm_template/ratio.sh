#!/bin/bash
# This script collects information about gcmc steps and energies
# It also gives out acceptance ratio for fragments

grep SOL acceptance.*.dat | grep insert | sort -n -k 7 > sol.insertion.txt
grep SOL acceptance.*.dat | grep delete | sort -n -k 7 > sol.deletion.txt
si=`grep success sol.insertion.txt | wc -l`
sd=`grep success sol.deletion.txt | wc -l`
i_sol=`echo "$si $tds" | awk '{print $2-$1}'`
echo "Frag   in      out     change  ratio" | awk '{print $1,"\t",$2,"\t",$3,"\t",$4,"  ",$5}'
echo "-----------------------------------------------"
echo "SOL    $si       $sd       $i_sol   -" | awk '{print $1,"\t",$2,"\t",$3,"\t",$4,"\t",$5}'
echo "-----------------------------------------------"
grep MG acceptance.*.dat | grep insert | sort -n -k 7 > mg.insertion.txt 
grep MG acceptance.*.dat | grep delete | sort -n -k 7 > mg.deletion.txt 
tid=`cat mg.insertion.txt mg.deletion.txt | wc -l`
si=`grep success mg.insertion.txt | wc -l`
sd=`grep success mg.deletion.txt | wc -l`
i_mg=`echo "$si $sd" | awk '{print $2-$1}'`
echo "MG     $si       $sd       $i_mg    $tid" | awk '{print $1,"\t",$2,"\t",$3,"\t",$4,"\t",(($2+$3)/$5)}'
echo "-----------------------------------------------"
grep POT acceptance.*.dat | grep insert | sort -n -k 7 > pot.insertion.txt
grep POT acceptance.*.dat | grep delete | sort -n -k 7 > pot.deletion.txt
tid=`cat pot.insertion.txt pot.deletion.txt | wc -l`
si=`grep success pot.insertion.txt | wc -l`
sd=`grep success pot.deletion.txt | wc -l`
i_pot=`echo "$si $sd" | awk '{print $2-$1}'`
echo "POT    $si       $sd       $i_pot   $tid" | awk '{print $1,"\t",$2,"\t",$3,"\t",$4,"\t",(($2+$3)/$5)}'
echo "-----------------------------------------------"
grep CLA acceptance.*.dat | grep insert | sort -n -k 7 > cla.insertion.txt
grep CLA acceptance.*.dat | grep delete | sort -n -k 7 > cla.deletion.txt
si=`grep success cla.insertion.txt | wc -l`
sd=`grep success cla.deletion.txt | wc -l`
i_cla=`echo "$si $sd" | awk '{print $2-$1}'`
charge=`echo "$i_mg $i_pot $i_cla" | awk '{print -$1*2-$2+$3}'`
echo "CLA    $si       $sd       $i_cla    -   CHARGE= ${charge}" | awk '{print $1,"\t",$2,"\t",$3,"\t",$4,"\t",$5,"\t","\t",$6,$7}'
echo "--------------------------------------------------------------"
ttr=`cat equi.accept.*.dat | wc -l`
str=`grep success equi.accept.*.dat | wc -l`
grep success equi.accept.*.dat | sort -k 5 > all.trans-rot.txt
st=`grep success equi.accept.*.dat | grep translate | wc -l`
sr=`grep success equi.accept.*.dat | grep rotate | wc -l`
ratio=`echo "$str $ttr" | awk '{print $1/$2}'`
echo "translations: $st rotations: $sr ratio = $ratio"
echo "--------------------------------------------------------------"
