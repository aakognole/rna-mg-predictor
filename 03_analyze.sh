#!/bin/bash
# This script screens the frames from dcds at each cycle in all the runs
# Collects the positions of MG atoms and performs clustering analysis. 
# Uses MDAnalysis module

cwd=`pwd`
source setenv

printf "\nEnter cutoff for clustering (Default: 2.5 Å) \n>>> "
read rep; if [ $rep ]; then cutoff=$rep; else cutoff=2.5; fi

# Number of non-hydrogen atoms in RNA
rna=`cat prep_system/natoms.txt`
# Number of lines for header in *nohyd.pdb files
header=2

mkdir -p cluster_analysis
mkdir -p cluster_analysis/max_cont_list
mkdir -p cluster_analysis/clusters

echo " Creating PDB with both RNA and MG positions..."
cd cluster_analysis
rm -rf all.mg.positions.pdb 2> /dev/null
for filename in ../*/solution.*.prod.*.nohyd.dcd
do
    echo $filename
    ${PYTHONDIR}/python frames.py "${filename}"
    for ts in $(seq 100 100 1001)
    do
	printf "... $ts"
	filename="frame.${ts}.pdb"
	head -n ${header} ${filename} > temp.pdb
	head -n $((rna+header)) ${filename} | tail -n ${rna} > temp1.pdb
	grep MG ${filename} >> temp1.pdb
	$convpdb -renumber 1 temp1.pdb > temp2.pdb
	cat temp.pdb temp2.pdb >> all.mg.positions.pdb
	rm -rf temp*pdb
    done
    printf "\n"
    rm frame*pdb
done

echo "Creating PDB with single RNA and combined MG positions..."
head -n $((rna+header)) all.mg.positions.pdb > temp.pdb
grep MG all.mg.positions.pdb >> temp.pdb
head -n ${header} all.mg.positions.pdb > mg.distribution.pdb
$convpdb -renumber 1 temp.pdb >> mg.distribution.pdb
rm temp.pdb

echo "Finding clusters of MG ions..."
rm -rf ./clusters/*pdb 2> /dev/null
${PYTHONDIR}/python find_clusters.py $cutoff
head -n ${header} ./clusters/rna.pdb > cluster_positions.pdb
grep ATOM ./clusters/rna.pdb > temp.pdb
for filename in ./clusters/mg*pdb
do
    grep ATOM ${filename} >> temp.pdb
done
$convpdb -renumber 1 temp.pdb >> cluster_positions.pdb
rm temp.pdb
cd ..

cp cluster_analysis/clusters/mg_predictor.pdb ./
if [ -e mg_predictor.pdb ]; then
    echo "Done! Check mg_predictor.pdb file results!"
else
    echo "Hmmm... something didn't work!"
fi

exit
