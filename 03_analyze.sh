#!/bin/bash
# This script screens the frames from dcds at each cycle in all the runs
# Collects the positions of MG atoms and performs clustering analysis. 
# Uses MDAnalysis module
python=${PYTHONDIR}/python
# Uses MMTSB Toolset (https://github.com/mmtsb/toolset)
convpdb=bin/toolset/perl/convpdb.pl

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
    $python frames.py "${filename}"
    for ts in $(seq 100 100 1001)
    do
	echo "--- $ts"
	filename="frame.${ts}.pdb"
	head -n ${header} ${filename} > temp.pdb
	head -n $((rna+header)) ${filename} | tail -n ${rna} > temp1.pdb
	grep MG ${filename} >> temp1.pdb
	$convpdb -renumber 1 temp1.pdb > temp2.pdb
	cat temp.pdb temp2.pdb >> all.mg.positions.pdb
	rm -rf temp*pdb
    done
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
$python find_clusters.py
head -n ${header} ./clusters/rna.pdb > cluster_positions.pdb
grep ATOM ./clusters/rna.pdb > temp.pdb
for filename in ./clusters/mg*pdb
do
    grep ATOM ${filename} >> temp.pdb
done
$convpdb -renumber 1 temp.pdb >> cluster_positions.pdb
rm temp.pdb
cd ..

echo "Done!"
