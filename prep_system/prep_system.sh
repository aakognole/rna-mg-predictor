#!/bin/bash
cwd=`pwd`
convpdb=$cwd/../bin/toolset/perl/convpdb.pl
if [ ! -e $convpdb ]; then
    echo -e "\n-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-\n"
    printf "Installing MMTSB Toolset..."
    cd ../bin
    ./install_mmtsb.sh
    cd $cwd
    if [ ! -e $convpdb ]; then
	echo "Something went wrong! MMTSB Toolset could not be installed..."
	echo -e "\n-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-\n"
    else
	printf "success\n"
	echo -e "\n-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-\n"
    fi
fi

$convpdb -renumber 1 -setchain ' ' only_rna.pdb > temp.0.pdb

GMXDIR=${GMXDIR}
if [ ! -e ${GMXDIR}/gmx ]; then
    echo " Need \${GMXDIR} path for gromacs: \${GMXDIR}/gmx"
    read rep
    GMXDIR=$rep
fi
echo -e "\n-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-\n"
${GMXDIR}/gmx pdb2gmx -f temp.0.pdb -o gmx.pdb -p gmx.top -water tip3p -ter -merge all
echo -e "\n-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-\n"

l=`wc -l gmx.top | awk '{print $1}'`
tl=$((l-22))
hl=$((tl-25))
tail -n ${tl} gmx.top | head -n ${hl} > ../toppar/RNAA.itp
sed -i -e "s~RNA ~RNAA~g" ../toppar/RNAA.itp

$convpdb -center -orient -setchain ' ' gmx.pdb > temp.1.pdb

#file='temp.1.pdb'
coorstat (){
    file=$1
    xmin=`grep ATOM $file | awk '{printf "%.5f\n",$6}' | sort -n | head -n 1`
    xmax=`grep ATOM $file | awk '{printf "%.5f\n",$6}' | sort -n | tail -n 1`
    ymin=`grep ATOM $file | awk '{printf "%.5f\n",$7}' | sort -n | head -n 1`
    ymax=`grep ATOM $file | awk '{printf "%.5f\n",$7}' | sort -n | tail -n 1`
    zmin=`grep ATOM $file | awk '{printf "%.5f\n",$8}' | sort -n | head -n 1`
    zmax=`grep ATOM $file | awk '{printf "%.5f\n",$8}' | sort -n | tail -n 1`
    xcen=`echo $xmin $xmax | awk '{printf "%.5f",$1+(($2-$1)/2)}'`
    ycen=`echo $ymin $ymax | awk '{printf "%.5f",$1+(($2-$1)/2)}'`
    zcen=`echo $zmin $zmax | awk '{printf "%.5f",$1+(($2-$1)/2)}'`
    xbox=`echo $xmin $xmax | awk '{printf "%.5f",$2-$1}'`
    ybox=`echo $ymin $ymax | awk '{printf "%.5f",$2-$1}'`
    zbox=`echo $zmin $zmax | awk '{printf "%.5f",$2-$1}'`
    echo "-----------------coorstat------------------"
    echo "XMIN $xmin XMAX $xmax XCEN $xcen XBOX $xbox"
    echo "YMIN $ymin YMAX $ymax YCEN $ycen YBOX $ybox"
    echo "ZMIN $zmin ZMAX $zmax ZCEN $zcen ZBOX $zbox"
    echo "-------------------------------------------"
}

coorstat 'temp.1.pdb'
echo -e "\n-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-\n"
box=$xbox
if awk 'BEGIN {exit !('$ybox' >= '$box')}'; then box=$ybox; fi
if awk 'BEGIN {exit !('$zbox' >= '$box')}'; then box=$zbox; fi
#if [[ $ybox -gt $box ]]; then box=$ybox; fi
#if [[ $zbox -gt $box ]]; then box=$zbox; fi
box=`echo $box | awk '{printf "%.5f",$1+20.0}'`
echo "Enter boxsize you want to set (in Angstroms) (recommended = $box rounded up) -"
read rep
box=`echo $rep | awk '{printf "%.3f",$1}'`
echo -e "\n-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-\n"
${GMXDIR}/gmx solvate -cp temp.1.pdb -cs -o temp.1.1.pdb -box 9 > gmx.log 2> error.log
echo -e "\n-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-\n"
nsol=`grep OW temp.1.1.pdb | wc -l`
rm temp.1.1.pdb
nmg=`echo $rep | awk '{printf "%.0f",$1*$1*$1*6.02*0.0001*0.05}'`
npot=`echo $nmg | awk '{print $1*3}'`
qtot=`grep -i qtot ../toppar/RNAA.itp | tail -n 1 | awk '{print $11}'`
ncla=`echo $nmg $npot $qtot | awk '{print 2*$1+$2+$3}'`
nsol=`echo $nsol $nmg $npot $cla | awk '{print $1-2*($2+$3+$4)}'`
echo "Will be adding $nsol waters, $nmg Mg2+, $npot K+ and $ncla Cl- ions..."
echo "proceed? (y) / Enter manually (n):"
read rep
if [ "$rep" == "y" ]; then
    echo "Adding water and ions..."
else
    echo "Enter number of water molecules to add (for best guess check step3_pbcsetup.pdb from CHARMM-GUI output)"
    read nsol
    echo "Enter number of MG to add e.g. = 6.02 * 0.0001 * (X*Y*Z) * Conc_in_M "
    read nmg
    echo "Enter number of POT to add e.g. = 6.02 * 0.0001 * (X*Y*Z) * Conc_in_M"
    read npot
    echo "Enter number of CLA to add e.g. = n(mg)*2 + n(pot) + charge_on_rna "
    read ncla
fi
sed -e "s~<N_SOL>~$nsol~g" -e "s~<N_MG>~$nmg~g" -e "s~<N_POT>~$npot~g" -e "s~<N_CLA>~$ncla~g" system.top.tmpl > system.top
echo -e "\n-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-\n"
dx=`echo $xmin $xmax $xcen $box | awk '{printf "%.5f",($4/2.0)-$3}'`
dy=`echo $ymin $ymax $ycen $box | awk '{printf "%.5f",($4/2.0)-$3}'`
dz=`echo $zmin $zmax $zcen $box | awk '{printf "%.5f",($4/2.0)-$3}'`
echo "Translating coordinates: dx $dx dy $dy dz $dz"

echo "CRYST1  $box  $box  $box  90.00  90.00  90.00 P 1           1" > temp.2.pdb
$convpdb -translate $dx $dy $dz temp.1.pdb >> temp.2.pdb

coorstat 'temp.2.pdb'
echo -e "\n-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-\n"

${GMXDIR}/gmx insert-molecules -seed 973475 -f temp.2.pdb -ci ../toppar/charmm36.ff/mol/sol.pdb -o temp.3.pdb -nmol ${nsol} > gmx.log 2> error.log

${GMXDIR}/gmx insert-molecules -seed 951573 -f temp.3.pdb -ci ../toppar/charmm36.ff/mol/mg.pdb -o temp.4.pdb -nmol ${nmg} > gmx.log 2> error.log

${GMXDIR}/gmx insert-molecules -seed 926651 -f temp.4.pdb -ci ../toppar/charmm36.ff/mol/pot.pdb -o temp.5.pdb -nmol ${npot} > gmx.log 2> error.log

${GMXDIR}/gmx insert-molecules -seed 982928 -f temp.5.pdb -ci ../toppar/charmm36.ff/mol/cla.pdb -o temp.6.pdb -nmol ${ncla} > gmx.log 2> error.log

echo -e "\n-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-\n"
cp temp.6.pdb system.pdb
rm temp.*.pdb
rm setup.gmx.log
python=${PYTHONDIR}/python
echo "using python = $python"
$python ./nohyd.py
echo -e "\n-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-\n"
echo "                            Completed system preparation!"
echo -e "\n-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-\n"
wait
exit
