#!/bin/bash

> setenv

cwd=`pwd`

gmx=`which gmx`
if [ $gmx ]; then
    printf "Found path to gromacs binary GMXDIR = ${gmx:0:-4}\n"
    printf "Press ENTER to continue or specify path \n>>> "
else
    printf "Enter path to gromacs binary i.e. \${GMXDIR}/gmx\n>>> "
fi
read rep; if [ $rep ]; then gmxdir=${rep}; else gmxdir=${gmx:0:-4}; fi
$gmxdir/gmx > tempo 2> /dev/null; G=`grep gmx tempo | wc -l`
if [ $G == 3 ]; then
    export GMXDIR=${gmxdir}; echo -e "export GMXDIR=${gmxdir}" >> setenv
else
    echo "Incorrect path!!! please try again"; exit
fi; rm tempo
echo -e "\n-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-\n"

python=`which python`
if [ $python ]; then
    printf "Found path to python binary PYTHONDIR = ${python:0:-7}\n"
    printf "Press ENTER to continue or specify path \n>>> "
else
    printf "Enter path to python binary i.e. \${PYTHONDIR}/python\n>>> "
fi
read rep; if [ $rep ]; then python=${rep}; else python=${python:0:-7}; fi
echo "print('3')" > tempi; $python/python tempi > tempo; P=`cat tempo`
if [ $P == 3 ]; then
    export PYTHONDIR=${python}; echo -e "export PYTHONDIR=${python}" >> setenv
else
    echo "Incorrect path!!! please try again"; exit
fi; rm tempi tempo
echo -e "\n-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-\n"

echo "export convpdb=${cwd}/bin/toolset/perl/convpdb.pl" >> setenv
export convpdb=${cwd}/bin/toolset/perl/convpdb.pl
if [ ! -e $convpdb ]; then
    printf "Installing MMTSB Toolset...\n---------------------------\n"
    cd bin
    ./install_mmtsb.sh
    cd ..
    if [ ! -e $convpdb ]; then
	echo "Something went wrong! MMTSB Toolset could not be installed. Exiting now..."
	echo -e "\n-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-\n"
	exit
    else
	printf "\n---------------------------\nsuccess"
	echo -e "\n-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-\n"
    fi
fi

nslots=`${PYTHONDIR}/python -c 'import multiprocessing as mp; print(mp.cpu_count())'`
if [ $nslots ]; then
    printf "Found $nslots cpu cores available for calculation \n"
    printf "Press ENTER to continue or specify number \n>>> "
else
    printf "Enter number of cpu cores available for calculation \n>>> "
fi
read rep; if [ $rep ]; then nslots=${rep}; fi
export NSLOTS=$nslots; echo -e "export NSLOTS=$nslots" >> setenv
echo -e "\n-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-\n"

source setenv
exit
