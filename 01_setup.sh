#!/bin/bash

echo "Enter job name:"
read jobname
echo -e "\n-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-\n"

GMX=`which gmx`
echo "Enter path to gmx binary i.e. \${GMXDIR}/gmx (Found: ${GMX:0:-4})"
echo "Press ENTER to continue OR specify path"
printf ">>> "
read rep; if [ $rep ]; then GMXDIR=${rep}; else GMXDIR=${GMX:0:-4}; fi
export GMXDIR=$GMXDIR; echo "GMXDIR set to $GMXDIR"
echo -e "\n-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-\n"

PYTHON=`which python`
echo "Enter path to python binary i.e. \${PYTHONDIR}/python (Found: $PYTHON)"
echo "Press ENTER to continue OR specify path:"
printf ">>> "
read rep; if [ $rep ]; then PYTHONDIR=${rep}; else PYTHONDIR=${PYTHON:0:-7}; fi
export PYTHONDIR=${PYTHONDIR}; echo "PYTHONDIR set to ${PYTHONDIR}"
echo -e "\n-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-\n"

nslots=`nproc`
echo "Enter number of cpu core available (Found: nproc = $nslots)"
echo "Press ENTER to continue OR specify number:"
printf ">>> "
read rep; if [ $rep ]; then nslots=${rep}; fi
export NSLOTS=$nslots
echo -e "\n-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-\n"

echo "Enter path to your only rna pdb file (e.g. rna.pdb"
printf ">>> "
read onlyrna
cp $onlyrna ./prep_system/only_rna.pdb
echo -e "\n-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-\n"

# Prepare system top and pdb files
cd prep_system
sed -i -e "s~<jobname>~rna_${jobname}~g" system.top.tmpl
./prep_system.sh
cd ..

echo "Now generating production runs..."
totruns=5
run=1
echo -e "\nWant to keep BB and SC restraints during production MD? (y/n)"
read rest
if [ ${rest} = y ] ; then
    rest='true'
    cd openmm_template/restraints
    ${PYTHONDIR}/python set_restraints.py
    cd ../..
else
    rest='false'
    echo -e "\n-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-\n"
    echo "Only center of mass restraint on RNA will be used."
    CENTEROFMASS=`cat prep_system/rna_com.txt`
    echo -e "Found center of mass of rna at: $CENTEROFMASS (values in nanometer)"
    echo -e "Press ENTER to continue OR specify location in same format"
    printf ">>> "
    read rep; if [ $rep ]; then CENTEROFMASS=${rep}; fi
    echo -e "\n-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-\n"
fi

echo -e "\nDelete existing setup?(y) / Overwrite?(n) / Cancel?(cancel) "
read clean
if [ ${clean} = y ] || [ ${clean} = n ] ; then
    printf "Building $totruns runs for $jobname ..."
else
    echo "Cancelling the operation..."
    exit
fi

#CENTEROFMASS="5.0,5.0,5.0" # of your RNA in nanometer to holdit at the center of the box
NUMBEROFRNAATOMS=`grep -i qtot toppar/RNAA.itp | tail -n 1 | awk '{printf "%d",$1}'` # all atoms in RNA
CHARGEONRNA=`grep -i qtot toppar/RNAA.itp | tail -n 1 | awk '{printf "%d",$11}'`
NWAT=`tail -n 10 prep_system/system.top | grep SOL | awk '{printf "%d",$2}'`
NMG=`tail -n 10 prep_system/system.top | grep MG | awk '{printf "%d",$2}'`
NPOT=`tail -n 10 prep_system/system.top | grep POT | awk '{printf "%d",$2}'`
NCLA=`tail -n 10 prep_system/system.top | grep CLA | awk '{printf "%d",$2}'`
NPOT2=`echo $NMG $NPOT | awk '{printf "%d",2*$1+$2-2}'`
echo "charge_on_rna=${CHARGEONRNA}" > openmm_template/composition
echo "# can not remove all MG because ParmEd gives error reading topology" >> openmm_template/composition
echo "declare -a max_ion_comp=('$NWAT' '$NMG' '$NPOT' '$NCLA');" >> openmm_template/composition
echo "declare -a min_ion_comp=('$NWAT' '1' '$NPOT2' '$NCLA');" >> openmm_template/composition

while [ $run -le $totruns ];do

	mkdir -p ${run}
	if [ ${clean} = y ] ;then
	   rm -rf ./${run}/*
	fi
        cp -rfp prep_system/system.pdb ./${run}/solution.${run}.pdb
	cp -rfp prep_system/system.top ./${run}/solution.${run}.top

	cp -rfp openmm_template/gcmc.tmpl ${run}
        cp -rfp openmm_template/gcmc.equi.tmpl ${run}

	if [ $1 ]; then
	    if [ -e openmm_template/job.tmpl.${1}.sh ]; then
		sed -e "s/<run>/$run/g" -e "s/<job>/$jobname/g" \
		    -e "s/<rest>/$rest/g" -e "s~<GMXDIR>~${GMXDIR}~g" \
		    -e "s~<PYTHONDIR>~${PYTHONDIR}~g" openmm_template/job.tmpl.${2}.sh > ${run}/job.sh
	    else
		echo "There is no job.tmpl.\${1}.sh file to modify!"
		echo "Copy job.tmpl.sh file to job.tmpl.\${1}.sh and update."
	    fi
	else
	    sed -e "s/<run>/$run/g" -e "s/<job>/$jobname/g" \
		-e "s/<rest>/$rest/g" -e "s~<GMXDIR>~${GMXDIR}~g" \
		-e "s~<PYTHONDIR>~${PYTHONDIR}~g" openmm_template/job.tmpl.sh > ${run}/job.sh
	fi
	sed -e "s/<CENTEROFMASS>/${CENTEROFMASS}/g" \
	    -e "s/<NUMBEROFRNAATOMS>/${NUMBEROFRNAATOMS}/g" openmm_template/openmm_run_parmed_plumed.py > ${run}/openmm_run_parmed_plumed.py
	sed -e "s/<run>/$run/g" openmm_template/nohyd.tmpl.py > ${run}/nohyd.tmpl.py
	
	cp -rfp openmm_template/spe_calc.mdp ${run}/
	cp -rfp openmm_template/ratio.* ${run}/
	cp -rfp openmm_template/*.inp ${run}/
        cp -rfp openmm_template/omm*.py* ${run}/
        cp -rfp openmm_template/restraints ${run}/
	cp -rfp openmm_template/composition ${run}/
	run=$((run+1))
done

printf "finished!"
echo -e "\n-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-\n"

