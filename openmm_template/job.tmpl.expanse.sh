#!/bin/bash
#SBATCH --job-name="<job>_<run>"
#SBATCH --output="out_<job>_<run>"
#SBATCH --partition=gpu-shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --gpus=1
#SBATCH --mem=2000
#SBATCH -A mda215
#SBATCH --no-requeue
#SBATCH -t 48:00:00

module load cpu
module load gpu
module load slurm
module load cuda10.2/toolkit
module load gcc
module load openmpi
#SET the number of openmp threads
export OMP_NUM_THREADS=8
export LD_LIBRARY_PATH=/cm/local/apps/cuda/libs/current/lib64:$LD_LIBRARY_PATH

echo $HOSTNAME `date`

ln -s ../toppar/charmm36.ff ./
ln -s ../toppar ./

# Setup environment vairables to use OpenMM with CUDA and GROMACS
export NSLOTS=8
python=<PYTHONDIR>/python
export GMX_MAXBACKUP=-1
export GMXDIR=<GMXDIR>
if [ ! -e ${GMXDIR}/gmx ]; then
    echo "Could not find ${GMXDIR}/gmx !! Make sure you have GROMACS installed."
    exit
else
    echo "GMXDIR set to $GMXDIR"
fi
mdrun="${GMXDIR}/gmx mdrun"

run=<run>
rest='<rest>'

if [ `ls -l *state 2> /dev/null | wc -l` -gt 0 ]; then
    curr=`wc -l *state | sort -k 2 -n | tail -n 1 | awk '{split ($2,a,"."); print a[1]}'`
    prev=$((curr-1))
    cnt=1
else
    curr=0
    prev=0
    cnt=1
    echo "cycle PE-delete PE-equi PE-insert PE-equi PE-md-equil PE-md-prod" | awk '{print $1,"\t",$2,"\t",$3,"\t",$4,"\t",$5,"\t",$6,"\t",$7}' > energy.dat
    echo "cycle istat #_sol #_mg #_pot #_cla charge" | awk '{print $1,"\t",$2,"\t",$3,"\t",$4,"\t",$5,"\t",$6,"\t",$7}' > checkpoint.dat
    touch gmx.log
fi

dr=0 ##cushion between simulation box and gcmc box
simulation_box_size=`cat solution.${run}.pdb | grep 'CRYST' | awk '{printf "%5f",$2}'`
echo "simulation box size: $simulation_box_size"
gcmc_box_size=`echo ${simulation_box_size} ${dr} | awk '{printf "%5f",($1-$2)}'`
echo "gcmc box size: $gcmc_box_size"

ion_status=1 # 0: more negative muex, i.e. delete  1: less negative muex, i.e insert

declare -a fragnamelist=('sol' 'mg' 'pot' 'cla');
declare -a resid_list=('SOL' 'MG' 'POT' 'CLA');
declare -a frag_charge=('0' '2' '1' '-1');

declare -a muexlist=('-5.6' '-437.38' '-70.51' '-81.26');
declare -a lhfe_list=('-5.6' '-837.38' '30.51' '-81.26');
declare -a uhfe_list=('-5.6' '-37.38' '-170.51' '-81.26');

source composition
#charge_on_rna=-76
## can not remove all MG because ParmEd gives error reading topology 
#declare -a max_ion_comp=('30085' '30' '90' '74');
#declare -a min_ion_comp=('30085' '1' '148' '74');

max_ions=${max_ion_comp[1]}
min_ions=${min_ion_comp[1]}

farraysize=${#fragnamelist[@]}
farraysize=$((farraysize-1))
echo "farray size is: ${farraysize}"

centerx=`echo ${simulation_box_size} | awk '{printf "%5f",$1/2}'`
centery=`echo ${simulation_box_size} | awk '{printf "%5f",$1/2}'`
centerz=`echo ${simulation_box_size} | awk '{printf "%5f",$1/2}'`
echo "centerx: $centerx"

while [ $curr -le 10 ];
do
        timer=`date`
	echo "Running cycle = ${curr} at ${timer}"
	if [ $curr -eq 0 ];then

	   iptop=solution.${run}.top
	   ippdb=solution.${run}.pdb
           x=0
           while [ $x -le  $farraysize ];do
                 ifrags[$x]=${max_ion_comp[$x]}
                 maxfrags[$x]=${max_ion_comp[$x]}
                 minfrags[$x]=${max_ion_comp[$x]}
                 x=$((x + 1))
           done
           x=0
           charge=${charge_on_rna}
           while [ $x -le  $farraysize ];do
                 charge=`echo ${charge} ${frag_charge[$x]} ${ifrags[$x]} | awk '{print $1+$2*$3}'`
                 x=$((x + 1))
           done
           #charge=`echo "${charge_on_rna} ${ifrags[1]} ${ifrags[2]} ${ifrags[3]}" | awk '{print $1+($2*2)+($3*1)-($4*1)}'`

           ${GMXDIR}/gmx grompp -f spe_calc.mdp -o spe.calc.tpr -c ${ippdb} -p ${iptop} >> gmx.log
           $mdrun -nt $NSLOTS -deffnm spe.calc -s spe.calc.tpr -rerun ${ippdb} >> gmx.log
           PE_md_prod=`grep -A 3 'Coulomb-14' spe.calc.log | tail -n 1 | awk '{printf "%d",$1}'`
	   echo "${curr}-${cnt} ------0 ------0 ------0 ------0 ------0 ${PE_md_prod}" | awk '{print $1,"\t",$2,"\t",$3,"\t",$4,"\t",$5,"\t",$6,"\t",$7}' >> energy.dat
	   PE_swap0="------0"; PE_equi0="------0"; PE_swap="------0"; PE_equi="------0"
	else
           if [ $cnt == 1 ]; then
               source $curr.state
           fi
	   iptop=solution.${run}.gc.${prev}.top
	   ippdb=solution.${run}.prod.${prev}.pdb

           simulation_box_size=`cat solution.${run}.prod.${prev}.pdb | grep 'CRYST' | awk '{printf "%5f",$2}'`
           gcmc_box_size=`echo ${simulation_box_size} ${dr} | awk '{printf "%5f",($1-$2)}'`
           centerx=`echo ${simulation_box_size} | awk '{printf "%5f",$1/2}'`
           centery=`echo ${simulation_box_size} | awk '{printf "%5f",$1/2}'`
           centerz=`echo ${simulation_box_size} | awk '{printf "%5f",$1/2}'`

	fi

	optop=solution.${run}.gc.${curr}.top
	oppdb=solution.${run}.gc.${curr}.pdb

	x=0
        muex=' '
        maxf=' '
        minf=' '
	while [ $x -le  $farraysize ];do
              muex=$muex' '${muexlist[$x]}
              maxf=$maxf' '${maxfrags[$x]}
              minf=$minf' '${minfrags[$x]}
              x=$((x + 1))
	done

	sed  -e "s~<iptop>~$iptop~g" -e "s~<ippdb>~$ippdb~g" \
	     -e "s/<centerx>/$centerx/g" -e "s/<centery>/$centery/g" -e "s/<centerz>/$centerz/g" \
             -e "s/<simulation_box_size>/$simulation_box_size/g" \
	     -e "s/<gcmc_box_size>/$gcmc_box_size/g" \
	     -e "s/<fragmuex>/$muex/g" -e "s/<maxfrags>/$maxf/g" -e "s/<minfrags>/$minf/g" \
	     -e "s~<optop>~$optop~g" -e "s~<oppdb>~$oppdb~g" gcmc.tmpl > gcmc.${curr}.${cnt}.inp

	if [ $curr == 0 ]
	then
	   cp -rf ${iptop} ${optop}
	   cp -rf ${ippdb} ${oppdb}
	else
	   sed -e "s/<cycle>/$curr/g" ratio.tmpl > ratio.${curr}.sh
	   if [ ${cnt} != 1 ]; then
              cp solution.${run}.gc.${curr}.top temp_gc_top
              cp solution.${run}.gc.${curr}.pdb temp_gc_pdb
              sed  -e "s~<iptop>~temp_gc_top~g" -e "s~<ippdb>~temp_gc_pdb~g" \
                   -e "s/<centerx>/$centerx/g" -e "s/<centery>/$centery/g" -e "s/<centerz>/$centerz/g" \
                   -e "s/<simulation_box_size>/$simulation_box_size/g" \
                   -e "s/<gcmc_box_size>/$gcmc_box_size/g" \
                   -e "s/<fragmuex>/$muex/g" -e "s/<maxfrags>/$maxf/g" -e "s/<minfrags>/$minf/g" \
                   -e "s~<optop>~$optop~g" -e "s~<oppdb>~$oppdb~g" gcmc.tmpl > gcmc.${curr}.${cnt}.inp
	      ../bin/mod_swap_gcmc gcmc.${curr}.${cnt}.inp gcmc.${curr}.${cnt}.out -nt 1 >> acceptance.${curr}.dat
	   else
	      ../bin/mod_swap_gcmc gcmc.${curr}.${cnt}.inp gcmc.${curr}.${cnt}.out -nt 1 >> acceptance.${curr}.dat
           fi

           ${GMXDIR}/gmx grompp -f spe_calc.mdp -o spe.calc.tpr -c ${oppdb} -p ${optop} >> gmx.log
           $mdrun -nt $NSLOTS -deffnm spe.calc -s spe.calc.tpr -rerun ${oppdb} >> gmx.log
           PE_swap=`grep -A 3 'Coulomb-14' spe.calc.log | tail -n 1 | awk '{printf "%d",$1}'`

	   cp solution.${run}.gc.${curr}.top temp_gc_top
	   cp solution.${run}.gc.${curr}.pdb temp_gc_pdb
           sed  -e "s~<iptop>~temp_gc_top~g" -e "s~<ippdb>~temp_gc_pdb~g" \
                -e "s/<centerx>/$centerx/g" -e "s/<centery>/$centery/g" -e "s/<centerz>/$centerz/g" \
                -e "s/<simulation_box_size>/$simulation_box_size/g" \
                -e "s/<gcmc_box_size>/$gcmc_box_size/g" \
                -e "s/<fragmuex>/$muex/g" -e "s/<maxfrags>/$maxf/g" -e "s/<minfrags>/$minf/g" \
                -e "s~<optop>~$optop~g" -e "s~<oppdb>~$oppdb~g" gcmc.equi.tmpl > gcmc.equi.${curr}.${cnt}.inp
	   ../bin/mod_equi_gcmc gcmc.equi.${curr}.${cnt}.inp gcmc.equi.${curr}.${cnt}.out -nt $NSLOTS >> equi.accept.${curr}.dat

           ${GMXDIR}/gmx grompp -f spe_calc.mdp -o spe.calc.tpr -c ${oppdb} -p ${optop} >> gmx.log
           $mdrun -nt $NSLOTS -deffnm spe.calc -s spe.calc.tpr -rerun ${oppdb} >> gmx.log
           PE_equi=`grep -A 3 'Coulomb-14' spe.calc.log | tail -n 1 | awk '{printf "%d",$1}'`

	   x=0
           while [ $x -le  $farraysize ];do
              ifrags[$x]=`grep ${resid_list[$x]} ${optop} | tail -n 1 | awk '{print $2}'`
                 x=$((x + 1))
           done
           x=0
           charge=${charge_on_rna}
           while [ $x -le  $farraysize ];do
                 charge=`echo $charge ${frag_charge[$x]} ${ifrags[$x]} | awk '{print $1+$2*$3}'`
                 x=$((x + 1))
           done
	   #charge=`echo "${charge_on_rna} ${ifrags[1]} ${ifrags[2]} ${ifrags[3]}" | awk '{print $1+($2*2)+($3*1)-($4*1)}'`
	fi

	status=1
	if [ ${charge} == 0 ]; then
	   if [ ${ifrags[1]} == ${max_ions} ] && [ ${curr} != 0 ]; then
	      opgro=solution.${run}.gc.${curr}.gro
              eqrst=solution.${run}.equil.${curr}.rst
              eqdcd=solution.${run}.equil.${curr}.dcd
              eqpdb=solution.${run}.equil.${curr}.pdb
              prodrst=solution.${run}.prod.${curr}.rst
	      proddcd=solution.${run}.prod.${curr}.dcd
              prodpdb=solution.${run}.prod.${curr}.pdb

              #Convert pdb to gro for ParmEd
	      ${GMXDIR}/gmx editconf -f ${oppdb} -o ${opgro} >> gmx.log

	      #Run minimization and equilibration
              $python -u openmm_run_parmed_plumed.py -i mini_equi.inp -p ${optop} -c ${opgro} -orst ${eqrst} -odcd ${eqdcd} > solution.${run}.equil.${curr}.out
	      rm -f cv_ZZZ.dat center.dat                                                                                                                                 
	      echo "TITLE     <job>_<run>_${curr}" > temp.pdb
              echo "REMARK    PDB after equilibration" >> temp.pdb
              grep "CRYST1" output.pdb >> temp.pdb
              echo "MODEL        1" >> temp.pdb
              #Using MMTSB PDB manipulation script
              ../bin/toolset/perl/convpdb.pl output.pdb | sed -e "s/HETATM/ATOM  /g" -e "s/END/ENDMDL/g" >> temp.pdb
              mv temp.pdb ${eqpdb}
              rm -rf output.pdb

              ${GMXDIR}/gmx grompp -f spe_calc.mdp -o spe.calc.tpr -c ${eqpdb} -p ${optop} >> gmx.log
              $mdrun -nt $NSLOTS -deffnm spe.calc -s spe.calc.tpr -rerun ${eqpdb} >> gmx.log
              PE_md_equil=`grep -A 3 'Coulomb-14' spe.calc.log | tail -n 1 | awk '{printf "%d",$1}'`
	      
	      #Run production
	      if [ $rest == "true" ]; then
		  $python -u openmm_run_parmed_plumed.py -i prod_rest.inp -p ${optop} -c ${opgro} -irst ${eqrst} -orst ${prodrst} -odcd ${proddcd} > solution.${run}.prod.${curr}.out
	      else
		  $python -u openmm_run_parmed_plumed.py -i prod.inp -p ${optop} -c ${opgro} -irst ${eqrst} -orst ${prodrst} -odcd ${proddcd} > solution.${run}.prod.${curr}.out
	      fi
	      mv cv_ZZZ.dat cv_ZZZ_${curr}.dat
	      mv center.dat center_${curr}.dat
              echo "TITLE     <job>_<run>_${curr}" > temp.pdb
              echo "REMARK    PDB after production" >> temp.pdb
              grep "CRYST1" output.pdb >> temp.pdb
              echo "MODEL        1" >> temp.pdb
              #Using MMTSB PDB manipulation script
              ../bin/toolset/perl/convpdb.pl output.pdb | sed -e "s/HETATM/ATOM  /g" -e "s/END/ENDMDL/g" >> temp.pdb
              mv temp.pdb ${prodpdb}
              rm -rf output.pdb
	      
	      sed -e "s/<curr>/${curr}/g" nohyd.tmpl.py > nohyd.py
	      $python nohyd.py > nohyd.out
	      if [ -e solution.${run}.prod.${curr}.nohyd.dcd ]; then
		  rm -rf solution.${run}.prod.${curr}.dcd temp.nohyd.dcd
	      fi

              ${GMXDIR}/gmx grompp -f spe_calc.mdp -o spe.calc.tpr -c ${prodpdb} -p ${optop} >> gmx.log
              $mdrun -nt $NSLOTS -deffnm spe.calc -s spe.calc.tpr -rerun ${prodpdb} >> gmx.log
              PE_md_prod=`grep -A 3 'Coulomb-14' spe.calc.log | tail -n 1 | awk '{printf "%d",$1}'`

              status=`grep "%" solution.${run}.prod.${curr}.out | tail -n 1 | awk '{printf "%d",$1}'`

	      echo "${curr}-${cnt} ${PE_swap0} ${PE_equi0} ${PE_swap} ${PE_equi} ${PE_md_equil} ${PE_md_prod}" | awk '{print $1,"\t",$2,"\t",$3,"\t",$4,"\t",$5,"\t",$6,"\t",$7}' >> energy.dat
	   elif [ ${ifrags[1]} == 1 ] && [ ${curr} != 0 ]; then
	      PE_swap0=${PE_swap}; PE_equi0=${PE_equi}
	      cp -rf solution.${run}.gc.${curr}.pdb solution.${run}.prod.${curr}.pdb
              status=100
	   elif [ ${curr} == 0 ]; then
	       PE_swap0=${PE_swap}; PE_equi0=${PE_equi}
	       cp -rf solution.${run}.gc.${curr}.pdb solution.${run}.prod.${curr}.pdb
	       status=100
	   fi
	fi

	echo "${curr}-${cnt} ${ion_status} ${ifrags[0]} ${ifrags[1]} ${ifrags[2]} ${ifrags[3]} ${charge}" | awk '{print $1,"\t",$2,"\t",$3,"\t",$4,"\t",$5,"\t",$6,"\t",$7}' >> checkpoint.dat

        # Check to see if things worked. If not, re-attempt the GC step with different muex
        if [ $status == 100 ]; then
           # everything is OK
	   echo "checkpoint-${curr}-${cnt} ${status} ${ion_status} ${ifrags[1]} ${charge}"
   	   
           if [ ${ion_status} -eq 0 ]; then
               if [ ${ifrags[1]} -gt ${min_ions} ]; then
                  ion_status=0
               else
                  ion_status=1
                  x=0
                  while [ $x -le  $farraysize ];do
                        muexlist[$x]=${uhfe_list[$x]}
			maxfrags[$x]=${max_ion_comp[$x]}
                        minfrags[$x]=${max_ion_comp[$x]}
                        x=$((x + 1))
                  done
               fi
           elif [ ${ion_status} -eq 1 ]; then
              if [ ${ifrags[1]} -lt ${max_ions} ]; then
                  ion_status=1
              else
                  ion_status=0
        	  x=0
        	  while [ $x -le  $farraysize ];do
            	  	muexlist[$x]=${lhfe_list[$x]}
			maxfrags[$x]=${min_ion_comp[$x]}
			minfrags[$x]=${min_ion_comp[$x]}
                	x=$((x + 1))
        	  done
              fi
           fi
           prev=$((curr))
           curr=$((curr+1))
           cnt=1

           echo "declare -a muexlist=('${muexlist[0]}' '${muexlist[1]}' '${muexlist[2]}' '${muexlist[3]}');" > ${curr}.state
           echo "declare -a maxfrags=('${maxfrags[0]}' '${maxfrags[1]}' '${maxfrags[2]}' '${maxfrags[3]}');" >> ${curr}.state
           echo "declare -a minfrags=('${minfrags[0]}' '${minfrags[1]}' '${minfrags[2]}' '${minfrags[3]}');" >> ${curr}.state

	else
	   # this round failed, so back up (don't increment)
           echo "checkpoint-${curr}-${cnt} ${status} ${ion_status} ${ifrags[1]} ${charge}"
           curr=$((curr))
           prev=$((prev))
	   cnt=$((cnt+1))

           if [ ${ion_status} -eq 0 ];then
              muexlist[1]=`echo "${muexlist[1]} ${cnt}" | awk '{print $1-(50*($2-1))}'`
	      muexlist[2]=`echo "${muexlist[2]} ${cnt}" | awk '{print $1+(20*($2-1))}'`
           elif [ ${ion_status} -eq 1 ];then
              muexlist[1]=`echo "${muexlist[1]} ${cnt}" | awk '{print $1+(20*($2-1))}'`
	      muexlist[2]=`echo "${muexlist[2]} ${cnt}" | awk '{print $1-(10*($2-1))}'`
           fi

	fi

	if [ ${cnt} == 20 ];then
	   exit
	fi
done
rm -rf solution.*.equil.*.dcd solution.*.equil.*.rst

