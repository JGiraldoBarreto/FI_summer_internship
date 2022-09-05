#!/bin/bash

## Sampling program using gmx with plumed.
## $1 ---> name of node of the pat, without the .pdb extension
## $2 ---> file with the rmsd open - closed restraint values.

## To run the program: bash run_gmx_script.sh

rm $1_folder/$2 $1_folder/plumed.dat $1_folder/$1_trajout.pdb
rm -r $1_folder/traj_frames/

cp New_restraint_files/$2 $1_folder/.

cd $1_folder/

var1=`awk '{print $1}' $2` ; var2=`awk '{print $2}' $2`

echo -e "
UNITS LENGTH=A

rmsd1: RMSD REFERENCE=../open_state.pdb TYPE=OPTIMAL
rmsd2: RMSD REFERENCE=../closed_state.pdb TYPE=OPTIMAL
restraint: RESTRAINT ARG=rmsd1,rmsd2 AT=$var1,$var2 KAPPA=50,50

PRINT STRIDE=25 ARG=rmsd1,rmsd2 FILE=COLVAR" > plumed.dat

#echo 6 | gmx pdb2gmx -f $1.pdb -o $1_processed.gro -water tip3p
#gmx editconf -f $1_processed.gro -o $1_newbox.gro -c -d $var3 -bt cubic
#gmx solvate -cp $1_newbox.gro -cs spc216.gro -o $1_solv.gro -p topol.top
#gmx grompp -f ions.mdp -c $1_solv.gro -p topol.top -o ions.tpr
#gmx genion -s ions.tpr -o $1_solv_ions.gro -p topol.top -pname NA -nname CL -neutral
#gmx grompp -f minim.mdp -c $1_solv_ions.gro -p topol.top -o em.tpr
#gmx mdrun -v -deffnm em
#gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
#gmx mdrun -deffnm nvt
#gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr
#gmx mdrun -deffnm npt
gmx mdrun -s npt.tpr -plumed plumed.dat -nb cpu -v

echo 14 | gmx trjconv -f traj.trr -s npt.tpr -o $1_trajout.pdb

bash frames_trajectories_extract.sh $1_trajout.pdb 

cd ../

echo "Done"
