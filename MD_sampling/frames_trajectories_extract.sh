#!/bin/bash

#Separar los frames de las trayectorias de MD para obtener los FEP.
#$1=trajectory file

m=1
mkdir -p traj_frames

### Bloque 1: En este bloque se especifica la cantidad de fotogramas en la trayectoria###

awk 'BEGIN{n=1}{if($1=="ENDMDL")n++; print >> "traj_frames/traj_"n".pdb"}' $1

sed -i '1d' traj_frames/traj_*   #Borra las filas con la palabra END.
find . -type f -empty -delete     #Borra los archivos sin info, para que no entren en el proceso.

while read line
do

        if [ "$line" == "ENDMDL" ]; then
                let "m=m+1"
                echo "traj_frames/traj_$m" >> traj_names
        fi

done<$1

for n in {1..100}; do rm traj_frames/traj_${n}.pdb; done

echo "$m"

###END Program###
