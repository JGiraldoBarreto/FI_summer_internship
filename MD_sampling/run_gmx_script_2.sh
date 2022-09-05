#!/bin/bash

## SYNTAXIS:
## gmx_sample_script.sh "name_of_the_node" "name_of_the_file_with_rmsd"
## e.g: bash gmx_sample_script.sh node-2 node_2_rmsd box_size_node_2

for n in {2..9}; do bash gmx_sample_script_2.sh node-${n} node_${n}_rmsd ; done

echo "END PROGRAM"
