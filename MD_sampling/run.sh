#!/bin/bash

rm rmsd_md_path.txt

cp ini_rmsd_md_path.txt rmsd_md_path.txt

bash run_gmx_script.sh

for i in {1..3}; do bash run_python_script.sh; bash run_gmx_script_2.sh ; done

bash run_python_script.sh

echo "End sampling optimization"
