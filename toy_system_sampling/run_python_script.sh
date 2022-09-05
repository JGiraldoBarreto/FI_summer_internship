#!/bin/bash

## SYNTAX:
## python3 cryoCV.py kappa_value images_file path_file
## e.g python3 cryoCV.py 2.5 Images_from_Grid_for_Hexapeptide.txt rmsd_md_path.txt

#for temp in {1.0,1.5,2.0,2.5,3.5}; do python3 cryoCV_diff_kappa.py 0.0 Images_from_Grid_for_Hexapeptide_${temp}.txt rmsd_md_path.txt ${temp} ; done

temp=1.0; python3 cryoCV.py 0.0 Images_from_Grid_for_Hexapeptide_${temp}.txt rmsd_md_path.txt ${temp}
