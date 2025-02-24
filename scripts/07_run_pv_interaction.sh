#!/bin/bash

# Set the Schrodinger installation path if it's not in your PATH
SCHRODINGER_PATH=/opt/schrodinger/suites2024-4
ROOT_DIR=./

# Loop through all _pv.maegz files in the current directory
for file in "${ROOT_DIR}"*_pv.maegz; do
    echo "Processing ${file} with poseviewer_interactions.py"

    
    # Run poseviewer_interactions.py on the input file
    "${SCHRODINGER_PATH}/run" poseviewer_interactions.py "${file}" -hbond_max_dist 4.1 -hbond_min_d_angle 100 -contacts hbond,salt -csv
    
done
