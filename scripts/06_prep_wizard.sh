#!/bin/bash

ROOT_DIR=./

# Loop through all .maegz files in the current directory
for file in ${ROOT_DIR}*_complex-out.maegz; do
    echo "Processing ${file} with prepwizard"
    
    # Construct the output file name
    output_file="${file%.maegz}_pv.maegz"
    
    # Run prepwizard to generate the Pose Viewer file
    ${SCHRODINGER}/utilities/prepwizard  -nopreprocess -noprotassign -noimpref "${file}" "${output_file}"
    
done
