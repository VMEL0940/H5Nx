#!/bin/bash

# Set the Schrodinger installation path if it's not in your PATH
MERGE_COMMAND="${SCHRODINGER}/utilities/glide_merge"
ROOT_DIR=./

# Reference files for merging
ligand="RBS_sln.maegz"

# Loop through each rot-*.maegz file in the directory
for file in ${ROOT_DIR}*_name.maegz; do
    # Extract the base filename without the .maegz extension
    base_name=$(basename "${file%.maegz}")
    
    echo "Merging ${file} with $ligand"

    # Merge each rot file with ligand_26 and ligand_23, saving each result
    $MERGE_COMMAND -epv -o "${base_name}_complex.maegz" "$file" "$ligand"
    
done

echo "Merging complete. Each merged file is saved with the suffix '_complex.maegz'."
