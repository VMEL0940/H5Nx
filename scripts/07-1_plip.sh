#!/bin/bash

# Loop through each PDB file in the current directory
for file in *.pdb; do
    # Extract the base name without the extension
    base_name=$(basename "$file" .pdb)
    
    # Run PLIP for each PDB file with the specified options
    plip -f "$file" -o "${base_name}_output" --rawstring --breakcomposite -tx

    # Check if PLIP ran successfully
    if [ $? -eq 0 ]; then
        echo "Successfully processed $file with PLIP"
    else
        echo "Error processing $file with PLIP" >&2
    fi
done

echo "PLIP analysis complete."
