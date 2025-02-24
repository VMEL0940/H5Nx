#!/bin/bash

ROOT_DIR=./

# Step 1: Merge structures using pv_convert.py
echo "Merging .maegz files using pv_convert.py..."
for file in ${ROOT_DIR}*pv.maegz; do
    echo "Processing ${file} with pv_convert.py"
    
    # Construct the output file name
    output_file="${file%.maegz}_complex.maegz"
    
    # Run pv_convert.py in merge mode
    ${SCHRODINGER}/run pv_convert.py -mode merge -o "${output_file}" "${file}"
    
    if [ $? -eq 0 ]; then
        echo "Merged ${file} into ${output_file}"
    else
        echo "Error merging ${file}" >&2
        continue
    fi
done

# Step 2: Convert merged .maegz files to .pdb
echo "Converting merged .maegz files to .pdb format..."
#!/bin/bash

# Loop through each .maegz file in the current directory
for file in *_complex.maegz; do
    # Extract the base name without the .maegz extension
    base_name=$(basename "$file" .maegz)
    
    # Convert .maegz to .pdb
    ${SCHRODINGER}/utilities/structconvert -imae "$file" -opdb "${base_name}.pdb"

    # Check if the conversion was successful
    if [ $? -eq 0 ]; then
        echo "Successfully converted $file to ${base_name}.pdb"
    else
        echo "Error converting $file" >&2
    fi
done

echo "Conversion complete."
