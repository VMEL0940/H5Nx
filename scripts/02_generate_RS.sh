#!/bin/bash

# Output file for the final .sh script
output_script="run_residue_scanning.sh"

# Header for the .sh script
echo "#!/bin/bash" > $output_script
echo "" >> $output_script

# Initialize a counter
counter=0
reference="4K63_apo.maegz"

# Loop over all .txt files in the folder, excluding the specific one
for txt_file in *.txt; do
    jobname="${txt_file%.txt}"

    # Write the SCHRODINGER command with the replaced mutations file to the output script
    echo -n "\"\${SCHRODINGER}/run\" residue_scanning_backend.py -fast -jobname ${jobname} -muts_file ${txt_file} -refine_mut prime_residue -calc hydropathy,rotatable,vdw_surf_comp,sasa_polar,sasa_nonpolar,sasa_total -dist 0.00 ${reference} -add_res_scan_wam -HOST localhost:1 -NJOBS 1 -TMPLAUNCHDIR" >> $output_script

    # Increment the counter
    counter=$((counter + 1))

    # Add -WAIT as a new line after every nth command
    if (( counter % 5 == 0 )); then
        echo " -WAIT" >> $output_script
    else
        echo "" >> $output_script
    fi
done

# Make the output script executable
chmod +x $output_script

echo "Script run_residue_scanning.sh has been created."
