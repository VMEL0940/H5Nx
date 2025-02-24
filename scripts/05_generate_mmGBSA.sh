#!/bin/bash

# Define the base command
BASE_COMMAND="${SCHRODINGER}/prime_mmgbsa -OVERWRITE -prime_opt OPLS_VERSION=S-OPLS"

# Define the root directory (current directory)
ROOT_DIR="./"  # You can adjust this if needed, but it points to the current directory

# Define the output file to store the generated commands
OUTPUT_FILE="mmgbsa_commands.sh"

# Clear the output file if it exists
> $OUTPUT_FILE

# Initialize a counter
counter=0

# Loop through each .maegz file in the current directory only
for file in ${ROOT_DIR}*.maegz; do
    echo -n "${BASE_COMMAND} ${file} -HOST localhost:18 -NJOBS 1 -TMPLAUNCHDIR" >> $OUTPUT_FILE

    # Increment the counter
    counter=$((counter + 1))

    # Add -WAIT after every nth line
    if (( counter % 5 == 0 )); then
        echo " -WAIT" >> $OUTPUT_FILE
    else
        echo "" >> $OUTPUT_FILE
    fi
done

# Make the output script executable
chmod +x $OUTPUT_FILE

echo "MMGBSA commands have been written to $OUTPUT_FILE"
