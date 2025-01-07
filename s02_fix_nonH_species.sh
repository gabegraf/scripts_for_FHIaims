# This script will selectively fix all non-Hydrogen atoms in your geometry.in file by appending "constrain_relaxation .true" in the following line

# A sample 



#to run this script,
########################################

#!/bin/bash

input_file="geometry.in"
output_file="geometry.tmp" #if you don't use a .tmp file it will delete all your coordinates when it clears the output file

# Clear the output file if it exists, or create a new one
> "$output_file"

# Read the input file line by line
while IFS= read -r line; do
    # Write the original line to the output file
    echo "$line" >> "$output_file"

    # Check if the line does NOT end with "H"
    if [[ ! "$line" =~ H$ ]] && [[ $line =~ ^atom ]]; then
        echo "constrain_relaxation .true." >> "$output_file"  # Append additional line
    fi
done < "$input_file"  # Reading input file here

# Move the temporary file to the original input file name
mv "$output_file" "$input_file"

echo "Processing complete. Results saved to $input_file"
