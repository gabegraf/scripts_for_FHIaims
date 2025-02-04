#!/bin/bash
# Script: Add Relaxation Constraints to Geometry File
# Description: This script processes a `geometry.in` file and adds a relaxation constraint (`constrain_relaxation .true.`) to lines that define atoms, 
#              except for those ending with "H" (hydrogen).
#        The script reads from `geometry.in` and writes the modified content back to the same file.

# Example: If `geometry.in` contains:
#          atom 0.0 0.0 0.0 O
#          atom 1.0 1.0 1.0 H
#          The script will append `constrain_relaxation .true.` to the oxygen line but not the hydrogen line.

# Dependencies:
# - The script assumes the input file (`geometry.in`) exists and is formatted correctly.
# - The script requires write permissions for `geometry.in` and the directory it resides in.

# Exit Codes:
# - None explicitly defined, but errors may occur if the input file is missing or permissions are insufficient.
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
