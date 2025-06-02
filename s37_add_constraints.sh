#!/bin/bash
# Script: Add or Remove Relaxation Constraints to Geometry File
# Description: This script processes a `geometry.in` file.
#              - By default, it adds `constrain_relaxation .true.` after atom lines
#                unless the element matches one of the arguments (e.g., H, C).
#              - With --remove or -r, it removes all `constrain_relaxation .true.` lines.

input_file="geometry.in"
output_file="geometry.tmp"

# Check if first argument is --remove or -r
if [[ "$1" == "--remove" || "$1" == "-r" ]]; then
    > "$output_file"
    while IFS= read -r line; do
        if [[ "$line" != "constrain_relaxation .true." ]]; then
            echo "$line" >> "$output_file"
        fi
    done < "$input_file"

    if cmp -s "$input_file" "$output_file"; then
        echo "No constraints found to remove."
    else
        mv "$output_file" "$input_file"
        echo "Processing complete. Constraints removed from $input_file"
    fi
else
    # List of elements to skip (all arguments)
    skip_elements=("$@")

    > "$output_file"
    while IFS= read -r line; do
        echo "$line" >> "$output_file"

        if [[ "$line" =~ ^atom ]]; then
            # Extract the element (last word)
            element=$(echo "$line" | awk '{print $NF}')

            # Check if element is in skip list
            skip=false
            for skip_element in "${skip_elements[@]}"; do
                if [[ "$element" == "$skip_element" ]]; then
                    skip=true
                    break
                fi
            done

            if [[ "$skip" == false ]]; then
                echo "constrain_relaxation .true." >> "$output_file"
            fi
        fi
    done < "$input_file"

    mv "$output_file" "$input_file"
    echo "Processing complete. Results saved to $input_file"
fi
