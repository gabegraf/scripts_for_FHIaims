#!/bin/bash

# Find all directories ending with '-displacement' and sort them
for dir in $(find . -type d | sort); do

    # Define the path to aims.out
    aims_out="$dir/aims.out"

    # Check if aims.out exists in the directory
    if [[ -f "$aims_out" ]]; then

        # Check if "Have a nice day" exists in the file
        if grep -q "Have a nice day" "$aims_out"; then
            echo "Success: 'Have a nice day' found in $aims_out"
        else
            echo "$dir"  # Print directory name if phrase is not found
        fi
    else
        echo "Error: $aims_out does not exist"
    fi
done

