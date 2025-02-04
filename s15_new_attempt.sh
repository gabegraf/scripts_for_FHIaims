#!/bin/bash
# This script is for documenting all your failures!
# When you run a calculation in a folder and it doesn't work, you can run this script
# to automatically put all those files in a folder called 01_attempt and copies the 
# geometry.in jobscript.sh control.in and any files you provide as arguments into a 
# new folder called 02_attempt. From there you can make the changes and resubmit.

# If you mess up in i_attempt (where i = 02, 03, 04, ...), running this script inside that attempt
# directory will make an i+1 attempt directory and copy the files above from the current
# attempt directory into the i+1 attempt directory. 

# List of mandatory files to copy
mandatory_files=("control.in" "geometry.in" "jobscript.sh")

# Check if *_attempt folder exists in the parent directory
highest_attempt=$(ls -d ../*_attempt 2>/dev/null | sort -V | tail -n 1)

# If no *_attempt folder exists, create 01_attempt in the current directory and move files there
if [ -z "$highest_attempt" ]; then
    new_dir="01_attempt"
    mkdir "$new_dir"
    echo "Created $new_dir in the current directory"
    
    # Move all files from the current directory to 01_attempt
    for file in *; do
        if [ -f "$file" ] && [ "$file" != "$0" ]; then
            mv "$file" "$new_dir"
            echo "Moved $file to $new_dir"
        fi
    done

    # After moving files to 01_attempt, create 02_attempt in the current directory
    next_dir="02_attempt"
    mkdir "$next_dir"
    echo "Created $next_dir in the current directory"
    
    # Copy mandatory files and any argument files from 01_attempt to 02_attempt
    for file in "${mandatory_files[@]}"; do
        if [ -f "$new_dir/$file" ]; then
            cp "$new_dir/$file" "$next_dir"
            echo "Copied $file from $new_dir to $next_dir"
        else
            echo "$file not found in $new_dir, skipping."
        fi
    done

    # Copy argument files (non-mandatory) to 02_attempt
    for file in "$@"; do
        if [ -f "$new_dir/$file" ]; then
            cp "$new_dir/$file" "$next_dir"
            echo "Copied $file from $new_dir to $next_dir"
        else
            echo "$file not found in $new_dir, skipping."
        fi
    done

else
    # Go up one directory
    cd ..

    # Extract the numeric part of the highest attempt folder and increment by 1 with leading zero
    base_num=$(echo "$highest_attempt" | grep -oE '[0-9]+' | tail -n 1)
    next_num=$((base_num + 1))
    new_dir=$(printf "%02d_attempt" $next_num)  # Use leading zero for single-digit numbers
    mkdir "$new_dir"
    echo "Created $new_dir in the parent directory"

    # Now copy mandatory files and argument files from the previous attempt folder
    prev_dir=$(printf "%02d_attempt" $base_num)
    for file in "${mandatory_files[@]}"; do
        if [ -f "$prev_dir/$file" ]; then
            cp "$prev_dir/$file" "$new_dir"
            echo "Copied $file from $prev_dir to $new_dir"
        else
            echo "$file not found in $prev_dir, skipping."
        fi
    done

    # Copy argument files (non-mandatory) from the previous attempt folder
    for file in "$@"; do
        if [ -f "$prev_dir/$file" ]; then
            cp "$prev_dir/$file" "$new_dir"
            echo "Copied $file from $prev_dir to $new_dir"
        else
            echo "$file not found in $prev_dir, skipping."
        fi
    done
fi

# Ensure the new directory is populated as requested
echo "Files copied to $new_dir"
