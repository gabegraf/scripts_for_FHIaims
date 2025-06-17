#!/bin/bash

# Handle flags (--be_nice, --no-sub)
be_nice=false
no_sub=false
for arg in "$@"; do
    if [[ "$arg" == "--be_nice" ]]; then
        be_nice=true
    elif [[ "$arg" == "--no-sub" ]]; then
        no_sub=true
    fi
done

# Process files in *natural order* (1, 2, 10, not 1, 10, 2)
for file in $(ls geometry.in-* | sort -V); do
    # Extract the exact number (no padding, no octal conversion)
    number=$(echo "$file" | grep -oP '(?<=geometry.in-)\d+')

    # Create folder like "displacement-1", "displacement-42", etc.
    folder="displacement-${number}"
    mkdir -p "$folder"
    mv "$file" "$folder/geometry.in"
    cp jobscript.sh control.in "$folder/"

    # Skip submission if --no-sub is set
    if $no_sub; then
        continue
    fi

    # Submit job (with --be_nice if set)
    if $be_nice; then
        cd "$folder" && sbatch --nice=100 jobscript.sh && cd ../
    else
        cd "$folder" && sbatch jobscript.sh && cd ../
    fi
done

if $no_sub; then
    echo "All folders prepared (no jobs submitted)."
else
    echo "All jobs submitted!"
fi
