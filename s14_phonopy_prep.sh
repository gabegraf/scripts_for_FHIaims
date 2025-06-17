#!/bin/bash

# Check if the --be_nice flag is passed
be_nice=false
no_sub=false
for arg in "$@"; do
    if [[ "$arg" == "--be_nice" ]]; then
        be_nice=true
    elif [[ "$arg" == "--no-sub" ]]; then
        no_sub=true
    fi
done

for file in geometry.in-*; do
    number=$(echo "$file" | grep -oP '(?<=geometry.in-)\d+')
    # Pad with zeros to 4 digits (e.g., 1 → 0001, 12 → 0012, 123 → 0123, 1234 → 1234)
    printf -v number "%04d" "$number"

    folder="${number}-displacement"

    mkdir -p "$folder"
    mv "$file" "$folder/geometry.in"

    cp jobscript.sh control.in "$folder/"

    # If --no_sub flag is passed, skip submission
    if $no_sub; then
        echo "Skipping submission for $folder (--no-sub flag set)"
        continue
    fi

    # If --be_nice flag is passed, run the job with nice 100
    if $be_nice; then
        cd "$folder" && sbatch --nice=100 jobscript.sh && cd ../
    else
        cd "$folder" && sbatch jobscript.sh && cd ../
    fi
done

if $no_sub; then
    echo "All folders have been prepared but no jobs were submitted!"
else
    echo "All jobs have been submitted!"
fi
