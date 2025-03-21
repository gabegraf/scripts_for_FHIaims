#!/bin/bash

# Check if the --be_nice flag is passed
be_nice=false
for arg in "$@"; do
    if [[ "$arg" == "--be_nice" ]]; then
        be_nice=true
        break
    fi
done

for file in geometry.in-*; do
    number=$(echo "$file" | grep -oP '(?<=geometry.in-)\d+')

    folder="${number}-displacement"

    mkdir -p "$folder"
    mv "$file" "$folder/geometry.in"

    cp jobscript.sh control.in "$folder/"

    # If --be_nice flag is passed, run the job with nice 100
    if $be_nice; then
        cd "$folder" && sbatch --nice=100 jobscript.sh && cd ../
    else
        cd "$folder" && sbatch jobscript.sh && cd ../
    fi
done

echo "All jobs have been submitted!"

