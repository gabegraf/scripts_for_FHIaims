#!/bin/bash

# Read the current z-coordinate from the first line of geometry.in (4th field)
current_z=$(awk 'NR==1 {print $4}' geometry.in)

# Calculate the maximum z value (current_z * 2) using awk
max_z=$(awk -v z="$current_z" 'BEGIN {print z * 2}')

# Calculate the step size for 19 intervals between 0.5 and max_z (20 total points) using awk
step=$(awk -v max="$max_z" 'BEGIN {print (max - 0.5) / 19}')

# Loop to create 20 files
for i in {01..20}; do
    # Calculate the new z-coordinate: 0.5 + (i-1)*step using awk
    new_z=$(awk -v step="$step" -v idx="$i" 'BEGIN {print 0.5 + (idx - 1) * step}')
    
    # Create the directory
    mkdir -p "$i-displacement"
    
    # Copy the original file and modify the first line
    awk -v z="$new_z" 'NR==1 {$4=sprintf("%.6f", z); print $1" "$2" "$3" "$4" "$5} NR>1 {print}' geometry.in > "$i-displacement/geometry.in"
    cp control.in jobscript.sh $i-displacement/
    cd $i-displacement
    sbatch jobscript.sh
    cd ..
done

echo "Created 20 modified geometry.in files in their respective folders."
