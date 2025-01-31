#!/usr/bin/env bash

# Check if no arguments are provided
if [ $# -eq 0 ]; then
    echo "Error: No arguments provided."
    echo "Usage: Add the species for your system, they are need for appending species defaults to the control.in. Ex: 'C N H Pb I'."
    exit 1
fi

# Exit immediately if any command fails
set -e

# Step 1: Create directories
mkdir -p ./01_light ./02_intermediate ./03_tight

# Step 2: Copy necessary files to each directory
cp ./control.in ./geometry.in ./jobscript.sh ./01_light/
cp ./control.in ./jobscript.sh ./02_intermediate/
cp ./control.in ./jobscript.sh ./03_tight/

# Step 3: Run `s01_add_species_defaults.sh` in each directory
cd 01_light && s01_add_species_defaults.sh light "$@" && cd ..
cd 02_intermediate && s01_add_species_defaults.sh intermediate "$@" && cd ..
cd 03_tight && s01_add_species_defaults.sh tight "$@" && cd ..

# Step 4: Append commands to `jobscript.sh` in `01_light` and `02_intermediate`
sed -i '/^srun/ s|$| \&\& if [[ $(tail -n 2 aims.out | head -n 1) == *"Have a nice day"* ]]; then \
cp geometry.in.next_step hessian.aims ../02_intermediate/ \&\& mv ../02_intermediate/geometry.in.next_step \
../02_intermediate/geometry.in \&\& cd ../02_intermediate \&\& sbatch jobscript.sh \&\& cd ..; fi|' ./01_light/jobscript.sh

sed -i '/^srun/ s|$| \&\& cp geometry.in.next_step hessian.out ../03_tight/ \&\& mv ../03_tight/geometry.in.next_step ../03_tight/geometry.in|' ./02_intermediate/jobscript.sh


