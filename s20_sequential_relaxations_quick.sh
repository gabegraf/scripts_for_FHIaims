#!/usr/bin/env bash

# Check if no arguments are provided
if [ $# -eq 0 ]; then
    echo "Error: No arguments provided."
    echo "Usage: Add the species for your system, they are needed for appending species defaults to the control.in. Ex: 'C N H Pb I'."
    exit 1
fi

# Exit immediately if any command fails
set -e

# Step 1: Create directories
mkdir -p ./01_light ./02_intermediate ./03_tight

# Step 2: Copy necessary files to each directory
cp ./control.in ./geometry.in ./01_light/
cp ./control.in ./02_intermediate/
cp ./control.in ./03_tight/

# Step 3: Run `s01_add_species_defaults.sh` in each directory
cd 01_light && s01_add_species_defaults.sh light "$@" && cd ..
cd 02_intermediate && s01_add_species_defaults.sh intermediate "$@" && cd ..
cd 03_tight && s01_add_species_defaults.sh tight "$@" && cd ..

# Append execution block to jobscript.sh
cat <<EOF >> jobscript.sh

# Move to the first directory
cd 01_light

# Run FHI-aims simulation
srun -n \${SLURM_NTASKS} ~/FHIaims/build/aims.250131.scalapack.mpi.x > aims.out

# Check for successful completion
if tail -n 2 aims.out | head -n 1 | grep -q "Have a nice day"; then
    if [[ -f geometry.in.next_step ]]; then
        cp geometry.in.next_step hessian.aims ../02_intermediate/
        mv ../02_intermediate/geometry.in.next_step ../02_intermediate/geometry.in
    fi
    cd ../02_intermediate
    srun -n \${SLURM_NTASKS} ~/FHIaims/build/aims.250131.scalapack.mpi.x > aims.out
fi

# Check again after second simulation
if tail -n 2 aims.out | head -n 1 | grep -q "Have a nice day"; then
    if [[ -f geometry.in.next_step ]]; then
        cp geometry.in.next_step hessian.aims ../03_tight/
        mv ../03_tight/geometry.in.next_step ../03_tight/geometry.in
    fi
    cd ../03_tight
    srun -n \${SLURM_NTASKS} ~/FHIaims/build/aims.250131.scalapack.mpi.x > aims.out
fi
EOF

