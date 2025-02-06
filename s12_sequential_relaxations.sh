#!/usr/bin/env bash
# This file is for setting up a series of relaxations for a specific material.
# The idea is to first relax with light settings, then intermediate, and then tight.
#
# To run this script you need to pass the species symbols as arguments, ie.
#
# s12_sequential_relaxations C H N Br Pb I
# 
# Additionally, you must have a jobscript.sh, a control.in, and a geometry.in in the
# working directory. The script will then make three sub-directories, called "01_light", 
# "02_intermediate", and "03_tight". It will move your files into the "01_light" directory
# and automatically append species defaults to your control.in. You must verify everything
# and then submit the job. After the job finishes, it should move your geometry.in.next_step
# and hessian.aims into the 02_intermediate folder automatically. The intermediate species
# defaults should already be added to your control.in here. This same process will happen
# after your job in the "02_intermediate" directory finishes.
#
# Currently, the script should not move your files into the next directory if the calculation
# failed. However, I would recommend checking this every time you do this.

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


