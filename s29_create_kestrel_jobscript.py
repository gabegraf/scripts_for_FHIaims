#!/usr/bin/env python
# Running this script will create a sample jobscript for SLURM 
# submission on Timewarp in the working directory. Please make 
# sure all the SBATCH settings are as you want them.

filename = "jobscript.sh"

# Define the content of the job script
job_script_content = """#!/bin/bash
#SBATCH -J forgot_to_name 
#SBATCH --ntasks-per-node=52
#SBATCH --nodes=2
#SBATCH --time=3:00:00
#SBATCH --account=hybridpero
#SBATCH --partition=standard
#SBATCH --mail-type=fail
#SBATCH --mail-type=end
#SBATCH --mail-user=gabriel.graf@duke.edu

ulimit -s unlimited
module load intel-oneapi
module load intel-oneapi-mkl
module load intel-oneapi-mpi
cd $SLURM_SUBMIT_DIR

srun -n 104 ~/software/FHIaims/build/aims.250320.scalapack.mpi.x > aims.out 2> aims.err
"""

# Write the content to the file
with open(filename, "w") as file:
    file.write(job_script_content)

print(f"Job script written to {filename}")

