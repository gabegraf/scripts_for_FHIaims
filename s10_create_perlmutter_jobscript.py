#!/usr/bin/env python
# Running this script will create a sample jobscript for SLURM 
# submission on Perlmutter in the working directory. Please make 
# sure all the SBATCH settings are as you want them.

import os

filename = "jobscript.sh"

# Get the current path and split into components
path_parts = os.path.abspath(os.getcwd()).split(os.sep)

# Safely get the last 3 directory names for the job name
if len(path_parts) >= 3:
    job_name = "_".join(path_parts[-3:])
else:
    job_name = "_".join(path_parts)

# Define the content of the job script
job_script_content = f"""#!/bin/bash
#SBATCH -J {job_name}
#SBATCH -A m3337
#SBATCH -C cpu
#SBATCH -q regular
#SBATCH -t 13:00:00
#SBATCH --nodes 2
#SBATCH --ntasks-per-node=128
#SBATCH -c 2
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=gabriel.graf@duke.edu

module load PrgEnv-intel
module load intel/2023.2.0
ulimit -s unlimited

# Load the working intel tuning file
export I_MPI_TUNING_BIN=/opt/intel/oneapi/mpi/2021.6.0/etc/tuning_generic_shm-ofi_mlx_hcoll.dat

srun -n ${{SLURM_NTASKS}} ~/FHIaims/build/aims.250312.scalapack.mpi.x > aims.out
"""

# Write the content to the file
with open(filename, "w") as file:
    file.write(job_script_content)

print(f"Job script written to {filename} with job name: {job_name}")

