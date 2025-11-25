#!/usr/bin/env python
# Running this script will create a sample jobscript for SLURM 
# submission on Timewarp in the working directory. Please make 
# sure all the SBATCH settings are as you want them.

import os

filename = "jobscript.sh"

# Get the current path and split into components
path_parts = os.path.abspath(os.getcwd()).split(os.sep)

# Safely get the last 3 directory names
if len(path_parts) >= 3:
    job_name = "_".join(path_parts[-3:])
else:
    job_name = "_".join(path_parts)

# Define the content of the job script with dynamic job name
job_script_content = f"""#!/bin/bash
#SBATCH -J {job_name}
#SBATCH --ntasks-per-node=52
#SBATCH --nodes=2
#SBATCH --time=3:00:00
#SBATCH --account=2d1dpero
#SBATCH --partition=standard

ulimit -s unlimited
module load intel-oneapi
module load intel-oneapi-mkl
module load intel-oneapi-mpi
cd $SLURM_SUBMIT_DIR

srun -n 104 ~/software/FHIaims/build/aims.x > aims.out 2> aims.err
"""

# Write the content to the file
with open(filename, "w") as file:
    file.write(job_script_content)

print(f"Job script written to {filename} with job name: {job_name}")


