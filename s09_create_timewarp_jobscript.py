#!/usr/bin/env python
# Running this script will create a sample jobscript for SLURM 
# submission on Timewarp in the working directory. Please make 
# sure all the SBATCH settings are as you want them.

filename = "jobscript.sh"

# Define the content of the job script
job_script_content = """#!/bin/bash
#SBATCH --job-name=Forgot_to_write_a_name
#SBATCH -N 6
#SBATCH --ntasks-per-node=20
#SBATCH --cpus-per-task=2
#SBATCH -p long
#SBATCH -t 30:00:00

# initialization
module load FHIaims

# For GPU (node18)
# module load cuda-12.3

# execution
cd $SLURM_SUBMIT_DIR/
srun -n $SLURM_NTASKS aims.x > aims.out 2>aims.err

# OR
# srun -n $SLURM_NTASKS aims.x > aims.out 2>aims.err
"""

# Write the content to the file
with open(filename, "w") as file:
    file.write(job_script_content)

print(f"Job script written to {filename}")

