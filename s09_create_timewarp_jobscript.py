#!/usr/bin/env python
filename = "jobscript.sh"

# Define the content of the job script
job_script_content = """#!/bin/bash
#SBATCH --job-name=r2MFSA_PbI_12_298K_red_h_relaxation
#SBATCH -N 4
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=2
#SBATCH -p small
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

