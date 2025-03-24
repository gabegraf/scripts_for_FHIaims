#!/usr/bin/env python
# Running this script will create a sample jobscript for SLURM 
# submission on Timewarp in the working directory. Please make 
# sure all the SBATCH settings are as you want them.

filename = "jobscript.sh"

# Define the content of the job script
job_script_content = """
#!/bin/bash
#SBATCH -p scavenger
#SBATCH --nodelist=dcc-courses-[1-50]
#SBATCH -N 2
#SBATCH --ntasks-per-node=42
#SBATCH -c 2
#SBATCH --mem-per-cpu=1700
#SBATCH --mail-type=end
#SBATCH --mail-user=gabriel.graf@duke.edu

# Initialization
source ~/.bashrc
module load FHIaims-intel
ulimit -s unlimited

# Execution
cd $SLURM_SUBMIT_DIR
mpirun -n $SLURM_NTASKS aims.x > aims.out 2> aims.err

### If you built from source, use the following instead.
# module load compiler/latest
# module load mkl/latest
# module load mpi/latest
# ulimit -s unlimited

# cd $SLURM_SUBMIT_DIR
# mpirun -n $SLURM_NTASKS /hpc/home/ukh/apps/FHIaims/build/aims.X.scalapack.mpi.x > aims.out 2> aims.err
"""

# Write the content to the file
with open(filename, "w") as file:
    file.write(job_script_content)

print(f"Job script written to {filename}")

