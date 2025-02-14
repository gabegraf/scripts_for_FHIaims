#!/usr/bin/env python
# THis file creates a sample control.in file which can be tailored
# for specificic calculations. 
# 
# IMPORTANT: DO NOT RUN THE SAMPLE SCRIPT BLINDLY!
#
# Do help avoid this the k-grid is intentionally left blank, so
# that one must fill it in before using the control.in file.
#
# A same k-grid line is: k_grid 5 5 3 
#
# There are many keywords here and not all of them need to be used.
# Comment and un-comment them as necessary for your specific job.
filename = "control.in"

# Define the content of the job script
job_script_content = """  xc                 pbe
  spin               none
  relativistic       atomic_zora scalar
  vdw_correction_hirshfeld

#Large calculation
  use_local_index .true.

# k-grid settings
  k_grid

# relaxation
   # relax_geometry bfgs 5E-3
   # relax_unit_cell none
   # write_restart_geometry .true.

   # computes_forces .true.
"""

# Write the content to the file
with open(filename, "w") as file:
    file.write(job_script_content)

print(f"Job script written to {filename}")

