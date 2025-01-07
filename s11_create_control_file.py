#!/usr/bin/env python
filename = "control.in"

# Define the content of the job script
job_script_content = """  xc                 pbe
  spin               none
  relativistic       atomic_zora scalar
  vdw_correction_hirshfeld

#Large calculation
  use_local_index .true.

# k-grid settings
  k_grid  5 5 2

# relaxation
   relax_geometry bfgs 5E-3
   relax_unit_cell none
   write_restart_geometry .true.
"""

# Write the content to the file
with open(filename, "w") as file:
    file.write(job_script_content)

print(f"Job script written to {filename}")

