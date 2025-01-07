#!/usr/bin/env python3

import os
import re

# Set the root directory containing the subdirectories
root_dir = os.getcwd()  # or specify the full path

# Regular expression to match the energy line in aims.out
energy_pattern = r"\|\s+Total energy of the DFT / Hartree-Fock s\.c\.f\. calculation\s+:\s+(-?\d+\.\d+)\s+eV"

# Dictionary to store the energy values for each subdirectory
energies = {}

# Loop over each subdirectory
for n in range(1, 17):
    sub_dir = f"{n}_structure_light_relax"
    aims_out_path = os.path.join(root_dir, sub_dir, "aims.out")
    
    if os.path.isfile(aims_out_path):
        with open(aims_out_path, "r") as file:
            for line in file:
                match = re.search(energy_pattern, line)
                if match:
                    energy_value = float(match.group(1))
                    energies[sub_dir] = energy_value
                    break

# Find the subdirectory with the lowest energy value
if energies:
    lowest_energy_subdir = min(energies, key=energies.get)
    lowest_energy_value = energies[lowest_energy_subdir]
    print(f"The subdirectory with the lowest energy is: {lowest_energy_subdir}")
    print(f"Lowest energy value: {lowest_energy_value} eV")
else:
    print("No energy values found.")
