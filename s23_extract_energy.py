#!/usr/bin/env python
import os
import re
import sys

def extract_energy(filename):
    energy_pattern = r"\|\s+Total energy of the DFT / Hartree-Fock s\.c\.f\. calculation\s+:\s+(-?\d+\.\d+)\s+eV"

    if os.path.exists(filename):
        with open(filename, "r") as file:
            for line in file:
                match = re.search(energy_pattern, line)
                if match:
                    energy_value = float(match.group(1))
                    return energy_value
    else:
        print(f"Error: File '{filename}' does not exist.")
        return None

if __name__ == "__main__":
    # Use default filename unless an argument is provided
    filename = "aims.out"
    if len(sys.argv) > 1:
        filename = sys.argv[1]

    energy = extract_energy(filename)
    if energy is not None:
        print(energy)
