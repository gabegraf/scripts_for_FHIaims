#!/usr/bin/env python
def extract_energy(filename)
    import os 
    import re

    energy_pattern = r"\|\s+Total energy of the DFT / Hartree-Fock s\.c\.f\. calculation\s+:\s+(-?\d+\.\d+)\s+eV"

        if os.path.exists(filename):
            with open(filename, "r") as file:
                for line in file:
                    match = re.search(energy_pattern, line)
                    if match:
                        energy_value = float(match.group(1))
                        break
    return energy_value


if __name__ == "__main__":
    filename = "aims.out"
    energy = extract_energy(filename)
    print(energy)
