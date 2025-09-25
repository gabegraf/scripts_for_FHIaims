#!/usr/bin/env python
# ASE
from ase.io import read

# Read input geometry, save the lattice constant and convert to phonopy atoms object
initial_atoms = read("geometry.in", format="aims")

scale_factors = [0.96, 0.98, 1.0, 1.02, 1.04]

print(f"Initial volume:           {initial_atoms.get_volume():.3f} \u212B^3")
for ii, factor in enumerate(scale_factors):
    unitcell = initial_atoms.copy()
    unitcell.set_cell(initial_atoms.cell * factor, scale_atoms=True)
    volume = unitcell.get_volume()
    unitcell.write(f"geometry.in-{ii:03d}", format="aims", scaled=True)

    # prompt what we've done:
    print(f"Volume of sample {ii:4d}:    {unitcell.get_volume():.3f} \u212B^3")
