#!/usr/bin/env python3

from ase.io import read, write

atoms = read("geometry.in", format="aims")

lower_cell, rotation = atoms.cell.standard_form("lower")
atoms.set_cell(lower_cell, scale_atoms=True)
atoms.wrap()

write("geometry2.in", atoms, format="aims")
