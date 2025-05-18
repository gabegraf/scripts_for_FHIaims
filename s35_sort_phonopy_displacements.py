#!/usr/bin/env python
import os
import shutil
import numpy as np
from ase.io.aims import read_aims
from glob import glob

def main(ref_file='geometry.in', pattern='geometry.in-*', round_to=4):
    # Load reference geometry
    ref_atoms = read_aims(ref_file)
    ref_positions = ref_atoms.get_positions()
    lattice = ref_atoms.cell.array.T  # 3x3 matrix (a, b, c as columns)

    # Normalize lattice vectors
    lattice_unit = lattice / np.linalg.norm(lattice, axis=0)

    files = sorted(glob(pattern))
    for file in files:
        disp_atoms = read_aims(file)
        disp_positions = disp_atoms.get_positions()

        displacement = disp_positions - ref_positions
        atom_index = np.argmax(np.linalg.norm(displacement, axis=1))
        disp_vec = displacement[atom_index]

        # Project displacement onto each lattice vector
        projections = [np.dot(disp_vec, lv) for lv in lattice_unit.T]
        abs_proj = np.abs(projections)
        i_coord = np.argmax(abs_proj)
        delta = projections[i_coord]
        delta_rounded = round(delta, round_to)

        # Format folder name without '+' and trailing zeros
        delta_str = f"{delta_rounded:.{round_to}f}".rstrip('0').rstrip('.')

        folder_name = f"aims.i_atom_{atom_index}.i_coord_{i_coord}.displ_{delta_str}"
        os.makedirs(folder_name, exist_ok=True)

        new_path = os.path.join(folder_name, "geometry.in")
        shutil.move(file, new_path)
        print(f"Moved {file} to {new_path}")

if __name__ == "__main__":
    main()
