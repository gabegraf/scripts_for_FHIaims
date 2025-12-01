#!/usr/bin/env python
from ase.io import read
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
import numpy as np
import click


@click.command()
@click.option("--tol", default=0.01, help="Tolerance for initial space group analysis")
@click.option("--output", default="geometry_refined.in", help="Name of the output file")
@click.option("--primitive", is_flag=True, help="Convert to primitive cell")
@click.option("--niggli", is_flag=True, help="Reduce the lattice to a standardized Niggli cell")

def symmetrization(output, tol, primitive, niggli):
    atoms = read("geometry.in", format="aims")
    structure = AseAtomsAdaptor.get_structure(atoms)
    sga_og = SpacegroupAnalyzer(structure, tol)
    print(
        f"Space Group Found! It's {sga_og.get_space_group_number()} ({sga_og.get_space_group_symbol()})"
    )

    structure = sga_og.get_refined_structure()
    sga = SpacegroupAnalyzer(structure, tol)
    if primitive:
        structure = sga.get_primitive_standard_structure()
        sga = SpacegroupAnalyzer(structure, tol)
    if niggli:
        structure = structure.get_reduced_structure()
        sga = SpacegroupAnalyzer(structure, tol)

    atoms = AseAtomsAdaptor.get_atoms(structure)
    atoms.center()
    atoms.wrap()
    atoms.write(output, format="aims")
    print(f"Symmetrized geometry written to {output}")
    write_metadata(sga_og, output, tol, primitive, niggli)


def write_metadata(sga_og, output, tol, primitive, niggli):
    with open("symmetry_refinement_metadata.txt", "a") as f:
        f.write("\n# Symmetry Refinement Metadata\n")
        f.write(
            f"# Space Group Found: {sga_og.get_space_group_number()} ({sga_og.get_space_group_symbol()})\n"
        )
        f.write(f"# Tolerance: {tol}\n")
        f.write(f'# Converted to Primitive Cell: {"Yes" if primitive else "No"}\n')
        f.write(f'# Converted to Niggli Reduced Lattice: {"Yes" if niggli else "No"}\n')


if __name__ == "__main__":
    symmetrization()
