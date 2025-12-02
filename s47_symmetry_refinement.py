#!/usr/bin/env python
from pymatgen.core.structure import Structure
from pymatgen.io.aims.inputs import AimsGeometryIn
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
import numpy as np
import click

@click.command()
@click.option("--tol", default=0.01, help="Tolerance for initial space group analysis")
@click.option("--output", default="geometry_refined.in", help="Name of the output file")
@click.option("--primitive", is_flag=True, help="Convert to primitive cell")

def symmetrization(output, tol, primitive):
    structure = AimsGeometryIn.from_file("geometry.in").structure
    sga = SpacegroupAnalyzer(structure, tol)
    print(
        f"Space Group Found! It's {sga.get_space_group_number()} ({sga.get_space_group_symbol()})")

    structure = sga.get_refined_structure()
    if primitive:
        structure = sga.get_primitive_standard_structure()
        sga = SpacegroupAnalyzer(structure, tol)

    AimsGeometryIn.from_structure(structure).write_file(overwrite=True)

    print(f"Symmetrized geometry written to {output}")
    write_metadata(sga, output, tol, primitive)

def write_metadata(sga_og, output, tol, primitive):
    with open("symmetry_refinement_metadata.txt", "w") as f:
        f.write("\n# Symmetry Refinement Metadata\n")
        f.write(f"# Space Group Found: {sga_og.get_space_group_number()} ({sga_og.get_space_group_symbol()})\n")
        f.write(f"# Tolerance: {tol}\n")
        f.write(f'# Converted to Primitive Cell: {"Yes" if primitive else "No"}\n')

if __name__ == "__main__":
    symmetrization()
