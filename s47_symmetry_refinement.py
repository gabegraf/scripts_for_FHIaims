#!/usr/bin/env python
from pymatgen.core.structure import Structure
from pymatgen.io.aims.inputs import AimsGeometryIn
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.analysis.structure_matcher import StructureMatcher
import numpy as np
import click

@click.command()
@click.option("--tol", default=0.01, help="Tolerance for initial space group analysis")
@click.option("--input", default="geometry.in", help="Input geometry file")
@click.option("--primitive", is_flag=True, help="Convert to primitive cell")
@click.option("--dry", is_flag=True, help="Just return the symmetry but don't refine the structure")
@click.option("--conventional", is_flag=True, help="Convert to conventional standard cell")
@click.option("--rebuild", is_flag=True, help="Rebuild a cell like the input cell after refinement by translations, rotations, and supercells")
def symmetrization(input, tol, primitive, dry, conventional, rebuild):
    structure = AimsGeometryIn.from_file(input).structure
    sga = SpacegroupAnalyzer(structure, tol)
    print(
        f"Space Group Found! It's {sga.get_space_group_number()} ({sga.get_space_group_symbol()})")
    if dry:
        return  
    structure_refined = sga.get_refined_structure()
    if primitive and not conventional and not rebuild:
        structure_final = sga.get_primitive_standard_structure()
    elif conventional and not primitive and not rebuild:
        structure_final = sga.get_conventional_standard_structure()
    elif rebuild and not conventional and not primitive:
        matcher = StructureMatcher(
            primitive_cell=False,
            scale=False,
            attempt_supercell=True)
        structure_final = matcher.get_s2_like_s1(structure, structure_refined)
        # supercell_transform = matcher.get_supercell_matrix(structure, structure_final)
        # supercell_to_primitive= np.linalg.inv(supercell_transform)
    else:
        raise ValueError("Only one of --primitive, --conventional, or --rebuild can be specified at a time.")


    AimsGeometryIn.from_structure(structure_final).write_file(overwrite=True)

    print(f"Symmetrized geometry written to geometry.in")
    write_metadata(sga, tol, primitive, conventional, rebuild)

def write_metadata(sga_og, tol, primitive, conventional, rebuild):
    with open("symmetry_refinement_metadata.txt", "w") as f:
        f.write("\n# Symmetry Refinement Metadata\n")
        f.write(f"# Space Group Found: {sga_og.get_space_group_number()} ({sga_og.get_space_group_symbol()})\n")
        f.write(f"# Tolerance: {tol}\n")
        f.write(f'# Converted to Primitive Cell: {"Yes" if primitive else "No"}\n')
        f.write(f'# Converted to Conventional Cell: {"Yes" if conventional else "No"}\n')
        f.write(f'# Rebuilt Structure: {"Yes" if rebuild else "No"}\n')

if __name__ == "__main__":
    symmetrization()