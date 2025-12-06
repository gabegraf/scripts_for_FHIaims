#!/usr/bin/env python
from pymatgen.core.structure import Structure
from pymatgen.io.aims.inputs import AimsGeometryIn
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.io.aims.inputs import AimsGeometryIn
import click
import numpy as np

@click.command()
@click.option("--primitive_file", default="geometry_primitive.in", help="Path to the primitive cell geometry file")
@click.option("--supercell_file", default="geometry_supercell.in", help="Path to the supercell geometry file")
def main(primitive_file, supercell_file):
    prim = AimsGeometryIn.from_file(primitive_file).structure
    sup = AimsGeometryIn.from_file(supercell_file).structure
    supercell_to_primitive = find_supercell_transform(prim, sup)
    print(supercell_to_primitive.flatten())


def find_supercell_transform(primitive_structure, supercell_structure):

    matcher = StructureMatcher(
        primitive_cell=False,
        scale=False,
        attempt_supercell=True)

    supercell_transform = matcher.get_supercell_matrix(supercell_structure, primitive_structure)
    supercell_to_primitive= np.linalg.inv(supercell_transform)
    return supercell_to_primitive

if __name__ == "__main__":
    main()