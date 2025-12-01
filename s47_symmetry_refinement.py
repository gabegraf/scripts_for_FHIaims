#!/usr/bin/env python
from ase.io import read
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
import numpy as np
import click

@click.command()
@click.option("--tol", default=0.001, help="Tolerance for initial space group analysis")
@click.option("--output", default="geometry_refined.in", help="Name of the output file")
@click.option("--no_primitive", is_flag=True, help="Do not convert to primitive cell")

def symmetrization(output, tol, no_primitive):
  atoms = read('geometry.in', format='aims')
  structure = AseAtomsAdaptor.get_structure(atoms)
  sga = SpacegroupAnalyzer(structure, tol)
  print(f"Space Group Found! It's {sga.get_space_group_number()} ({sga.get_space_group_symbol()})")
  
  structure_2 = sga.get_refined_structure()
  sga2 = SpacegroupAnalyzer(structure_2, tol)
  if not no_primitive:
    structure_2 = sga2.get_primitive_standard_structure()

  atoms = AseAtomsAdaptor.get_atoms(structure_2)
  atoms.center()
  atoms.wrap()

  atoms.write(output, format='aims')
  print(f"Symmetrized geometry written to {output}")
  write_metadata(sga, output, tol, no_primitive)

def write_metadata(sga, output, tol, no_primitive):
  with open(output, 'a') as f:
    f.write('\n# Symmetry Refinement Metadata\n')
    f.write(f"# Space Group Found: {sga.get_space_group_number()} ({sga.get_space_group_symbol()})\n")
    f.write(f'# Tolerance: {tol}\n')
    f.write(f'# Converted to Primitive Cell: {"No" if no_primitive else "Yes"}\n')

if __name__ == '__main__':
  symmetrization()
