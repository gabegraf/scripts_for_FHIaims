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

def symmetrization(output, tol):
  atoms = read('geometry.in', format='aims')
  structure = AseAtomsAdaptor.get_structure(atoms)
  sga = SpacegroupAnalyzer(structure, tol)
  print(f"Space Group Found! It's {sga.get_space_group_number()} ({sga.get_space_group_symbol()})")

  structure_2 = sga.get_refined_structure()
  sga2 = SpacegroupAnalyzer(structure_2, tol)
  structure_3 = sga2.find_primitive()

  atoms = AseAtomsAdaptor.get_atoms(structure_3)
  atoms.center()
  atoms.wrap()

  atoms.write(output, format='aims')
  print(f"Symmetrized geometry written to {output}")

if __name__ == '__main__':
  symmetrization()
