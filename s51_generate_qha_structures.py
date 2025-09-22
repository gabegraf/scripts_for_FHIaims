from ase.io import read, write
import pymatgen
from pymatgen.io.ase import AseAtomsAdaptor
import numpy as np
import click

@click.command()
@click.argument("input", type=click.Path(exists=True))
@click.option("--min_strain", default=-0.05)
@click.option("--max_strain", default=0.05)
@click.option("--n_outputs", default=10)
def create_qha_structures(input, min_strain, max_strain, n_outputs):
  atoms = read(input, format='aims')
  structure = AseAtomsAdaptor.get_structure(atoms)

  strain_list = np.linspace(min_strain, max_strain, n_outputs) #these are percentages around 0 e.g. -0.05 to 0.05
  with open("generated_strains.csv", 'w') as f:
    f.write("Percent of Volume, Volume, Geometry File\n")
    for i, value in enumerate(strain_list):
      strained_structure = structure.apply_strain(value, inplace=False)
      volume = strained_structure.volume
      strained_atoms = AseAtomsAdaptor.get_atoms(strained_structure)
      write(f"geometry.in-{i:03d}", strained_atoms, format='aims')
      f.write(f"{1+value}, {volume}, geometry.in-{i:03d}\n")

if __name__ == '__main__':
  create_qha_structures()
