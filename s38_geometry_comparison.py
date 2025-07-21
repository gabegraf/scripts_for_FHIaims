#!/usr/bin/env python3

import click
import ase
from ase.io import read
import numpy as np

def compare_lattice_vectors(atoms1, atoms2, tol=1e-8):
    lattice1 = atoms1.get_cell()
    lattice2 = atoms2.get_cell()
    return np.allclose(lattice1, lattice2, atol=tol)

def get_cell_properties(atoms):
    cell = atoms.get_cell()
    lengths = cell.lengths()
    volume = cell.volume
    angles = cell.angles()
    return lengths, volume, angles

def check_atoms_compatibility(atoms1, atoms2):
    if len(atoms1) != len(atoms2):
        raise ValueError("Number of atoms differ between the two structures.")
    symbols1 = atoms1.get_chemical_symbols()
    symbols2 = atoms2.get_chemical_symbols()
    if symbols1 != symbols2:
        raise ValueError("Atom order or species differ between the two structures.")

def compute_displacements(atoms1, atoms2, same_lattice):
    positions1 = atoms1.get_positions()
    positions2 = atoms2.get_positions()
    symbols = atoms1.get_chemical_symbols()

    if same_lattice:
        # Use PBC-aware distance calculation
        _, distances = ase.geometry.get_distances(
            p1=positions1,
            p2=positions2,
            cell=atoms1.get_cell(),
            pbc=atoms1.get_pbc()
        )
        # Get diagonal elements (distance between corresponding atoms)
        displacements = np.diag(distances)
    else:
        # Just Euclidean norm difference (no PBC correction)
        displacements = np.linalg.norm(positions2 - positions1, axis=1)

    return symbols, displacements

@click.command()
@click.argument('file1', type=click.Path(exists=True))
@click.argument('file2', type=click.Path(exists=True))
@click.option('--print-all', is_flag=True, help="Print displacement of all atoms.")
def main(file1, file2, print_all):
    try:
        atoms1 = read(file1, format='aims')
        atoms2 = read(file2, format='aims')
    except Exception:
        raise click.ClickException("Failed to read one or both geometry.in files.")

    same_lattice = compare_lattice_vectors(atoms1, atoms2)

    if same_lattice:
        click.echo("Lattice vectors are the same.")
        lengths, volume, angles = get_cell_properties(atoms1)
        click.echo("\nCell properties:")
        click.echo(f"Volume: {volume:.6f} Å³")
        click.echo(f"Lattice lengths (Å): {lengths}")
        click.echo(f"Angles (degrees): α={angles[0]:.3f}, β={angles[1]:.3f}, γ={angles[2]:.3f}")
    else:
        click.echo("Lattice vectors are different.")
        lengths1, volume1, angles1 = get_cell_properties(atoms1)
        lengths2, volume2, angles2 = get_cell_properties(atoms2)

        click.echo("\nCell 1 properties:")
        click.echo(f"Volume: {volume1:.6f} Å³")
        click.echo(f"Lattice lengths (Å): {lengths1}")
        click.echo(f"Angles (degrees): α={angles1[0]:.3f}, β={angles1[1]:.3f}, γ={angles1[2]:.3f}")

        click.echo("\nCell 2 properties:")
        click.echo(f"Volume: {volume2:.6f} Å³")
        click.echo(f"Lattice lengths (Å): {lengths2}")
        click.echo(f"Angles (degrees): α={angles2[0]:.3f}, β={angles2[1]:.3f}, γ={angles2[2]:.3f}")

        length_diff = np.array(lengths2) - np.array(lengths1)
        angle_diff = np.array(angles2) - np.array(angles1)
        volume_diff = volume2 - volume1

        click.echo("\nDifferences (Second - First):")
        click.echo(f"Volume difference: {volume_diff:.6f} Å³")
        click.echo(f"Lattice lengths difference (Å): {length_diff}")
        click.echo(f"Angles difference (degrees): α={angle_diff[0]:.3f}, β={angle_diff[1]:.3f}, γ={angle_diff[2]:.3f}")

    # Check atoms compatibility
    try:
        check_atoms_compatibility(atoms1, atoms2)
    except ValueError as e:
        raise click.ClickException(str(e))

    symbols, displacements = compute_displacements(atoms1, atoms2, same_lattice)

    # Collect displacement statistics per atom type
    unique_symbols = sorted(set(symbols))
    click.echo("\nDisplacement statistics per atom type:")
    for sym in unique_symbols:
        disp_sym = displacements[np.array(symbols) == sym]

        avg_disp = np.mean(disp_sym)
        std_disp = np.std(disp_sym)
        click.echo(f"  {sym}: {avg_disp:.6f} ± {std_disp:.6f} Å")

    # Atom with largest displacement
    max_idx = np.argmax(displacements)
    click.echo(f"\nAtom with largest displacement: index {max_idx} ({symbols[max_idx]}) moved {displacements[max_idx]:.6f} Å")

    # Optionally print all atomic displacements
    if print_all:
        click.echo("\nAll atomic displacements:")
        for i, (sym, disp) in enumerate(zip(symbols, displacements)):
            click.echo(f"  Atom {i} ({sym}): {disp:.6f} Å")

if __name__ == "__main__":
    main() 

