import re
import numpy as np
import os
from pymatgen.io.aims.input import AimsGeometryIn
from pymatgen.electronic_structure.bandstructure import BandStructureSymmLine, Kpoint
from pymatgen.io.aims.outputs import AimsOutput
from pymatgen.io.aims.parsers import (
    get_lines,
    get_header_chunk,
    get_aims_out_chunks,
    AimsOutCalcChunk
)
import click
from pymatgen.electronic_structure.core import Spin

# If there is mulliken data, I want to process it from the bandmlk<symmline_index>.out files
# but if not I can pull the energies from the band<symmline_index>.out files. Let's start with mulliken data.

#eigenvals is a dictionary {spin: 2D array [band_index][kpoint_index]}
# BandStructureSymmLine(kpoints: ArrayLike, eigenvals: Mapping[Spin, ArrayLike], lattice: Lattice, efermi: float, labels_dict: Mapping[str, Kpoint], coords_are_cartesian: bool = False, structure: Structure | None = None, projections: Mapping[Spin, NDArray] | None = None)
@click.command()
@click.option("--geometry_file", default="geometry.in", help="Path to the AIMS geometry file")
@click.option("--outfile", default="aims.out", help="Path to the AIMS output file")
def main(geometry_file, outfile):
    """
    This script reads bandmlk files in the current directory, processes their contents
    and builds BandStructureSymmLine objects for each symmline and spin channel.
    """
    structure = AimsGeometryIn.from_file(geometry_file).structure
    reciprocal_lattice = structure.lattice.reciprocal_lattice
    efermi = process_aimsout(outfile)
    kpoints_on_symmlines = find_k_points(structure, reciprocal_lattice)
    symmlines = build_symmlines_from_bandmlk_files(kpoints_on_symmlines, structure, reciprocal_lattice, efermi)
    print(f"Processed {len(symmlines)} symmlines from bandmlk files.")

def build_symmlines_from_bandmlk_files(kpoints_on_symmlines, structure, lattice, efermi) -> list[BandStructureSymmLine]:
    """
    This function scans the current directory for bandmlk files, reads their contents,
    and builds BandStructureSymmLine objects for each symmline and spin channel.
    """
    symmlines = []
    for filename in os.listdir("."):
        match = re.match(r'bandmlk([12])(\d{3})\.out', filename)
        if not match:
            continue
        spin_channel = int(match.group(1))
        symmline_index = int(match.group(2))
        kpoints_on_symmline = kpoints_on_symmlines.get(symmline_index, [])
        if not kpoints_on_symmline:
            continue
        mlk_filename = filename
        eigenvals, projections = read_mlk(mlk_filename, kpoints_on_symmline, structure, spin_channel)

        # build labels_dict from labeled kpoints (if any)
        labels_dict = {kp.label: kp for kp in kpoints_on_symmline if getattr(kp, "label", None)}

        symmline = BandStructureSymmLine(
            kpoints_on_symmline,
            eigenvals,
            lattice,
            efermi,
            labels_dict,
            coords_are_cartesian=False,
            structure=structure,
            projections=projections,
        )
        symmlines.append(symmline)

    return symmlines

def read_mlk(mlk_filename, kpoints_on_symmline, structure, spin_channel) -> tuple[dict, dict]:
    """
    Parse a bandmlk file and return eigenvalues and projections already as
    {Spin: ndarray(...)} structures expected by BandStructureSymmLine.

    This does a single pass to collect entries, determines array sizes, then
    allocates numpy arrays and fills them.
    """
    kpoint_list = kpoints_on_symmline
    nk = len(kpoint_list)

    entries = []
    state_set = set()
    max_orb = 0
    max_ion = 0

    # Read file and collect relevant data rows
    with open(mlk_filename, "r") as f:
        current_k = None
        current_state = None
        for raw in f:
            line = raw.strip()

            k_match = re.match(r"k point number:\s*(\d+):\s*\(\s*([\-\d\.Ee+]+)\s+([\-\d\.Ee+]+)\s+([\-\d\.Ee+]+)\s*\)", line)
            if k_match:
                if current_state is not None:
                    current_state = None
                k_idx = int(k_match.group(1)) - 1
                frac_k_coord = np.array([float(x) for x in k_match.groups()[1:]])
                kpoint = kpoint_list[k_idx]
                if not np.allclose(frac_k_coord, kpoint.frac_coords):
                    raise ValueError("K-point fractional coordinates do not match!")
                current_k = k_idx
                continue

            state_match = re.match(r"State\s+(\d+)", line)
            if state_match:
                current_state = int(state_match.group(1))
                state_set.add(current_state)
                continue

            if current_state is not None and current_k is not None and line:
                parts = line.split()
                # first_parts = parts (State line already handled); parts format: eigenvalue occ.number atom spin total s p d f ...
                eigenvalue = float(parts[1])
                atom_number = int(parts[3])
                line_spin = int(parts[4])
                ion_index = 2 * atom_number - (2 if line_spin == 1 else 1)
                orbital_projections = [float(x) for x in parts[6:]]

                # update maxima
                max_orb = max(max_orb, len(orbital_projections))
                max_ion = max(max_ion, ion_index + 1)
                entries.append((current_k, current_state, eigenvalue, ion_index, orbital_projections))

    if len(state_set) == 0:
        return {}, {}

    # build a stable ordering of unique states and map to contiguous rows
    state_list = sorted(state_set)
    state_to_idx = {s: i for i, s in enumerate(state_list)}
    n_states = len(state_list)

    # allocate arrays (bands x kpoints) and (bands x kpoints x orb x ion)
    eig_arr = np.full((n_states, nk), np.nan, dtype=float)
    proj_arr = np.zeros((n_states, nk, max_orb, max_ion), dtype=float)

    # fill arrays using state->row mapping
    for k_idx, state, eigenvalue, ion_index, orbital_projections in entries:
        state_idx = state_to_idx[int(state)]
        eig_arr[state_idx, int(k_idx)] = eigenvalue
        for orb_idx, val in enumerate(orbital_projections):
            proj_arr[state_idx, int(k_idx), orb_idx, int(ion_index)] = val

    spin_key = Spin.up if int(spin_channel) == 1 else Spin.down
    return {spin_key: eig_arr}, {spin_key: proj_arr}

def find_k_points(structure, reciprocal_lattice) -> dict[int, list[Kpoint]]:
    """
    Input: Pymatgen structure and its reciprocal lattice
    Output: Dictionary mapping symmline index to list of kpoints along that symmline

    Example if the G -> X symmline is index 1 and has 10 kpoints, then
    kpoints_on_symmlines[1] = [Kpoint1, Kpoint2, ..., Kpoint10]

    Indices are base zero.
    """
    kpoints_on_symmlines = {}
    for filename in os.listdir("."):
        match = re.match(r'band([12])(\d{3})\.out', filename)
        if match:
            spin_channel = int(match.group(1))
            symmline_index = int(match.group(2))
            kpoints_on_symmlines[symmline_index] = []
            with open(filename, 'r') as f:
                for line in f:
                    line = line.strip().split()
                    frac_k_coords = np.array([float(x) for x in line[1:4]])
                    kpoints_on_symmlines[symmline_index].append(
                    Kpoint(
                        frac_k_coords, reciprocal_lattice, to_unit_cell=False,
                        coords_are_cartesian=False, label=None
                    ))
    return kpoints_on_symmlines

def get_species(structure) -> tuple [list[str], list[str]]:
    """
    Input: Pymatgen structure
    Output: (atom_species, unique_species)
    where atom_species is a list of species symbols for each atom in the structure in order
    and unique_species is a list of unique species symbols in the structure

    Example:
    atom_species = ['H', 'O', 'H']
    unique_species = ['H', 'O']
    """
    atom_species = [x.symbol for x in structure.species]
    unique_species = [x.symbol for x in structure.elements]
    print(f"Species retrieved from the Pymatgen structure!")
    return atom_species, unique_species

def process_aimsout(outfile) -> float:
    """
    Input: Path to AIMS output file
    Output: Fermi energy as a float
    
    This will likely be expanded.
    """
    aimsout = AimsOutput.from_outfile(outfile)
    efermi = aimsout.efermi
    #lines = get_lines("aims.out")
    #header_chunk = get_header_chunk(lines)
    #calc_chunks = get_aims_out_chunks(lines, header_chunk)
    return efermi

if __name__ == "__main__":
    main()