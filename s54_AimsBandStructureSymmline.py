#!/usr/bin/env python3
import re
import numpy as np
import os
import json
from typing import Self
from monty.json import MontyEncoder, MontyDecoder
from pymatgen.io.aims.inputs import AimsGeometryIn
from pymatgen.electronic_structure.bandstructure import LobsterBandStructureSymmLine, Kpoint
from pymatgen.io.aims.outputs import AimsOutput
from pymatgen.io.aims.parsers import (
    get_lines,
    get_header_chunk,
    get_aims_out_chunks,
    AimsOutCalcChunk
)
from pymatgen.electronic_structure.core import Spin

class AimsBandStructureSymmline(LobsterBandStructureSymmLine):
    """
    Subclass of LobsterBandStructureSymmLine to add AIMS-specific methods if needed.
    """

    def ion_index_to_atom_map(self, structure) -> dict[int, int]:
        """
        Input: Pymatgen structure
        Output: Dictionary mapping ion index to its atomic species in the structure

        Ion index is defined as:
        ion_index = 2 * atom_index + spin_offset
        where spin_offset is 0 for spin up and 1 for spin down.

        Thus, there are two ion indices per atom index.

        This dictionary maps ion indices to their corresponding species in the structure.
        0 -> 'H' (for atom index 0, spin up)
        1 -> 'H' (for atom index 0, spin down)
        2 -> 'O' (for atom index 1, spin up)
        3 -> 'O' (for atom index 1, spin down)
        """
        atom_species, unique_species = get_species(structure)
        ion_index_to_atom_map = {}
        for atom_index, species in enumerate(atom_species):
            ion_index_to_atom_map[2 * atom_index] = species      # spin up
            ion_index_to_atom_map[2 * atom_index + 1] = species  # spin down
        return ion_index_to_atom_map

    @classmethod
    def from_bandmlk(cls, mlk_filename, geometry_filename="geometry.in", aims_out_filename="aims.out") -> Self:
        """
        Parse a bandmlk file and return eigenvalues and projections already as
        {Spin: ndarray(...)} structures expected by AimsBandStructureSymmline.

        This does a single pass to collect entries, determines array sizes, then
        allocates numpy arrays and fills them.
        """
        structure = AimsGeometryIn.from_file(geometry_filename).structure
        reciprocal_lattice = structure.lattice.reciprocal_lattice
        efermi = process_aimsout(aims_out_filename)

        match = re.match(r'bandmlk([12])(\d{3})\.out', mlk_filename)
        if not match:
            raise ValueError(f"Filename {mlk_filename} does not match expected pattern for bandmlk files.")
        spin_channel = int(match.group(1))

        kpoint_list = []

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
                    kpoint_list.append(frac_k_coord)
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
            raise ValueError(f"No band states found in {mlk_filename}.")

        # build a stable ordering of unique states and map to contiguous rows
        state_list = sorted(state_set)
        state_to_idx = {s: i for i, s in enumerate(state_list)}
        n_states = len(state_list)
        nk = len(kpoint_list)

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
        eigenvals = {spin_key: eig_arr}
        projections = {spin_key: proj_arr}
        kpt_coords = np.array(kpoint_list)
        return cls(
                    kpt_coords,
                    eigenvals,
                    reciprocal_lattice,
                    efermi,
                    {},                     # labels_dict
                    coords_are_cartesian=False,
                    structure=structure,
                    projections=projections,
                )

    @classmethod
    def from_json(cls, json_filename) -> Self:
        """
        Open a AimsBandStructureSymmline object from a JSON file.
        """
        try:
            with open(json_filename) as f:
                data = json.load(f, cls=MontyDecoder)
            if isinstance(data, cls):
                return data
            if isinstance(data, dict):
                return cls.from_dict(data)
            else:
                raise ValueError(f"JSON file {json_filename} does not contain a valid AimsBandStructureSymmline object.")
        except Exception as e:
            raise ValueError(f"Error reconstructing AimsBandStructureSymmline from dict: {e}")
                   
        
    def to_json(self, json_filename="symmline.json") -> None:
        """
        Save a AimsBandStructureSymmline object to a JSON file.
        Accepts either a dict mapping keys->objects or a single object.
        """
        dict_rep = self.as_dict() if hasattr(self, "as_dict") else self
        try:
            with open(json_filename, "w") as f:
                json.dump(dict_rep, f, indent=2, cls=MontyEncoder)
            print(f"Saved symmline to {json_filename}")
        except Exception as e:
            print(f"Error saving symmline to {json_filename}: {e}")

    @classmethod
    def load_all_symmlines_from_bandmlk(cls, geometry_filename="geometry.in", aims_out_filename="aims.out") -> dict[str, AimsBandStructureSymmline]:
        """
        Load all bandmlk files in the current directory and return a dictionary
        mapping the index to AimsBandStructureSymmline objects.

        Example:
        bandmlk1002.out, bandmlk1008.out -> {'1002': AimsBandStructureSymmline object, '1008': AimsBandStructureSymmline object}
        """
        symmlines = {}
        for filename in os.listdir("."):
            if re.match(r'bandmlk\d{4}\.out', filename):
                match = re.match(r'bandmlk(\d{4})\.out', filename)
                index = match.group(1)
                symmline = AimsBandStructureSymmline.from_bandmlk(filename, geometry_filename, aims_out_filename)
                symmlines[index] = symmline
        return symmlines

    @classmethod
    def write_all_symmlines_to_json(cls, symmlines: dict[str, AimsBandStructureSymmline]) -> None:
        """
        Input: Dictionary mapping symmline keys to AimsBandStructureSymmline objects
        Output: None
        Saves each symmline to a JSON file named symmline_{key}.json
        """
        for key, symmline in symmlines.items():
            json_filename = f"symmline_{key}.json"
            symmline.to_json(json_filename)
    
    @classmethod
    def load_all_symmlines_from_json(cls, geometry_filename="geometry.in", aims_out_filename="aims.out") -> dict[str, AimsBandStructureSymmline]:
        """
        Load all symmline JSON files in the current directory and return a dictionary
        mapping the index to AimsBandStructureSymmline objects.

        Example:
        symmline_1002.json, symmline_1008.json -> {'1002': AimsBandStructureSymmline object, '1008': AimsBandStructureSymmline object}
        """
        symmlines = {}
        for filename in os.listdir("."):
            if re.match(r'symmline_\d{4}\.json', filename):
                match = re.match(r'symmline_(\d{4})\.json', filename)
                index = match.group(1)
                symmline = AimsBandStructureSymmline.from_json(filename)
                symmlines[index] = symmline
        return symmlines

    @classmethod
    def loading_logic(cls, geometry_file, outfile, symmline, write_jsons=True) -> Self | dict[str, Self]:
        if symmline is not None:
            json_filename = f"symmline_{symmline}.json"
            if os.path.exists(json_filename):
                print(f"Found existing symmline JSON file: {json_filename}! Loading...")
                symmline_obj = AimsBandStructureSymmline.from_json(json_filename)
                print(f"Loaded symmline from {json_filename}")
                return symmline_obj
            else:
                print(f"No existing symmline JSON file: {json_filename}. Generating from bandmlk file...")
                mlk_filename = f"bandmlk{symmline}.out"
                symmline_obj = AimsBandStructureSymmline.from_bandmlk(mlk_filename, geometry_file, outfile)
                symmline_obj.to_json(json_filename)
                print(f"Generated and saved symmline to {json_filename}")
            return symmline_obj
        else:
            if any(re.match(r'symmline_\d{4}\.json', f) for f in os.listdir(".")):
                print("There is at least on symmline JSON file in the current directory.")
                print("Since no specific symmline was provided, we will load all of them as a dictionary.")
                symmlines = AimsBandStructureSymmline.load_all_symmlines_from_json(geometry_file, outfile)
                return symmlines
            else:
                print("No symmline JSON files found.")
                print("Processing all bandmlk files in the current directory to generate symmlines.")
                symmlines = AimsBandStructureSymmline.load_all_symmlines_from_bandmlk(geometry_file, outfile)
                print(f"Processed all bandmlk files to generate symmlines.")
                if write_jsons:
                    print("Saving all symmlines to JSON files...")
                    AimsBandStructureSymmline.write_all_symmlines_to_json(symmlines)
                    print(f"Done!")
                return symmlines

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
    efermi = aimsout.fermi_energy
    #lines = get_lines("aims.out")
    #header_chunk = get_header_chunk(lines)
    #calc_chunks = get_aims_out_chunks(lines, header_chunk)
    return efermi
