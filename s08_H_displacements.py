#!/usr/bin/env python

########## H DISPLACEMENTS AND BOND LENGTH DIFFERENCES ##############

#####################################################################

# If you run a hydrogen-only DFT relaxation, this script will
# tell you how far your hydrogens moved and how much their 
# bond lengths to their nearest atoms changed.

#####################################################################

# IMPORTANT:
# Currently, this script only works with cartesian coordinate
# input files (although I hope to extend it to fractional soon).
# To convert fractional geometry.in files to catesian coordinates,
# use the frac2cart function in the "utilities" directory of FHI-aims.
# 
# The script requires two files to be present in the directory that 
# it is run from, an input geometry (ie. a starting structure before
# H-relaxation) and a final geometry from after the H-relaxation has
# been done. You can provide these file names as arguments when you 
# run the script if you prefer:
#
# s08_H_displacements.py samplename1.in samplename2.in
#
# If no arguments are given, the script will look for the files:
# 
# initial_geometry.in
# final_geometry.in
# 
# This script uses the minimum image convention (MIC) to treat the 
# periodic boundary conditions of the system. The mathematics are
# discussed in the script body.

####################################################################

# NOT IMPORTANT:
# I originally got stuck while making the MIC method and tried a 
# separate method in which I did not use fractional coordinates,
# but instead made a supercell of cartesian coordinates. It is
# significantly slower and I ran into some problems with hydrogens
# moving through the unit cell walls (even though this should not 
# be a problem). I then found the issue with this MIC script,
# (I needed to transpose the lattice vector matrix),
# and therefore, continued the development of this one. The supercell
# method and the MIC method seem to agree on systems that the
# supercell method performs well on. As such, the MIC method (which
# is faster and not restricted) is the only one present here.

####################################################################

import numpy as np
import argparse
import math

# Read the geometry from the geometry.in file and break it into parts
def parse_geometry_file(geometry_file):

    lattice_vectors = []
    atomic_coords = []

    with open(f"{geometry_file}", "r") as file:
        lines = file.readlines()

        for line in lines:
            parts = line.split()
            if not line.strip(): # ignores empty lines
                continue
            elif parts[0] == "lattice_vector":
                lattice_vectors.append([float(x) for x in parts[1:4]])
            elif parts[0] == "atom":
                atomic_coords.append(parts)

    return lattice_vectors, atomic_coords

def ensure_all_atoms_in_unit_cell(atomic_coords, lattice_vectors):
    
    # This works by transforming the cartesian coordinates to fractional
    # coordinates, and moves every coordinate greater that 1 or less than
    # 0 to its appropriate value 0 =< x =< 1. While I don't think you'll 
    # find coordinates outside of 0 and 1 in the input file, the hydrogens
    # that are near the edge can definitely move outside of 0 and 1 
    # during the relaxation.
    #
    # Technically, I think that the MIC method can deal with this, but
    # I'd need to think a little harder about it. 


    # It's important to note that the intial matrix must be transposed,
    # so that the lattice matrix represents the correct transformation.

    lattice_matrix = (np.array(lattice_vectors, dtype=np.float64)).T
    inv_lattice_matrix = np.linalg.inv(lattice_matrix)

    for atom in atomic_coords:
        cartesian_coords = np.array(atom[1:4], dtype=np.float64)

        frac_atom = np.dot(inv_lattice_matrix, cartesian_coords)

        for i, frac_coord in enumerate(frac_atom):
            if frac_coord > 1.0 or frac_coord < 0.0:
                frac_atom[i] = frac_coord - np.floor(frac_coord)
                print(f"Adjusted {atom[4]} atom fractional coordinate {frac_coord} to {frac_atom[i]}")

    atom[1:4] = np.dot(lattice_matrix, frac_atom)
    print("All atoms are in the unit cell")
    return atomic_coords

def calculate_distance_with_mic(atom_1, atom_2, lattice_vectors):

    lattice_matrix = (np.array(lattice_vectors, dtype=np.float64)).T
    inv_lattice_matrix = np.linalg.inv(lattice_matrix)

    cart_diff_pre_mic = np.array(atom_2[1:4], dtype = np.float64) - np.array(atom_1[1:4], dtype = np.float64)
    frac_diff = np.dot(inv_lattice_matrix, cart_diff_pre_mic)

    frac_diff = frac_diff - np.round(frac_diff, 0)
    cart_diff_post_mic = np.dot(lattice_matrix, frac_diff)
    distance = np.linalg.norm(cart_diff_post_mic)

    return distance

def find_nearest_non_H_atom(atom_1, atomic_coords, lattice_vectors):

    # There is an implict assumption that there are no H-H bonds in 
    # the structure. This is reasonable every system this script is 
    # for, but make sure to verify this is correct for yours.

    min_distance = float('inf')
    closest_species = None
    closest_atom_coords = None

    for atom in atomic_coords:
        if atom[4] != "H":
            distance = calculate_distance_with_mic(atom_1, atom, lattice_vectors)

            if distance < min_distance:
                min_distance = distance
                closest_species = atom[4]
                closest_atom_coords = atom[1:4]

    atom_1.append(closest_species)
    atom_1.append(closest_atom_coords)
    atom_1.append(min_distance)

    return min_distance, closest_species, closest_atom_coords

def check_bond_type_agreement(initial_geometry, final_geometry, tolerance=1e-4):

    # There is an implicit assumption that no atoms move so much that
    # bond to a different atom now. This is reasonable for every system
    # this script is for, but make sure to verify this is correct for yours.
    # If it's not, the script will yell at you.

    for i, (atom_1, atom_2) in enumerate(zip(initial_geometry, final_geometry)):
        if len(atom_1) == len(atom_2) > 6:
            if atom_1[5] != atom_2[5]:
                print(f"Error: Atom {i} has bond type {atom_1[5]} in initial and {atom_2[5]} in final geometry!")

            # Check bond coordinate agreement within tolerance
            atom_1_vector = np.array(atom_1[6], dtype=np.float64)
            atom_2_vector = np.array(atom_2[6], dtype=np.float64)
            distance = np.linalg.norm(atom_1_vector - atom_2_vector)

            if distance > tolerance:
                print(f"Error: Atom {i} is bonded to {atom_1[5]} at {atom_1[6]} in the initial, "
                      f"but to {atom_2[5]} at {atom_2[6]} in the final geometry! "
                      f"Displacement: {distance:.4f} exceeds tolerance {tolerance:.4f}.")

    print("\nBond type agreement has been checked.")

def calculate_H_displacement_with_mic(initial_geometry, final_geometry, lattice_vectors):

    displacements = { "C-H": [], "N-H": []}

    for i, (atom_1, atom_2) in enumerate(zip(initial_geometry, final_geometry)):
        if atom_1[4] == atom_2[4] and atom_1[4] == "H":
            displacement = calculate_distance_with_mic(atom_1, atom_2, lattice_vectors)

            if atom_1[5] == "C":
                displacements["C-H"].append(displacement)
            elif atom_1[5] == "N":
                displacements["N-H"].append(displacement)
            else:
                print("Error: Wack-ass bond type")
    
    print(displacements["C-H"])
    with open("H_relaxation_postprocessing.out", 'w') as file:
        file.write("Hydrogen Atom Displacement Data for C-H and N-H Bonds\n")
        
        for bond_type, displacements in displacements.items():
            avg_displacement = np.mean(displacements)
            std_displacement = np.std(displacements)

            rounded_std_displacement = float(f"{std_displacement:.1g}")
            order = math.floor(math.log10(abs(std_displacement)))  # Order of magnitude of ref
            factor = 10**order  # Scaling factor
            rounded_avg_displacement = round(avg_displacement / factor) * factor  # Round and scale back

            count = len(displacements)

            file.write(f"{bond_type} Bonds: {count}\n")
            file.write(f"{bond_type} Average Displacement: {rounded_avg_displacement} ± {rounded_std_displacement} Å\n\n")

            # Print to console
            print(f"{bond_type} H atom Displacement (Å) = {rounded_avg_displacement} ± {rounded_std_displacement}")
            print(f"{bond_type} Bonds: {count}\n")

        print(f"H atom displacement data written to H_relaxation_postprocessing.out")

def calculate_bond_length_changes(initial_geometry_cleaned, final_geometry_cleaned):

    bond_length_differences = { "C-H": [], "N-H": []}

    for i, (atom_1, atom_2) in enumerate(zip(initial_geometry_cleaned, final_geometry_cleaned)):
        if atom_1[4] == atom_2[4] and atom_1[4] == "H":
            bond_length_change = atom_2[7] - atom_1[7]

            if atom_1[5] == "C":
                bond_length_differences["C-H"].append(bond_length_change)
            elif atom_1[5] == "N":
                bond_length_differences["N-H"].append(bond_length_change)
            else:
                print("Error: Wack-ass bond type")
    
    with open("H_relaxation_postprocessing.out", 'a') as file:
        file.write("Bond Length Difference Data for C-H and N-H Bonds\n")
        
        for bond_type, differences in bond_length_differences.items():
            avg_difference = np.mean(differences)
            std_difference = np.std(differences)
            
            rounded_std_difference = float(f"{std_difference:.1g}")
            order = math.floor(math.log10(abs(std_difference)))  # Order of magnitude of ref
            factor = 10**order  # Scaling factor
            rounded_avg_difference = round(avg_difference / factor) * factor  # Round and scale back
            
            count = len(differences)

            file.write(f"{bond_type} Bonds: {count}\n")
            file.write(f"{bond_type} Average Length Difference: {rounded_avg_difference} ± {rounded_std_difference} Å\n\n")

            # Print to console
            print(f"{bond_type} Length Difference (Å) = {rounded_avg_difference} ± {rounded_std_difference}")
            print(f"{bond_type} Bonds: {count}\n")

    print(f"Bond length difference data written to H_relaxation_postprocessing.out")

############################# MAIN PROGRAM ##########################

parser = argparse.ArgumentParser()

parser.add_argument(
        "initial_geometry",
        type=str,
        nargs="?",  # Makes it optional
        default="initial_geometry.in",  # Default value if not provided
        help="Argument 1: Enter the name of your geometry.in file from before H-relaxation."
    )

parser.add_argument(
        "final_geometry",
        type=str,
        nargs="?",  # Makes it optional
        default="final_geometry.in",  # Default value if not provided
        help="Argument 2: Enter the name of your geometry.in file from after H-relaxation."
    )

args = parser.parse_args()

#####################################################################

print("Initial geometry processing.")

lattice_vectors, atomic_coords_initial = parse_geometry_file(args.initial_geometry)
atomic_coords_initial = ensure_all_atoms_in_unit_cell(atomic_coords_initial, lattice_vectors)

for idx in atomic_coords_initial:
    if idx[4] == "H":
        find_nearest_non_H_atom(idx, atomic_coords_initial, lattice_vectors)

initial_geometry_cleaned = atomic_coords_initial
print("Initial geometry complete, final geometry processing.")

#####################################################################

lattice_vectors, atomic_coords_final = parse_geometry_file(args.final_geometry)
atomic_coords_final = ensure_all_atoms_in_unit_cell(atomic_coords_final, lattice_vectors)

for idx in atomic_coords_final:
    if idx[4] == "H":
        find_nearest_non_H_atom(idx, atomic_coords_final, lattice_vectors)

final_geometry_cleaned = atomic_coords_final
print("Final geometry complete.")

#####################################################################

check_bond_type_agreement(initial_geometry_cleaned, final_geometry_cleaned)
calculate_H_displacement_with_mic(initial_geometry_cleaned, final_geometry_cleaned, lattice_vectors)
calculate_bond_length_changes(initial_geometry_cleaned, final_geometry_cleaned)
