#!/usr/bin/env python
import numpy as np

def parse_geometry_file(geometry_file):
  
    lattice_vectors = []
    atomic_coords = []
    
    with open(f"{geometry_file}", "r") as file:
        lines = file.readlines()

        for line in lines:
            parts = line.split()
            if not line.strip():
                continue
            elif parts[0] == "lattice_vector":
                lattice_vectors.append([float(x) for x in parts[1:4]])
            elif parts[0] == "atom":
                atomic_coords.append(parts)


    # print(lattice_vectors)
    # print(atomic_coords)

    return lattice_vectors, atomic_coords

def calculate_distance_with_mic(atom_1, atom_2, lattice_vectors):

    lattice_matrix = np.array(lattice_vectors, dtype = np.float64)
    inv_lattice_matrix = np.linalg.inv(lattice_vectors)
    

    cart_diff_pre_mic = np.array(atom_2[1:4], dtype = np.float64) - np.array(atom_1[1:4], dtype = np.float64)
    frac_diff = np.dot(inv_lattice_matrix, cart_diff_pre_mic)
    
    if debug:
        print(f"Lattice Matrix: \n{lattice_matrix}")
        print(f"Inverse Lattice Matrix: \n{inv_lattice_matrix}")
        print(f"Cartesian Difference Before MIC: {cart_diff_pre_mic}")
        print(f"Fractional Difference Before MIC: {frac_diff}")

    frac_diff = frac_diff - np.round(frac_diff, 0)
    cart_diff_post_mic = np.dot(lattice_matrix, frac_diff)
    distance = np.linalg.norm(cart_diff_post_mic)

    if debug:
        print(f"Fractional Difference After MIC: {frac_diff}")
        print(f"Cartesian Difference After MIC: {cart_diff_post_mic}")
        print(f"Calculated Distance: {distance}")

    return distance

def calculate_distance_with_supercell(atom_1, atomic_coords, lattice_vectors):
    
    supercell_coords = []
    
    for x in [-1, 0, 1]:
        for y in [-1, 0, 1]:
            for z in [-1, 0, 1]:
                for atom in atomic_coords:
                    # Shift atom by (x, y, z) in the supercell
                    shift_vector = shift_vector = np.dot(x, lattice_vectors[0]) + np.dot(y, lattice_vectors[1]) + np.dot(z, lattice_vectors[2])
                    shifted_atom = np.array(atom[1:4], dtype=np.float64) + shift_vector
                    new_atom = [atom[0]] + list(shifted_atom) + [atom[4]]  # Include species
                    supercell_coords.append(new_atom)

    # Now calculate the distance between atom_1 and every other atom in the supercell
    min_distance = float('inf')
    closest_atom = None

    atom_1_coords = np.array(atom_1[1:4], dtype=np.float64)  # Coordinates of atom_1

    for atom in supercell_coords:
        if atom[4] != "H":
            atom_coords = np.array(atom[1:4], dtype=np.float64)  # Coordinates of the atom in the supercell
        
            # Calculate the distance using Euclidean norm
            distance = np.linalg.norm(atom_coords - atom_1_coords)

            # Update if the distance is smaller than the current minimum distance
            if distance < min_distance:
                min_distance = distance
                closest_atom = atom[4]  # Species of the closest atom
                closest_coords = atom[1:4]

    atom_1.append(closest_atom)
    atom_1.append(closest_coords)
    atom_1.append(min_distance)

    return min_distance, closest_atom, closest_coords

def find_nearest_non_H_atom(atom_1, atomic_coords, lattice_vectors):
    
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

    if debug:
        print(f"Minimum Distance: {min_distance}")
        print(f"Closest Species: {closest_species}")
        print(f"Closest Atom Coordinates: {closest_atom_coords}")
    
    return min_distance, closest_species, closest_atom_coords

def check_bond_type_agreement(initial_geometry, final_geometry, tolerance=1e-4):
 
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
    
    with open("bond_length_differences.out", 'w') as file:
        file.write("Bond Length Difference Data for C-H and N-H Bonds\n")

        for bond_type, differences in bond_length_differences.items():
            avg_difference = np.mean(differences)
            std_difference = np.std(differences)
            count = len(differences)

            file.write(f"{bond_type} Bonds: {count}\n")
            file.write(f"{bond_type} Average Length Difference: {avg_difference:.4f} Å\n")
            file.write(f"{bond_type} Standard Deviation: {std_difference:.4f} Å\n\n")

            # Print to console
            print(f"{bond_type} Length Difference (Å) = {avg_difference:.4f} ± {std_difference:.4f}")

    print(f"Bond length difference data written to bond_length_differences.out")

# main program
debug = False
lattice_vectors, atomic_coords = parse_geometry_file("initial_geometry.in")

for idx in atomic_coords:
    if idx[4] == "H":
        #find_nearest_non_H_atom(idx, atomic_coords, lattice_vectors)
        calculate_distance_with_supercell(idx, atomic_coords, lattice_vectors)

initial_geometry_cleaned = atomic_coords
print(initial_geometry_cleaned)
print("Initial geometry finished, now loading final geometry. \n")

lattice_vectors, atomic_coords = parse_geometry_file("final_geometry.in")

for idx in atomic_coords:
    if idx[4] == "H":
        #find_nearest_non_H_atom(idx, atomic_coords, lattice_vectors)
        calculate_distance_with_supercell(idx, atomic_coords, lattice_vectors)
        
final_geometry_cleaned = atomic_coords
print(final_geometry_cleaned)

check_bond_type_agreement(initial_geometry_cleaned, final_geometry_cleaned)
calculate_bond_length_changes(initial_geometry_cleaned, final_geometry_cleaned)
