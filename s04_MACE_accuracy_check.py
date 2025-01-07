#!/usr/bin/env python3
import os
import numpy as np
import csv
import matplotlib.pyplot as plt

def parse_geometry_file(filename):
    """Parse a geometry file to extract atomic coordinates."""
    atom_coords = []
    with open(filename, 'r') as file:
        for line in file:
            if line.startswith('atom_frac'):
                parts = line.split()
                coords = list(map(float, parts[1:4]))
                atom_coords.append(coords)
    return np.array(atom_coords)

def calculate_displacements(coords1, coords2):
    """Calculate displacements between two sets of coordinates."""
    displacements = np.linalg.norm(coords1 - coords2, axis=1)
    avg_displacement = np.mean(displacements)
    std_displacement = np.std(displacements)
    return avg_displacement, std_displacement, displacements

def plot_displacement_data(results, overall_avg):
    """Plot the average displacement data."""
    structures = [f"Structure {i}" for i, _, _ in results]
    avg_displacements = [avg for _, avg, _ in results]
    std_devs = [std for _, _, std in results]

    # Create the bar chart
    plt.figure(figsize=(10, 6))
    plt.bar(structures, avg_displacements, yerr=std_devs, capsize=5, color='skyblue', edgecolor='black', label="Average Displacement")
    plt.axhline(y=overall_avg, color='red', linestyle='--', label="Overall Average")
    
    # Labeling
    plt.title("Average Atomic Displacement for Each Structure", fontsize=14)
    plt.xlabel("Structures", fontsize=12)
    plt.ylabel("Average Displacement (Å)", fontsize=12)
    plt.xticks(rotation=45, ha='right')
    plt.legend()
    plt.tight_layout()

    # Save and show the plot
    plt.savefig("accuracy_check_plot.png")

def main():
    folder1 = "04_light_DFT_mace_converged_structures"
    folder2 = "08_MACE_accuracy_check"
    subfolder_template = "{}_structure_light_relax"
    filename = "geometry.in.next_step"

    all_displacements = []
    pair_results = []

    # Text file output
    with open("accuracy_check.txt", "w") as txt_out, open("accuracy_check.csv", "w", newline="") as csv_out:
        # Set up CSV writer
        csv_writer = csv.writer(csv_out)
        csv_writer.writerow(["Structure", "Average Displacement", "Standard Deviation"])

        for i in range(1, 17):
            subfolder1 = os.path.join(folder1, subfolder_template.format(i))
            subfolder2 = os.path.join(folder2, subfolder_template.format(i))
            file1 = os.path.join(subfolder1, filename)
            file2 = os.path.join(subfolder2, filename)

            if not os.path.exists(file1) or not os.path.exists(file2):
                msg = f"Missing geometry files for structure {i}. Skipping...\n"
                txt_out.write(msg)
                print(msg.strip())
                continue

            # Parse geometry files
            coords1 = parse_geometry_file(file1)
            coords2 = parse_geometry_file(file2)

            if coords1.shape != coords2.shape:
                msg = f"Mismatch in number of atoms for structure {i}. Skipping...\n"
                txt_out.write(msg)
                print(msg.strip())
                continue

            # Calculate displacements
            avg_disp, std_disp, displacements = calculate_displacements(coords1, coords2)
            all_displacements.extend(displacements)

            # Store results for the pair
            pair_results.append((i, avg_disp, std_disp))
            txt_out.write(f"Structure {i}: Avg displacement = {avg_disp:.6f}, Std deviation = {std_disp:.6f}\n")
            csv_writer.writerow([i, avg_disp, std_disp])
            print(f"Structure {i}: Avg displacement = {avg_disp:.6f}, Std deviation = {std_disp:.6f}")

        # Overall statistics
        overall_avg_disp = np.mean(all_displacements)
        overall_std_disp = np.std(all_displacements)
        txt_out.write("\nOverall Results:\n")
        txt_out.write(f"Average displacement: {overall_avg_disp:.6f}\n")
        txt_out.write(f"Standard deviation: {overall_std_disp:.6f}\n")
        print("\nOverall Results:")
        print(f"Average displacement: {overall_avg_disp:.6f}")
        print(f"Standard deviation: {overall_std_disp:.6f}")

        # Write overall results to CSV
        csv_writer.writerow([])
        csv_writer.writerow(["Overall", overall_avg_disp, overall_std_disp])

    # Plot the results
    plot_displacement_data(pair_results, overall_avg_disp)

if __name__ == "__main__":
    main()

