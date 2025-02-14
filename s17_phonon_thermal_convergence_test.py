#!/usr/bin/env python3
import yaml
import numpy as np
import matplotlib.pyplot as plt
import argparse

# Use non-GUI backend
import matplotlib
matplotlib.use("Agg")  # Ensures compatibility with headless systems

# Function to load YAML data from a file and extract thermal properties
def load_thermal_data(filename):
    with open(filename, "r") as file:
        data = yaml.safe_load(file)
    thermal_properties = data["thermal_properties"]
    return np.array([[entry["temperature"], entry["free_energy"], entry["entropy"], entry["heat_capacity"], entry["energy"]]
                     for entry in thermal_properties])

# Set up argument parser
parser = argparse.ArgumentParser(description="Plot and save thermal properties from two YAML files.")
parser.add_argument("file1", type=str, help="First YAML file")
parser.add_argument("file2", type=str, help="Second YAML file")
parser.add_argument("--output", type=str, default="thermal_plots.png", help="Output filename for the figure (default: thermal_plots.png)")
args = parser.parse_args()

# Load data from the input files
data1 = load_thermal_data(args.file1)
data2 = load_thermal_data(args.file2)

# Extract columns
T1, F1, S1, C1, E1 = data1.T  # Transpose to get separate arrays
T2, F2, S2, C2, E2 = data2.T

J_to_ev = (1.602176634 * 10 ** -19) ** -1
mol_to_atom = (6.02214076 * 10 ** 23)

F1 = F1 * J_to_ev * 1000 / mol_to_atom
F2 = F2 * J_to_ev * 1000 / mol_to_atom
S1 = S1 * J_to_ev / mol_to_atom
S2 = S2 * J_to_ev / mol_to_atom
E1 = E1 * J_to_ev * 1000 / mol_to_atom
E2 = E2 * J_to_ev * 1000 / mol_to_atom

# Create subplots
fig, axes = plt.subplots(2, 2, figsize=(12, 8))
labels = ["Free Energy (eV/unit cell)", "Entropy (eV/K/unit cell)", "Heat Capacity (J/K/mol)", "Energy (eV/unit cell)"] # mol is mol of unit cells, not of atoms
y_values = [(F1, F2), (S1, S2), (C1, C2), (E1, E2)]
colors = ["blue", "red"]
files = [args.file1, args.file2]

# Iterate through subplots
for ax, (y1, y2), label in zip(axes.flat, y_values, labels):
    # Filter data to keep only values where 250 ≤ T ≤ 500
    mask1 = (T1 >= 250) & (T1 <= 500)
    mask2 = (T2 >= 250) & (T2 <= 500)

    T1_filtered, y1_filtered = T1[mask1], y1[mask1]
    T2_filtered, y2_filtered = T2[mask2], y2[mask2]

    # Plot the filtered data
    ax.plot(T1_filtered, y1_filtered, label="4 x 4 x 2", color=colors[0])
    ax.plot(T2_filtered, y2_filtered, label="5 x 5 x 3", color=colors[1], linestyle="dashed")

    # Set labels
    ax.set_xlabel("Temperature (K)")
    ax.set_ylabel(label)
    ax.set_xlim(250, 500)

    # Dynamically adjust the y-axis limits based on filtered data
    ymin = min(y1_filtered.min(), y2_filtered.min())  # Find the lowest value in the filtered range
    ymax = max(y1_filtered.max(), y2_filtered.max())  # Find the highest value in the filtered range
    ax.set_ylim(ymin, ymax)

    ax.legend()

# Save the figure
plt.tight_layout()
plt.savefig(args.output, dpi=300, bbox_inches="tight")
print(f"Plot saved as {args.output}")

