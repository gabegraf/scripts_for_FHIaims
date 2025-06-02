#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np

# Use precise conversion factor: 1 cm⁻¹ = 0.0299792458 THz
CM_TO_THz = 0.0299792458

plt.rcParams.update({
    "font.family": "serif",
    "font.size": 12,
    "axes.labelsize": 12,
    "axes.titlesize": 14,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "legend.fontsize": 12,
    "figure.dpi": 300
})

file_name = "aims.RAMAN"
written_file = "aims_cm.RAMAN"

mode = []
frequency_cm = []
zpe = []
intensity = []

# Read and parse input file
with open(file_name, "r") as file:
    for line in file:
        if line.strip().startswith("#") or line.strip() == "":
            continue
        parts = line.split()
        if len(parts) == 4:
            mode.append(float(parts[0]))
            frequency_cm.append(float(parts[1]))
            zpe.append(float(parts[2]))
            intensity.append(float(parts[3]))

# Convert to NumPy arrays
frequency_cm = np.array(frequency_cm)
#frequency_thz = frequency_cm * CM_TO_THz
intensity = np.array(intensity)

# Sort by intensity (descending)
sorted_indices = np.argsort(intensity)[::-1]
frequency_cm_sorted = frequency_cm[sorted_indices]
intensity_sorted = intensity[sorted_indices]

# Plotting
plt.figure(figsize=(10, 6))
plt.bar(frequency_cm_sorted, intensity_sorted, width=3, color="#3366cc", edgecolor="black", zorder=2)

plt.xlabel("Frequency (cm⁻¹)")
plt.ylabel("Raman Intensity (a.u.)")
plt.title("Simulated Raman Spectrum")
plt.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.7, zorder=1)

plt.xlim(left=0)
plt.ylim(bottom=0,top=4000)

plt.tight_layout()
plt.savefig("raman_spectrum.png", dpi=300)
plt.show()

# Write converted data to new file
with open(written_file, "w") as out_file:
    out_file.write("# Mode    Frequency (cm⁻¹)    ZPE (eV)    Intensity\n")
    for m, f_cm, z, i in zip(mode, frequency_cm, zpe, intensity):
        f_thz = f_cm * CM_TO_THz
        out_file.write(f"{m:8.1f}  {f_cm:12.6f}  {z:12.6f}  {i:12.6f}\n")
