#!/usr/bin/env python3
import os
import re
import subprocess
import csv

output_csv = "ExtractTotalEnergies.csv"
pattern = re.compile(r"displacement-(\d+)")

data = []

# Loop through directories
for dirname in os.listdir("."):
    match = pattern.fullmatch(dirname)
    if match and os.path.isdir(dirname):
        displacement_num = int(match.group(1))
        aims_out = os.path.join(dirname, "aims.out")

        if not os.path.exists(aims_out):
            print(f"⚠️ Skipping {dirname}, no aims.out found")
            continue

        try:
            # Run the external script and capture output
            result = subprocess.run(
                ["s23_extract_energy.py", "aims.out"],
                cwd=dirname,
                capture_output=True,
                text=True,
                check=True
            )
            energy = float(result.stdout.strip())
            data.append((displacement_num, energy))
            print(f"✓ {dirname}: {energy}")
        except Exception as e:
            print(f"⚠️ Failed to extract energy from {dirname}: {e}")

# Sort by displacement number
data.sort(key=lambda x: x[0])

# Write CSV
with open(output_csv, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["displacement", "energy"])
    writer.writerows(data)

print(f"\nSaved results to {output_csv}")

