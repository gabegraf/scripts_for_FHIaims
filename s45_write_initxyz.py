#!/usr/bin/env python3

import sys
import numpy as np

input = "geometry.in"
output = "init.xyz"

if len(sys.argv) == 1:
  print("Defaulting to opening geometry.in and writing init.xyz")
elif len(sys.argv) == 3:
  input = sys.argv[1]
  output = sys.argv[2]
else:
  print("Usage: python3 script.py <input_filename> <output_filename>")
  sys.exit(1)

print("Reading file", input)
print("writing file", output)

data = []
cell = []
try:
  with open(input, "r") as f:
    for line in f:
      parts = line.split()
      if not parts:
        continue

      if parts[0] != "atom" and parts[0] != "atom_frac" and parts[0] != "lattice_vector":
        continue
      elif parts[0] == "lattice_vector":
          if len(parts) != 4:
              print(f"Skipping line with insufficient data: {line.strip()}")
              continue
          try:
              x_comp = parts[1]
              y_comp = parts[2]
              z_comp = parts[3]
              cell.append((x_comp, y_comp, z_comp))
          except ValueError:
              print(f"Skipping line with invalid numerical data: {line.strip()}")

      elif parts[0] == "atom_frac":
        print("Please convert atom_frac to atom first!")
        sys.exit(1)
      elif parts[0] == "atom":
        if len(parts) < 5:
          print(f"Skipping line with insufficient data: {line.strip()}")
          continue
        try:
          x = float(parts[1])
          y = float(parts[2])
          z = float(parts[3])
          species = parts[4]
          data.append((x, y, z, species))
        except ValueError:
          print(f"Skipping line with invalid numerical data: {line.strip()}")
except FileNotFoundError:
  print(f"Error: {input} not found!")
  sys.exit(1)
except Exception as e:
  print(f"An unexpected error occurred: {e}")
  sys.exit(1)

count = len(data)
comment_line = "#  Initial atomic positions in angstroms"

if len(cell) != 3:
    print("Error: You should have 3 lattice vectors, not {len(cell)}")
    sys.exit(1)
else:
    print(f"[ {cell[0][0]}, {cell[1][0]}, {cell[2][0]}, {cell[0][1]}, {cell[1][1]}, {cell[2][1]}, {cell[0][2]}, {cell[1][2]}, {cell[2][2]} ]")
    if float(cell[0][1]) != 0.0 or float(cell[0][2]) != 0.0 or float(cell[1][2]) != 0.0:
        print("Make sure your cell matrix is upper diagonal first!")
        sys.exit(1)

try:
  with open("init.xyz", "w") as f:
    f.write(f"{count}\n")
    f.write(f"{comment_line}\n")
    for x, y, z, species in data:
      f.write(f"{species} {x} {y} {z}\n")
  print(f"Successfully wrote data to {output}!")
  sys.exit(0)
except Exception as e:
  print(f"Error writing to {output}: {e}")
  sys.exit(1)
