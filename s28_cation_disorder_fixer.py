#!/usr/bin/env python
import os
import glob
import argparse
import pandas as pd

parser = argparse.ArgumentParser(description="Open a CIF file in the current working directory.")
parser.add_argument("cif_file", nargs="?", help="Specify the CIF file you would like to open")
args = parser.parse_args()

cif_files = glob.glob("*.cif")

if len(cif_files) == 0:
	print("There are no CIF files in the current working directory")
	exit(1)
elif len(cif_files) == 1:
	cif_file = cif_files[0]
else:
	if args.cif_file:
		if args.cif_file in cif_files:
			cif_file = args.cif_file
			print(f"Now opening {cif_file}")
		else:
			print(f"No CIF file found in the working directory called {args.cif_file}")
			exit(1)
	else: 
		print("Multiple CIF files found in current working directory, please pass one as an argument.")
		exit(1)

df = {}

with open(cif_file, "r") as file:

    for line_number, line in enumerate(file, start=0):
        if "atom_site_disorder_group" in line:
            disorder_group_number = line_number
            break

    for i in range(label_number -1, 0, -1):
        if "loop_" in line:
            loop_line_number = i
            break

    occupancy_index = loop_line_number - disorder_group_number

	for line_number, line in enumerate(file, start=loop_line_number):
		if not line.strip().startswith("_"):
			first_atom_line_number = line_number
			break
		else:
			continue

	number_of_attributes = ###################################################################################
		
    for line_number, line in enumerate(file, start=disorder_group_number):
        if line.strip() == "":
            continue
        elif line.strip().startswith("_"):
            continue
        else:
            row = line.strip().split()
            # Use the first element as the key and the rest as the value
            key = row[0]
            value = row[1:]  # All the other elements
            df[key] = value
