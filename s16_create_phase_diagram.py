#!/usr/bin/env python3
import os
import argparse
import yaml
import numpy as np
import re
import csv

# Create an argument parser
parser = argparse.ArgumentParser(description="Process aims.out, YAML, and vib.out files.")
parser.add_argument("file_name", type=str, nargs='+', help="Path to the input files")

def read_aims_out(file_path):
    """Reads an aims.out file and extracts the total energy."""
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"File not found: {file_path}")

    with open(file_path, 'r') as f:
        lines = f.readlines()

    if len(lines) < 2 or "Have a nice day." not in lines[-2]:
        raise ValueError('The calculation did not finish successfully.')

    for line in reversed(lines):
        if "Total energy of the DFT / Hartree-Fock s.c.f. calculation" in line:
            try:
                match = re.search(r":\s+([-+]?\d*\.\d+|\d+)", line)
                if match:
                    return float(match.group(1))
            except Exception as e:
                raise ValueError(f"Error extracting energy: {e}")

    raise ValueError('Total energy line not found in the output file.')

def extract_temperature_free_energy(file_path):
    """Extracts temperature and free energy into an Nx2 NumPy array."""
    with open(file_path, "r") as file:
        data = yaml.safe_load(file)

    thermal_properties = data.get("thermal_properties", [])
    temp_energy_list = [(entry["temperature"], entry["free_energy"]) for entry in thermal_properties]
    return np.array(temp_energy_list)

def process_vib_out(file_path):
    """Processes a .vib.out file and writes vibrational data to a CSV."""
    header_elements = [
        "temperature [K]",
        "vibrational free energy [eV]",
        "rotational free energy [eV]",
        "pressure [Pa]",
        "translational free energy [eV]"
    ]
    header_pattern = r'\s+'.join(re.escape(element) for element in header_elements)
    number_pattern = r'-?\d+\.\d+'
    output_file = "rovibe.csv"

    try:
        with open(file_path, "r") as infile, open(output_file, "w", newline='') as outfile:
            found_header = False
            writer = csv.writer(outfile)

            for line in infile:
                line = line.replace("-Infinity", "-99999.99")

                if re.search(header_pattern, line):
                    if found_header:
                        continue
                    else:
                        found_header = True
                        writer.writerow(header_elements)
                elif "---------" in line:
                    continue
                else:
                    if found_header:
                        if "Leaving the Hessian diagonalizer." in line:
                            print("Completed first pass. Filling in missing values now.")
                            break
                        numbers = re.findall(number_pattern, line)
                        if len(numbers) == 5:
                            writer.writerow(numbers)
                        elif len(numbers) == 2:
                            writer.writerow(['', '', '', numbers[0], numbers[1]])
    except Exception as e:
        print(f"An error occurred: {e}")

    try:
        with open(output_file, "r+") as outfile:
            rows = list(csv.reader(outfile))
            header = rows[0]
            data = rows[1:]
            current_value = None

            for i in [0, 1, 2]:
                for row in data:
                    if re.match(number_pattern, row[i]):
                        current_value = float(row[i])
                    elif row[i] == '':
                        row[i] = str(current_value)

            outfile.seek(0)
            writer = csv.writer(outfile)
            writer.writerow(header)
            writer.writerows(data)
            print("Missing data added. Calculating total energetic contributions now.")

    except Exception as e:
        print(f"An error occurred: {e}")

    try:
        with open(output_file, "r+") as outfile:
            rows = list(csv.reader(outfile))
            header = rows[0]
            data = rows[1:]
            header.append('Total Energetic Contribution [eV]')

            for row in data:
                row.append(str(float(row[1]) + float(row[2]) + float(row[4])))

            outfile.seek(0)
            writer = csv.writer(outfile)
            writer.writerow(header)
            writer.writerows(data)
    except Exception as e:
        print(f"An error occurred: {e}")
    
    print("Finished. Take a look!")

    try:
        # Read the CSV to extract relevant columns
        data_array = []
        with open(output_file, 'r') as infile:
            reader = csv.reader(infile)
            header = next(reader)  # Skip the header

            for row in reader:
                # Extract the columns: temperature, pressure, and total energy
                try:
                    temperature = float(row[0])
                    pressure = float(row[3])
                    total_energy = float(row[5])
                    data_array.append([temperature, pressure, total_energy])
                except ValueError:
                    # Skip rows with missing or invalid data
                    continue

        # Convert the list to a numpy array
        return np.array(data_array)
    except Exception as e:
        print(f"An error occurred while reading the CSV: {e}")
    return np.array([])  # Return an empty array if there is an error

def calculate_total_free_energy(vibration_check, phonon_free_energies, energy):
    """Calculates the total free energy at each temperature and pressure point."""
    total_data = []

    for temperature, pressure, rovibrational_energy in vibration_check:
        # Find the corresponding temperature correction from the .yaml data
        yaml_free_energy = None
        for temp, free_energy in phonon_free_energies:
            if abs(temp - temperature) < 1e-3:  # Tolerance for temperature match
                yaml_free_energy = free_energy
                break

        if yaml_free_energy is None:
            raise ValueError(f"No matching temperature found in YAML for T = {temperature}")

        # Calculate the total free energy: f(T,P) + f(T) + internal energy
        total_free_energy = rovibrational_energy + yaml_free_energy + energy

        # Append the result
        total_data.append([temperature, pressure, rovibrational_energy, yaml_free_energy, energy, total_free_energy])

    # Convert the list to a numpy array and return
    return np.array(total_data)

# Example usage within the main script
if __name__ == "__main__":
    args = parser.parse_args()

    # Initialize these variables outside the loop
    energy = None
    phonon_free_energies = None
    vibration_check = None

    # Loop through each file provided
    for file_path in args.file_name:
        print(f"Processing file: {file_path}")

        if file_path.endswith("aims.out"):
            energy = read_aims_out(file_path)
            print(f"Total Energy: {energy} eV")
        elif file_path.endswith(".yaml"):
            phonon_free_energies = extract_temperature_free_energy(file_path)
            print("Extracted Temperature & Free Energy Data:\n", phonon_free_energies)
        elif file_path.endswith(".vib.out"):
            vibration_check = process_vib_out(file_path)
            print("Vibration Data:", vibration_check)
        else:
            print(f"Error: Unsupported file type for {file_path}")

    # Calculate the total free energy array only after all files have been processed
    if energy is not None and phonon_free_energies is not None and vibration_check is not None:
        total_free_energy_array = calculate_total_free_energy(vibration_check, phonon_free_energies, energy)
        print("Total Free Energy Data:\n", total_free_energy_array)

        # Save the total free energy array to CSV
        with open("total_free_energy.csv", "w", newline='') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(["Temperature [K]", "Pressure [Pa]", "Rovibrational Free Energy [eV]",
                             "YAML Free Energy [eV]", "Internal Energy [eV]", "Total Free Energy [eV]"])
            writer.writerows(total_free_energy_array)
    else:
        print("Error: Missing data. Ensure all necessary files are provided.")

