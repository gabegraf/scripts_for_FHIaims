#!/usr/bin/env python3

import csv
import argparse

def process_file(input_file, output_file):
    # Read the file and process the data
    with open(input_file, 'r') as infile, open(output_file, 'w', newline='') as outfile:
        reader = csv.reader(infile)
        writer = csv.writer(outfile)

        # Skip the header row
        next(reader)

        # Collect positive energy states
        positive_states = []
        for row in reader:
            state = int(row[0])
            x_value = float(row[1])
            energy = float(row[2])
            if energy > 0:
                positive_states.append((state, x_value, energy))

        # Write only the rows corresponding to the first two states with positive energy
        if len(positive_states) >= 2:
            # Determine the first two positive states
            lowest_state = positive_states[0][0]
            next_state = positive_states[1][0]

            # Collect matching rows
            infile.seek(0)  # Reset the input file reader to process all rows again
            next(reader)  # Skip the header row again
            matching_rows = []
            for row in reader:
                state = int(row[0])
                if state in {lowest_state, next_state}:
                    matching_rows.append(row)

            # Sort rows by state (first column)
            matching_rows.sort(key=lambda x: int(x[0]))

            # Write sorted rows to the output file
            writer.writerow(["State", "x-value", "y-value"])
            writer.writerows(matching_rows)

def main():
    parser = argparse.ArgumentParser(description="Process spin splitting files.")
    parser.add_argument(
        "numbers",
        nargs="+",  # Accept multiple numbers as arguments
        type=int,
        help="Numbers to replace '*' in the file name pattern, e.g., 1001 1002.",
    )
    args = parser.parse_args()

    for number in args.numbers:
        input_file = f"stupid_file_for_bandmlk{number}_spin_splitting.csv"
        output_file = f"spin_splitting_input_for_bandmlk{number}.csv"
        print(f"Processing {input_file} -> {output_file}")
        process_file(input_file, output_file)

if __name__ == "__main__":
    main()

