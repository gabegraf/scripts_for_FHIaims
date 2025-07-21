#!/usr/bin/env python
import numpy as np
import sys

def read_forces(filename):
    """
    Reads forces from the custom format in the specified file.
    """
    forces = []
    with open(filename, 'r') as f:
        lines = f.readlines()

    idx = 0
    n_lines = len(lines)
    while idx < n_lines:
        line = lines[idx].strip()
        if line.isdigit():
            # Number signals start of a new group
            n = int(line)
            idx += 2  # Skip the following blank/parameter row
            group = []
            for _ in range(n):
                coords = [float(x) for x in lines[idx].strip().split()]
                group.append(coords)
                idx += 1
            forces.append(np.array(group))
        else:
            idx += 1

    return forces


def compute_displacements(file1, file2):
    """
    Computes displacements between two files with matching force groups.
    """
    forces1 = read_forces(file1)
    forces2 = read_forces(file2)

    displacements = []

    for group1, group2 in zip(forces1, forces2):
        group_displacements = np.linalg.norm(group1 - group2, axis=1)
        displacements.append(group_displacements)
    return displacements


def main():
    if len(sys.argv) < 3:
        print("Usage: python force_displacement.py <file1> <file2>")
        sys.exit(1)

    file1 = sys.argv[1]
    file2 = sys.argv[2]

    displacements = compute_displacements(file1, file2)

    for i, group in enumerate(displacements, start=1):
        print(f"Group {i} displacements:")
        for j, d in enumerate(group, start=1):
            print(f"  Force {j}: {d}")

if __name__ == "__main__":
    main()

