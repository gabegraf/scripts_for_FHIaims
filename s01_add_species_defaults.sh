#!/bin/bash
# Script: Append Species Defaults to control.in
# Description: This script appends species-specific default data to a `control.in` file based on the provided tightness level and species symbols. 
#              It validates the tightness level and species symbols, looks up the corresponding atomic numbers, and appends the relevant data from predefined files.
# Usage: ./s01_add_species_defaults <tightness> <species symbols...>
#        <tightness>: One of "light", "intermediate", "tight", or "really_tight".
#        <species symbols...>: One or more valid chemical symbols (e.g., H, He, C, O).

# Example: ./s01_add_species_defaults tight H O C
#          This appends the default data for Hydrogen (H), Oxygen (O), and Carbon (C) using the "tight" tightness level.

# Dependencies:
# - The script assumes the existence of a directory structure: $SPECIES_DEFAULTS/defaults_2020/<tightness>/<atomic_number>_<species>_default.
# - The `control.in` file must be writable by the script.

# Exit Codes:
# - 1: Invalid number of arguments, invalid tightness level, or unknown species symbol.

# Check that at least two arguments are provided
if [ "$#" -lt 2 ]; then
  echo "Usage: $0 <tightness> <species symbols...>"
  exit 1
fi

# Set up variables
tightness="$1"
shift  # Shift arguments to process species symbols

# Valid tightness levels
valid_tightness_levels=("light" "intermediate" "tight" "really_tight")

# Check if the provided tightness level is valid
if [[ ! " ${valid_tightness_levels[@]} " =~ " ${tightness} " ]]; then
  echo "Invalid tightness level. Options are: ${valid_tightness_levels[*]}"
  exit 1
fi

# Define an associative array for species atomic numbers
declare -A atomic_numbers=(
  [H]="01" [He]="02" [Li]="03" [Be]="04" [B]="05" [C]="06" [N]="07" [O]="08" [F]="09" [Ne]="10"
  [Na]="11" [Mg]="12" [Al]="13" [Si]="14" [P]="15" [S]="16" [Cl]="17" [Ar]="18" [K]="19" [Ca]="20"
  [Sc]="21" [Ti]="22" [V]="23" [Cr]="24" [Mn]="25" [Fe]="26" [Co]="27" [Ni]="28" [Cu]="29" [Zn]="30"
  [Ga]="31" [Ge]="32" [As]="33" [Se]="34" [Br]="35" [Kr]="36" [Rb]="37" [Sr]="38" [Y]="39" [Zr]="40"
  [Nb]="41" [Mo]="42" [Tc]="43" [Ru]="44" [Rh]="45" [Pd]="46" [Ag]="47" [Cd]="48" [In]="49" [Sn]="50"
  [Sb]="51" [Te]="52" [I]="53" [Xe]="54" [Cs]="55" [Ba]="56" [La]="57" [Hf]="72" [Ta]="73" [W]="74"
  [Re]="75" [Os]="76" [Ir]="77" [Pt]="78" [Au]="79" [Hg]="80" [Tl]="81" [Pb]="82" [Bi]="83" [Po]="84"
  [At]="85" [Rn]="86" [Fr]="87" [Ra]="88"
)

# Iterate over each species symbol provided as arguments
for species in "$@"; do
  # Lookup atomic number in associative array
  atomic_number="${atomic_numbers[$species]}"

  # Check if atomic number is found
  if [ -z "$atomic_number" ]; then
    echo "Unknown species symbol: $species"
    exit 1
  fi

  # Construct the command and execute it
  file_path="$SPECIES_DEFAULTS/defaults_2020/${tightness}/${atomic_number}_${species}_default"
  cat "$file_path" >> control.in
done

echo "Data appended to control.in successfully."

