#!/bin/bash

# Default: print everything
error_only=false

# Parse flags
while getopts "e" opt; do
  case $opt in
    e) error_only=true ;;
    *) echo "Usage: $0 [-e]"; exit 1 ;;
  esac
done

# Find all directories and sort them
for dir in $(find . -type d | sort); do

    aims_out="$dir/aims.out"

    if [[ -f "$aims_out" ]]; then
        if grep -q "Have a nice day" "$aims_out"; then
            # Only print on success if not in error-only mode
            if ! $error_only; then
                echo "Success: 'Have a nice day' found in $aims_out"
            fi
        else
            # Print directories missing the phrase
            if ! $error_only; then
                echo "$dir"
            fi
        fi
    else
        # Print only if -e flag is active, or always otherwise
        if $error_only; then
            echo "$dir"
        else
            echo "Error: $aims_out does not exist"
        fi
    fi
done


