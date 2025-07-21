#!/usr/bin/env python3
import sys

def extract_section(input_file, start_marker, end_marker, section_name):
    try:
        with open(input_file, 'r') as f:
            content = f.read()
    except FileNotFoundError:
        sys.exit(f"Error: Input file '{input_file}' not found.")

    # Find start and end indices
    start_idx = content.find(start_marker)
    if start_idx == -1:
        sys.exit(f"Error: Start marker for {section_name} not found in {input_file}.")

    start_idx += len(start_marker)  # Move past the marker

    end_idx = content.find(end_marker, start_idx)
    if end_idx == -1:
        sys.exit(f"Error: End marker for {section_name} not found in {input_file}.")

    # Extract content between markers (excluding markers themselves)
    section_content = content[start_idx:end_idx].strip()

    # Remove leading whitespace from each line
    section_content = '\n'.join(line.lstrip() for line in section_content.split('\n'))

    return section_content

# Define the exact markers (as they appear in aims.out)
control_markers = {
    'start': """-----------------------------------------------------------------------
  Parsing control.in (first pass over file, find array dimensions only).
  The contents of control.in will be repeated verbatim below
  unless switched off by setting 'verbatim_writeout .false.' .
  in the first line of control.in .
  -----------------------------------------------------------------------""",
    'end': """-----------------------------------------------------------------------
  Completed first pass over input file control.in .
  -----------------------------------------------------------------------"""
}

geometry_markers = {
    'start': """-----------------------------------------------------------------------
  Parsing geometry.in (first pass over file, find array dimensions only).
  The contents of geometry.in will be repeated verbatim below
  unless switched off by setting 'verbatim_writeout .false.' .
  in the first line of geometry.in .
  -----------------------------------------------------------------------""",
    'end': """-----------------------------------------------------------------------
  Completed first pass over input file geometry.in .
  -----------------------------------------------------------------------"""
}

# Extract and write files
try:
    # Extract control.in
    control_content = extract_section(
        'aims.out',
        control_markers['start'],
        control_markers['end'],
        'control.in'
    )
    with open('control.in', 'w') as f:
        f.write(control_content)

    # Extract geometry.in
    geometry_content = extract_section(
        'aims.out',
        geometry_markers['start'],
        geometry_markers['end'],
        'geometry.in'
    )
    with open('geometry.in', 'w') as f:
        f.write(geometry_content)

    print("Successfully extracted control.in and geometry.in")
except Exception as e:
    sys.exit(f"Error: {str(e)}")
