#!/usr/bin/env python3
"""
AIMS Output Analysis Script

This script analyzes aims.out files in displacement-* directories to extract 
and compute statistics on self-consistency cycles, SCF reinitializations, 
and total CPU times.
"""

import os
import re
import glob
import statistics
from datetime import datetime


def natural_sort_key(text):
    """
    Generate a key for natural sorting of displacement directories.
    Converts 'displacement-001' to a sortable format.
    """
    # Extract the number from displacement-XXX
    match = re.search(r'displacement-(\d+)', text)
    if match:
        return int(match.group(1))
    return 0


def process_aims_file(file_path, cycles_list, reinit_list, times_list):
    """
    Process a single aims.out file to extract required values.
    
    Args:
        file_path: Path to the aims.out file
        cycles_list: List to append self-consistency cycles to
        reinit_list: List to append SCF reinitializations to
        times_list: List to append CPU times to
    
    Raises:
        Exception: If file cannot be read or required patterns not found
    """
    found_cycles = False
    found_reinit = False
    found_time = False
    
    # Read file from bottom to top
    try:
        with open(file_path, 'r') as file:
            all_lines = file.readlines()
    except Exception as e:
        raise Exception(f"Could not read file: {e}")
    
    # Process lines from bottom to top
    for line in reversed(all_lines):
        # Search for self-consistency cycles (if not found yet)
        if not found_cycles:
            cycles_match = re.search(r'\|\s*Number of self-consistency cycles\s*:\s*(\d+)', line)
            if cycles_match:
                cycles_value = int(cycles_match.group(1))
                cycles_list.append(cycles_value)
                found_cycles = True
        
        # Search for SCF reinitializations (if not found yet)
        if not found_reinit:
            reinit_match = re.search(r'\|\s*Number of SCF \(re\)initializations\s*:\s*(\d+)', line)
            if reinit_match:
                reinit_value = int(reinit_match.group(1))
                reinit_list.append(reinit_value)
                found_reinit = True
        
        # Search for total time (if not found yet)
        if not found_time:
            time_match = re.search(r'\|\s*Total time\s*:\s*(\d+\.\d+)\s*s\s*\d+\.\d+\s*s', line)
            if time_match:
                cpu_time = float(time_match.group(1))
                times_list.append(cpu_time)
                found_time = True
        
        # If all values found, break early
        if found_cycles and found_reinit and found_time:
            break
    
    # Check if all required values were found
    if not found_cycles:
        raise Exception("Self-consistency cycles pattern not found")
    if not found_reinit:
        raise Exception("SCF reinitializations pattern not found")
    if not found_time:
        raise Exception("Total time pattern not found")


def main():
    """Main function to process all displacement directories and analyze results."""
    
    # Initialize data storage
    self_consistency_cycles = []
    scf_reinitializations = []
    cpu_times = []
    processed_files_count = 0
    error_count = 0
    
    # Find all displacement directories matching pattern
    directory_pattern = "displacement-*"
    all_directories = glob.glob(directory_pattern)
    
    # Sort directories to process in numerical order
    all_directories.sort(key=natural_sort_key)
    
    print(f"Found {len(all_directories)} displacement directories")
    
    # Process each directory
    for directory in all_directories:
        aims_file_path = os.path.join(directory, "aims.out")
        
        # Check if aims.out exists
        if not os.path.exists(aims_file_path):
            print(f"ERROR: aims.out not found in {directory}")
            error_count += 1
            continue
        
        # Try to process the file
        try:
            process_aims_file(aims_file_path, self_consistency_cycles, scf_reinitializations, cpu_times)
            processed_files_count += 1
            print(f"Successfully processed: {directory}")
        except Exception as e:
            print(f"ERROR processing {aims_file_path}: {e}")
            error_count += 1
    
    # Print processing summary
    print(f"\nProcessing Summary:")
    print(f"Successfully processed: {processed_files_count} files")
    print(f"Errors encountered: {error_count} files")
    
    # Check if we have data to analyze
    if len(self_consistency_cycles) == 0:
        print("No valid data found to analyze")
        return
    
    # Calculate statistics for self-consistency cycles
    cycles_mean = statistics.mean(self_consistency_cycles)
    cycles_stdev = statistics.stdev(self_consistency_cycles) if len(self_consistency_cycles) > 1 else 0
    
    # Calculate statistics for SCF reinitializations
    reinit_mean = statistics.mean(scf_reinitializations)
    reinit_stdev = statistics.stdev(scf_reinitializations) if len(scf_reinitializations) > 1 else 0
    
    # Calculate statistics for CPU times
    time_mean = statistics.mean(cpu_times)
    time_stdev = statistics.stdev(cpu_times) if len(cpu_times) > 1 else 0
    total_time_seconds = sum(cpu_times)
    total_time_minutes = total_time_seconds / 60
    
    # Output results to console
    print("\n" + "="*60)
    print("ANALYSIS RESULTS")
    print("="*60)
    
    print(f"\nSelf-Consistency Cycles:")
    print(f"  Average: {cycles_mean:.2f}")
    print(f"  Standard Deviation: {cycles_stdev:.2f}")
    
    print(f"\nSCF Reinitializations:")
    print(f"  Average: {reinit_mean:.2f}")
    print(f"  Standard Deviation: {reinit_stdev:.2f}")
    
    print(f"\nCPU Times:")
    print(f"  Average: {time_mean:.2f} seconds")
    print(f"  Standard Deviation: {time_stdev:.2f} seconds")
    print(f"  Total Time: {total_time_seconds:.2f} seconds ({total_time_minutes:.2f} minutes)")
    
    # Save results to file
    output_filename = "aims_analysis_results.txt"
    current_timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    
    try:
        with open(output_filename, 'w') as output_file:
            output_file.write("AIMS Output Analysis Results\n")
            output_file.write(f"Generated on: {current_timestamp}\n")
            output_file.write("="*60 + "\n\n")
            
            output_file.write("Processing Summary:\n")
            output_file.write(f"Successfully processed: {processed_files_count} files\n")
            output_file.write(f"Errors encountered: {error_count} files\n\n")
            
            output_file.write("Self-Consistency Cycles:\n")
            output_file.write(f"  Average: {cycles_mean:.2f}\n")
            output_file.write(f"  Standard Deviation: {cycles_stdev:.2f}\n\n")
            
            output_file.write("SCF Reinitializations:\n")
            output_file.write(f"  Average: {reinit_mean:.2f}\n")
            output_file.write(f"  Standard Deviation: {reinit_stdev:.2f}\n\n")
            
            output_file.write("CPU Times:\n")
            output_file.write(f"  Average: {time_mean:.2f} seconds\n")
            output_file.write(f"  Standard Deviation: {time_stdev:.2f} seconds\n")
            output_file.write(f"  Total Time: {total_time_seconds:.2f} seconds ({total_time_minutes:.2f} minutes)\n")
        
        print(f"\nResults saved to: {output_filename}")
        
    except Exception as e:
        print(f"ERROR: Could not save results to file: {e}")


if __name__ == "__main__":
    main()
