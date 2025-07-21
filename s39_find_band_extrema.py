#!/usr/bin/env python3
"""
Band Structure Analyzer - Finds absolute extrema along k-path
"""

import numpy as np
import matplotlib.pyplot as plt
import os
import re
import csv
import sys
from scipy.interpolate import CubicSpline
from scipy.optimize import root_scalar
from typing import Dict, List, Tuple, Optional

# Constants
TOLERANCE = 1e-5
DERIV_TOL = 1e-2  # Threshold for flat derivative
OCCUPANCY_VALUES = {0.0, 1.0}
SUMMARY_CSV = "band_analysis_summary.csv"

def find_band_files() -> List[str]:
    """Find and sort band{number}.out files"""
    pattern = re.compile(r'band(\d+)\.out')
    matched_files = []
    for f in os.listdir('.'):
        match = pattern.fullmatch(f)
        if match:
            try:
                matched_files.append((int(match.group(1)), f))
            except ValueError:
                continue
    if not matched_files:
        raise FileNotFoundError("No band{number}.out files found")
    return [f for _, f in sorted(matched_files)]

def parse_file(filename: str) -> Tuple[Dict[Tuple[float, float, float], Dict[str, float]], List[Tuple[float, float, float]]]:
    """Parse file with validation"""
    k_point_data = {}
    all_k_points = []
    
    with open(filename, 'r') as f:
        for line_idx, line in enumerate(f, 1):
            if not line.strip():
                continue
                
            elements = line.split()
            try:
                k_point = tuple(map(float, elements[1:4]))
                all_k_points.append(k_point)
                
                rest = list(map(float, elements[4:]))
                if len(rest) % 2 != 0:
                    raise ValueError(f"Line {line_idx}: Odd number of values")
                    
                occupancy = rest[::2]
                energy = rest[1::2]
                
                invalid_occ = [o for o in occupancy if o not in OCCUPANCY_VALUES]
                if invalid_occ:
                    raise ValueError(f"Line {line_idx}: Invalid occupancies {invalid_occ}")
                
                pt_homo = max(e for o, e in zip(occupancy, energy) if o == 1.0)
                pt_lumo = min(e for o, e in zip(occupancy, energy) if o == 0.0)
                
                k_point_data[k_point] = {'homo': pt_homo, 'lumo': pt_lumo}
                
            except (ValueError, IndexError) as e:
                raise ValueError(f"Line {line_idx}: {str(e)}")
    
    if not k_point_data:
        raise ValueError("No valid data in file")
    return k_point_data, all_k_points

def compute_k_distance(k_points: List[Tuple[float, float, float]]) -> np.ndarray:
    """Convert 3D k-points to 1D path distance"""
    distances = [0.0]
    for i in range(1, len(k_points)):
        delta = np.array(k_points[i]) - np.array(k_points[i-1])
        distances.append(distances[-1] + np.linalg.norm(delta))
    return np.array(distances)

def find_absolute_extrema(k_distances: np.ndarray, 
                         energies: np.ndarray,
                         extremum_type: str) -> Tuple[float, float]:
    """
    Finds absolute extrema, comparing critical points with endpoints
    
    Args:
        k_distances: Path distances
        energies: Band energies
        extremum_type: 'max' or 'min'
        
    Returns:
        (k_position, energy) of absolute extremum
    """
    spline = CubicSpline(k_distances, energies)
    derivative = spline.derivative(1)
    
    # Find all critical points
    critical_points = []
    try:
        # Search between all consecutive k-points
        for i in range(len(k_distances)-1):
            try:
                sol = root_scalar(derivative,
                                bracket=[k_distances[i], k_distances[i+1]],
                                method='brentq')
                if abs(derivative(sol.root)) < DERIV_TOL:
                    critical_points.append((sol.root, spline(sol.root)))
            except ValueError:
                continue
    except Exception as e:
        print(f"Warning: Critical point search failed: {str(e)}")
    
    # Always include endpoints
    critical_points.extend([
        (k_distances[0], energies[0]),
        (k_distances[-1], energies[-1])
    ])
    
    # Find absolute extremum
    if extremum_type == 'max':
        idx = np.argmax([e for _, e in critical_points])
    else:  # 'min'
        idx = np.argmin([e for _, e in critical_points])
    
    return critical_points[idx]

def analyze_single_band(filename: str,
                       global_homo: Dict,
                       global_lumo: Dict) -> Dict:
    """Process single band file"""
    try:
        k_point_data, all_k_points = parse_file(filename)
        k_distances = compute_k_distance(all_k_points)
        
        # Extract energies
        homo_energies = np.array([data['homo'] for data in k_point_data.values()])
        lumo_energies = np.array([data['lumo'] for data in k_point_data.values()])
        
        # Find absolute extrema
        homo_k, homo_energy = find_absolute_extrema(k_distances, homo_energies, 'max')
        lumo_k, lumo_energy = find_absolute_extrema(k_distances, lumo_energies, 'min')
        
        # Convert k-distance back to 3D k-point
        def distance_to_kpoint(distance):
            for i in range(1, len(k_distances)):
                if distance <= k_distances[i]:
                    t = (distance - k_distances[i-1])/(k_distances[i] - k_distances[i-1])
                    return tuple(np.array(all_k_points[i-1]) + t*(np.array(all_k_points[i]) - np.array(all_k_points[i-1])))
            return all_k_points[-1]
        
        homo_kpoint = distance_to_kpoint(homo_k)
        lumo_kpoint = distance_to_kpoint(lumo_k)
        
        # Validate on path
        k_start, k_end = all_k_points[0], all_k_points[-1]
        if not is_on_path(k_start, k_end, homo_kpoint):
            print(f"Warning: HOMO point not on k-path in {filename}")
        if not is_on_path(k_start, k_end, lumo_kpoint):
            print(f"Warning: LUMO point not on k-path in {filename}")
        
        # Update global extrema
        if homo_energy > global_homo["energy"]:
            global_homo.update({
                "energy": homo_energy,
                "k_point": homo_kpoint,
                "band_file": filename
            })
        if lumo_energy < global_lumo["energy"]:
            global_lumo.update({
                "energy": lumo_energy,
                "k_point": lumo_kpoint,
                "band_file": filename
            })
        
        # Generate plot
        plot_filename = plot_band_results(
            k_distances, homo_energies, lumo_energies,
            (homo_k, homo_energy), (lumo_k, lumo_energy),
            filename
        )
        
        return {
            "band_file": filename,
            "homo_energy": homo_energy,
            "homo_k_point": homo_kpoint,
            "lumo_energy": lumo_energy,
            "lumo_k_point": lumo_kpoint,
            "k_separation": np.linalg.norm(np.array(homo_kpoint) - np.array(lumo_kpoint)),
            "energy_gap": homo_energy - lumo_energy,
            "plot_file": plot_filename
        }
        
    except Exception as e:
        print(f"ERROR processing {filename}: {str(e)}", file=sys.stderr)
        raise

def is_on_path(k_start: Tuple[float, float, float], 
              k_end: Tuple[float, float, float], 
              k_test: Tuple[float, float, float]) -> bool:
    """Check if point lies on k-path segment"""
    k_start, k_end, k_test = map(np.array, (k_start, k_end, k_test))
    line_vec = k_end - k_start
    test_vec = k_test - k_start
    t = np.dot(test_vec, line_vec)/np.dot(line_vec, line_vec)
    return (0-TOLERANCE <= t <= 1+TOLERANCE and 
            np.linalg.norm(k_test - (k_start + t*line_vec)) < TOLERANCE)

def plot_band_results(k_distances, homo_energies, lumo_energies,
                     homo_critical, lumo_critical, filename):
    """Generate band structure plot"""
    plt.figure(figsize=(10,6))
    k_fine = np.linspace(k_distances[0], k_distances[-1], 500)
    
    # Plot splines
    homo_spline = CubicSpline(k_distances, homo_energies)
    lumo_spline = CubicSpline(k_distances, lumo_energies)
    plt.plot(k_fine, homo_spline(k_fine), 'b-', label='HOMO spline')
    plt.plot(k_fine, lumo_spline(k_fine), 'r-', label='LUMO spline')
    
    # Plot critical points
    plt.plot(homo_critical[0], homo_critical[1], 'g*', markersize=12, label='HOMO')
    plt.plot(lumo_critical[0], lumo_critical[1], 'm*', markersize=12, label='LUMO')
    
    plt.xlabel('k-path distance (1/Å)')
    plt.ylabel('Energy (eV)')
    plt.title(f'Band Structure: {filename}')
    plt.legend()
    plt.grid(True)
    
    plot_filename = f"{os.path.splitext(filename)[0]}_analysis.png"
    plt.savefig(plot_filename, dpi=300, bbox_inches='tight')
    plt.close()
    return plot_filename

def write_summary_csv(results, global_homo, global_lumo):
    """Write results to CSV"""
    with open(SUMMARY_CSV, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow([
            "Band", "HOMO (eV)", "HOMO k-point", 
            "LUMO (eV)", "LUMO k-point", 
            "k-Separation", "Gap (eV)", "Plot"
        ])
        for r in results:
            writer.writerow([
                r["band_file"],
                f"{r['homo_energy']:.6f}",
                str(r["homo_k_point"]),
                f"{r['lumo_energy']:.6f}",
                str(r["lumo_k_point"]),
                f"{r['k_separation']:.4f}",
                f"{r['energy_gap']:.4f}",
                r["plot_file"]
            ])
        writer.writerow([])
        writer.writerow(["Global HOMO:", f"{global_homo['energy']:.6f}", f"from {global_homo['band_file']}"])
        writer.writerow(["Global LUMO:", f"{global_lumo['energy']:.6f}", f"from {global_lumo['band_file']}"])
        writer.writerow(["Fundamental Gap:", f"{global_homo['energy']-global_lumo['energy']:.6f}"])

def main():
    """Main execution"""
    try:
        global_homo = {"energy": -np.inf, "k_point": None, "band_file": None}
        global_lumo = {"energy": np.inf, "k_point": None, "band_file": None}
        
        band_files = find_band_files()
        print(f"Found {len(band_files)} bands to analyze")
        
        results = []
        for band_file in band_files:
            results.append(analyze_single_band(band_file, global_homo, global_lumo))
            r = results[-1]
            print(f"\n{band_file}:")
            print(f"HOMO: {r['homo_energy']:.4f} eV at {r['homo_k_point']}")
            print(f"LUMO: {r['lumo_energy']:.4f} eV at {r['lumo_k_point']}")
            print(f"Separation: {r['k_separation']:.4f} 1/Å, Gap: {r['energy_gap']:.4f} eV")
        
        print("\n==== Global Results ====")
        print(f"Global HOMO: {global_homo['energy']:.4f} eV")
        print(f"  From: {global_homo['band_file']} at {global_homo['k_point']}")
        print(f"Global LUMO: {global_lumo['energy']:.4f} eV")
        print(f"  From: {global_lumo['band_file']} at {global_lumo['k_point']}")
        print(f"Fundamental Gap: {global_homo['energy']-global_lumo['energy']:.4f} eV")
        
        write_summary_csv(results, global_homo, global_lumo)
        print(f"\nSummary saved to {SUMMARY_CSV}")
        
    except Exception as e:
        print(f"FATAL ERROR: {str(e)}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
