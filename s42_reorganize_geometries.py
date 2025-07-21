#!/usr/bin/env python3
"""
Geometry File Organization Script

This script organizes displacement geometry files by matching them with reference files
and creating a properly structured directory hierarchy.
"""

import os
import shutil
import glob
import numpy as np
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import logging
from tqdm import tqdm

try:
    from ase.io import read
    from ase import Atoms
except ImportError:
    print("Error: ASE (Atomic Simulation Environment) is required but not installed.")
    print("Install it with: pip install ase")
    exit(1)

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class GeometryOrganizer:
    def __init__(self, correct_assignments_dir: str = "correct_assignments", 
                 output_dir: str = "correct_order_2"):
        self.correct_assignments_dir = Path(correct_assignments_dir)
        self.output_dir = Path(output_dir)
        self.reference_geometries = {}
        self.matched_pairs = []
        self.empty_directories = []
        self.unmatched_files = []
        
    def round_coordinates(self, atoms: Atoms, decimals: int = 10) -> np.ndarray:
        """Round atomic coordinates to specified decimal places."""
        return np.round(atoms.get_positions(), decimals)
    
    def load_reference_geometries(self) -> Dict[str, np.ndarray]:
        """Load all reference geometry files from correct_assignments directory."""
        logger.info(f"Loading reference geometries from {self.correct_assignments_dir}")
        
        if not self.correct_assignments_dir.exists():
            raise FileNotFoundError(f"Reference directory {self.correct_assignments_dir} not found")
        
        reference_files = list(self.correct_assignments_dir.glob("geometry.in-*"))
        logger.info(f"Found {len(reference_files)} reference files")
        
        reference_geometries = {}
        failed_loads = []
        
        for ref_file in tqdm(reference_files, desc="Loading reference files"):
            try:
                # Explicitly specify the format as 'aims' for FHI-aims geometry.in files
                atoms = read(ref_file, format='aims')
                rounded_coords = self.round_coordinates(atoms)
                # Extract number from filename (e.g., geometry.in-001 -> 001)
                file_number = ref_file.name.split('-')[1]
                reference_geometries[file_number] = rounded_coords
                logger.debug(f"Successfully loaded {ref_file} with {len(atoms)} atoms")
            except Exception as e:
                logger.error(f"Failed to load reference file {ref_file}: {e}")
                failed_loads.append(ref_file)
        
        if failed_loads:
            logger.warning(f"Failed to load {len(failed_loads)} reference files")
        
        logger.info(f"Successfully loaded {len(reference_geometries)} reference geometries")
        return reference_geometries
    
    def find_displacement_directories(self) -> List[Path]:
        """Find all directories matching the *-displacement pattern."""
        displacement_dirs = []
        
        # Look for directories matching the pattern
        for item in Path(".").iterdir():
            if item.is_dir() and item.name.endswith("-displacement"):
                displacement_dirs.append(item)
        
        # Sort by the numeric part for consistent processing
        displacement_dirs.sort(key=lambda x: x.name.split('-')[0])
        
        logger.info(f"Found {len(displacement_dirs)} displacement directories")
        return displacement_dirs
    
    def compare_geometries(self, coords1: np.ndarray, coords2: np.ndarray, 
                          tolerance: float = 1e-10) -> bool:
        """Compare two sets of coordinates with given tolerance."""
        if coords1.shape != coords2.shape:
            return False
        return np.allclose(coords1, coords2, atol=tolerance)
    
    def find_best_match(self, target_coords: np.ndarray, 
                       reference_geometries: Dict[str, np.ndarray],
                       first_guess: str = None) -> Optional[str]:
        """Find the best matching reference geometry, optionally trying first_guess first."""
        
        # Try the first guess if provided
        if first_guess and first_guess in reference_geometries:
            if self.compare_geometries(target_coords, reference_geometries[first_guess]):
                return first_guess
        
        # Search through all reference geometries
        for ref_id, ref_coords in reference_geometries.items():
            if self.compare_geometries(target_coords, ref_coords):
                return ref_id
        
        return None
    
    def process_displacement_directory(self, disp_dir: Path, 
                                     reference_geometries: Dict[str, np.ndarray]) -> Optional[Tuple[str, str]]:
        """Process a single displacement directory and find its match."""
        geometry_file = disp_dir / "geometry.in"
        
        # Check if geometry.in exists
        if not geometry_file.exists():
            logger.warning(f"No geometry.in found in {disp_dir}")
            self.empty_directories.append(str(disp_dir))
            return None
        
        # Check if directory is effectively empty (only contains geometry.in)
        dir_contents = list(disp_dir.iterdir())
        if len(dir_contents) == 0:
            logger.warning(f"Directory {disp_dir} is empty")
            self.empty_directories.append(str(disp_dir))
            return None
        
        try:
            # Load and round the geometry - explicitly specify format for FHI-aims
            atoms = read(geometry_file, format='aims')
            rounded_coords = self.round_coordinates(atoms)
            
            # Extract number from directory name for first guess
            dir_number = disp_dir.name.split('-')[0]
            # Pad to 3 digits to match reference file format
            first_guess = dir_number.zfill(3)
            
            # Find the best match
            match = self.find_best_match(rounded_coords, reference_geometries, first_guess)
            
            if match:
                return (str(disp_dir), match)
            else:
                logger.warning(f"No match found for {disp_dir}")
                self.unmatched_files.append(str(disp_dir))
                return None
                
        except Exception as e:
            logger.error(f"Error processing {disp_dir}: {e}")
            self.unmatched_files.append(str(disp_dir))
            return None
    
    def create_output_structure(self):
        """Create the correct_order directory structure."""
        if self.output_dir.exists():
            logger.info(f"Output directory {self.output_dir} already exists")
            response = input("Do you want to overwrite it? (y/n): ")
            if response.lower() != 'y':
                logger.info("Aborted by user")
                return False
            shutil.rmtree(self.output_dir)
        
        self.output_dir.mkdir(exist_ok=True)
        logger.info(f"Created output directory: {self.output_dir}")
        return True
    
    def copy_matched_directories(self):
        """Copy matched directories to the correct_order structure."""
        logger.info("Copying matched directories to correct_order structure")
        
        for original_dir, match_id in tqdm(self.matched_pairs, desc="Copying directories"):
            # Create new directory name: displacement + match_id (removing leading zeros)
            new_dir_name = f"displacement{match_id.lstrip('0') or '0'}"
            new_dir_path = self.output_dir / new_dir_name
            
            try:
                shutil.copytree(original_dir, new_dir_path)
                logger.debug(f"Copied {original_dir} -> {new_dir_path}")
            except Exception as e:
                logger.error(f"Failed to copy {original_dir} to {new_dir_path}: {e}")
    
    def print_summary(self):
        """Print a summary of the organization process."""
        print("\n" + "="*50)
        print("GEOMETRY ORGANIZATION SUMMARY")
        print("="*50)
        print(f"Total displacement directories found: {len(self.find_displacement_directories())}")
        print(f"Successfully matched: {len(self.matched_pairs)}")
        print(f"Empty or missing geometry.in: {len(self.empty_directories)}")
        print(f"Unmatched files: {len(self.unmatched_files)}")
        
        if self.empty_directories:
            print(f"\nEmpty directories:")
            for empty_dir in self.empty_directories:
                print(f"  - {empty_dir}")
        
        if self.unmatched_files:
            print(f"\nUnmatched files:")
            for unmatched in self.unmatched_files:
                print(f"  - {unmatched}")
        
        print(f"\nOutput directory created: {self.output_dir}")
        print("="*50)
    
    def run(self):
        """Main execution method."""
        logger.info("Starting geometry organization process")
        
        # Load reference geometries
        reference_geometries = self.load_reference_geometries()
        
        # Find displacement directories
        displacement_dirs = self.find_displacement_directories()
        
        if not displacement_dirs:
            logger.error("No displacement directories found")
            return
        
        # Process each displacement directory
        logger.info("Processing displacement directories")
        for disp_dir in tqdm(displacement_dirs, desc="Processing directories"):
            result = self.process_displacement_directory(disp_dir, reference_geometries)
            if result:
                self.matched_pairs.append(result)
        
        # Create output structure
        if not self.create_output_structure():
            return
        
        # Copy matched directories
        if self.matched_pairs:
            self.copy_matched_directories()
        
        # Print summary
        self.print_summary()
        
        logger.info("Geometry organization process completed")

def main():
    """Main entry point."""
    organizer = GeometryOrganizer()
    organizer.run()

if __name__ == "__main__":
    main()
