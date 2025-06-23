#!/usr/bin/env python3
"""
Script to process KofamScan output files and extract unique, sorted KO IDs
from significant hits (lines marked with an asterisk *).
"""

import os
import sys

class KofamOutputProcessor:
    """
    Class to process KofamScan output files and extract unique KO IDs
    from significant hits.
    """
    
    def __init__(self, input_dir='input_files', output_dir='significant_ko_lists'):
        """
        Initialize the processor with input and output directories.
        
        Args:
            input_dir: Directory containing KofamScan output files
            output_dir: Directory to save processed files in a significant_ko_lists folder 
        """
        self.input_dir = input_dir
        self.output_dir = output_dir
        
    def process_directory(self):
        """
        Process all files in the input directory and generate output files
        with unique, sorted KO IDs.
        """
        # Ensure output directory exists
        os.makedirs(self.output_dir, exist_ok=True)
        
        # Get list of files in input directory
        input_files = [f for f in os.listdir(self.input_dir) if os.path.isfile(os.path.join(self.input_dir, f))]
        
        processed_count = 0
        for input_filename in input_files:
            input_path = os.path.join(self.input_dir, input_filename) 
            output_filename = f"processed_KO_list_{input_filename}"
            output_path = os.path.join(self.output_dir, output_filename)
            
            try:
                self.process_file(input_path, output_path)
                processed_count += 1
                print(f"Processed {input_filename} -> {output_filename}")
            except Exception as e:
                print(f"Error processing {input_filename}: {str(e)}")
                
        print(f"Completed processing {processed_count} files.")
        
    def process_file(self, input_file, output_file):
        """
        Process a single KofamScan output file to extract unique, sorted KO IDs.
        
        Args:
            input_file: Path to the input file
            output_file: Path to the output file
        """
        ko_ids = self.extract_ko_ids(input_file)
        
        # Sort and make unique
        unique_ko_ids = sorted(set(ko_ids))
        
        # Write to output file
        with open(output_file, 'w') as outfile:
            for ko_id in unique_ko_ids:
                outfile.write(f"{ko_id}\n")
                
    def extract_ko_ids(self, input_file):
        """
        Extract KO IDs from significant hits in a KofamScan output file.
        
        Args:
            input_file: Path to the input file
            
        Returns:
            List of KO IDs
        """
        ko_ids = []
        with open(input_file, 'r') as infile:
            for line in infile:
                # Check if line starts with asterisk (significant hit)
                if line.strip().startswith('*'):
                    # Split line by whitespace and get the KO ID (second field)
                    parts = line.strip().split()
                    if len(parts) >= 2:
                        ko_id = parts[2]  # Assuming KO ID is the third column
                        ko_ids.append(ko_id)
        
        return ko_ids

if __name__ == "__main__":
    # Can be run with optional input and output directory arguments
    input_dir = 'input_files'  # Default
    output_dir = 'significant_ko_lists'
    
    if len(sys.argv) >= 2:
        input_dir = sys.argv[1]
    if len(sys.argv) >= 3:
        output_dir = sys.argv[2]
        
    processor = KofamOutputProcessor(input_dir, output_dir)
    processor.process_directory()
