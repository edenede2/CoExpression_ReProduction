#!/usr/bin/env python3
"""
Fix CSV file with unquoted commas in text fields.
Reads a malformed CSV and writes a properly quoted version.
"""

import csv
import sys

def fix_csv_file(input_path, output_path):
    """Read CSV with potential issues and write properly quoted version."""
    
    # Read the file line by line to handle malformed rows
    with open(input_path, 'r', encoding='utf-8') as infile:
        lines = infile.readlines()
    
    # Parse header
    header_line = lines[0].strip()
    # Count expected fields from header
    # We'll use csv module with proper quoting
    
    with open(output_path, 'w', encoding='utf-8', newline='') as outfile:
        writer = csv.writer(outfile, quoting=csv.QUOTE_NONNUMERIC)
        
        # Write header
        header_reader = csv.reader([header_line])
        header = next(header_reader)
        writer.writerow(header)
        
        expected_fields = len(header)
        print(f"Expected {expected_fields} fields per row")
        
        # Process data rows
        for line_num, line in enumerate(lines[1:], start=2):
            line = line.strip()
            if not line:
                continue
                
            # Split by comma but be smart about it
            # The issue is in the Description field (column 6)
            # Format: Cluster,Cluster.ID,category,subcategory,ID,Description,GeneRatio,...
            
            parts = line.split(',')
            
            if len(parts) == expected_fields:
                # Line is OK, just write it
                reader = csv.reader([line])
                row = next(reader)
                writer.writerow(row)
            elif len(parts) > expected_fields:
                # Too many commas - likely in Description field (index 5)
                # Reconstruct the row by merging Description parts
                
                # First 5 fields are safe: Cluster,Cluster.ID,category,subcategory,ID
                fixed_row = parts[0:5]
                
                # Calculate how many extra commas we have
                extra_commas = len(parts) - expected_fields
                
                # Description field plus extra parts (should be 1 + extra_commas parts)
                description_parts = parts[5:5+extra_commas+1]
                description = ','.join(description_parts)
                fixed_row.append(description)
                
                # Add remaining fields
                remaining = parts[5+extra_commas+1:]
                fixed_row.extend(remaining)
                
                if len(fixed_row) == expected_fields:
                    writer.writerow(fixed_row)
                    print(f"Fixed line {line_num}: merged Description field")
                else:
                    print(f"WARNING: Line {line_num} still has {len(fixed_row)} fields after fix")
                    writer.writerow(fixed_row)
            else:
                print(f"WARNING: Line {line_num} has only {len(parts)} fields (expected {expected_fields})")
                writer.writerow(parts)
    
    print(f"\nFixed CSV written to: {output_path}")

if __name__ == "__main__":
    input_file = "/Users/edeneldar/Library/CloudStorage/GoogleDrive-edenede2@gmail.com/My Drive/ROSMAP_Data/RNAseq_RISK_BrainRegions_STasaki/outputs/kegg_rosmap_constBeta_CT2_TS3.csv"
    output_file = "/Users/edeneldar/Library/CloudStorage/GoogleDrive-edenede2@gmail.com/My Drive/ROSMAP_Data/RNAseq_RISK_BrainRegions_STasaki/outputs/kegg_rosmap_constBeta_CT2_TS3_fixed.csv"
    
    fix_csv_file(input_file, output_file)
    print("\nDone! Use the '_fixed.csv' file instead of the original.")
