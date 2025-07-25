"""
Debug script to examine sample ID format
"""

from data_reader import GTExDataReader

def debug_sample_ids():
    reader = GTExDataReader()
    sample_attrs = reader.read_sample_attributes()
    
    print("First 10 sample IDs:")
    for i, sample_id in enumerate(sample_attrs.index[:10]):
        print(f"  {i+1}: {sample_id}")
    
    # Test regex extraction
    import re
    pattern = r'(GTEX-[^-]+)'
    
    print(f"\nTesting regex pattern: {pattern}")
    for sample_id in sample_attrs.index[:5]:
        match = re.search(pattern, sample_id)
        if match:
            print(f"  {sample_id} -> {match.group(1)}")
        else:
            print(f"  {sample_id} -> NO MATCH")

if __name__ == "__main__":
    debug_sample_ids()
