import parmed as pmd
import pandas as pd
import os
import sys

def change_radii_to_mbondi3(top_file, output_file):
    parm = pmd.load_file(top_file)
    parm.radii = 'mbondi3'
    parm.save(output_file, overwrite=True)

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Usage: python create_mbondii3_top.py <directory> <csv_file>")
        sys.exit(1)

    top_dir = sys.argv[1]
    csv_file = sys.argv[2]

    # Load the list of PDB codes
    df = pd.read_csv(csv_file)
    if 'PDB codes' not in df.columns:
        print("CSV file must contain a column named 'PDB codes'")
        sys.exit(1)

    pdb_codes = df['PDB codes'].dropna().astype(str)

    for code in pdb_codes:
        input_top = os.path.join(top_dir, f"{code}.top")
        output_top = f"{code}.top"

        if not os.path.exists(input_top):
            print(f"WARNING: {input_top} not found. Skipping.")
            continue

        print(f"Updating radii to mbondi3 for {code}...")
        change_radii_to_mbondi3(input_top, output_top)
        print(f"Saved to {output_top}")
