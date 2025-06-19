import os
import sys
import subprocess
import pandas as pd

def change_radii(top_file, modified_top_file):
    parmed_commands = f"""
changeRadii mbondi3
outparm {modified_top_file}
quit
"""
    with open("parmed.in", "w") as f:
        f.write(parmed_commands)

    subprocess.run(f"parmed -O -p {top_file} -i parmed.in", shell=True)
    os.remove("parmed.in")
    print(f"radii changed to mbondi3: {modified_top_file}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python change_radii_script.py <top_dir> <csv_file>")
        sys.exit(1)

    top_dir = sys.argv[1]
    csv_file = sys.argv[2]

    if not os.path.isdir(top_dir):
        print(f"Error: Directory '{top_dir}' not found.")
        sys.exit(1)

    if not os.path.isfile(csv_file):
        print(f"Error: CSV file '{csv_file}' not found.")
        sys.exit(1)

    # Load CSV
    df = pd.read_csv(csv_file)
    if "PDB codes" not in df.columns:
        print("Error: CSV file must contain a column named 'PDB codes'")
        sys.exit(1)

    # Process each PDB code
    for pdb_code in df["PDB codes"]:
        top_file = os.path.join(top_dir, f"{pdb_code}.top")
        modified_top_file = os.path.join(os.getcwd(), f"{pdb_code}.top")

        if not os.path.isfile(top_file):
            print(f"Warning: File '{top_file}' not found. Skipping.")
            continue

        change_radii(top_file, modified_top_file)
