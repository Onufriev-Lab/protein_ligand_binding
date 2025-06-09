import os
import pandas as pd

# File and directory locations
csv_path = "pdb_codes.csv"
pdb_input_dir = "/home/johann/vt_ml/pdb txt and LIGAND files"
output_dir = "./divided-files"

# Ensure output directory exists
os.makedirs(output_dir, exist_ok=True)

# Load the PDB codes and ligand names
df = pd.read_csv(csv_path)

# Process each PDB file
for _, row in df.iterrows():
    pdb_code = row["PDB code"]
    ligand_name = row["Ligand Name"]

    pdb_filename = os.path.join(pdb_input_dir, f"{pdb_code}.pdb")
    output_filename = os.path.join(output_dir, f"{pdb_code}_ligand_only.pdb")

    if ligand_name == "Not found":
        print(f"{pdb_code}: Ligand name not found — skipping.")
        continue

    if not os.path.isfile(pdb_filename):
        print(f"{pdb_code}: File not found — skipping.")
        continue

    with open(pdb_filename, "r") as f:
        lines = f.readlines()

    ligand_lines = [
        line for line in lines
        if line.startswith(("HETATM")) and line[17:20].strip() == ligand_name
    ]

    if ligand_lines:
        with open(output_filename, "w") as out_f:
            out_f.writelines(ligand_lines)
        print(f"{pdb_code}: Ligand written to {output_filename}")
    else:
        print(f"{pdb_code}: Ligand {ligand_name} not found in file.")
