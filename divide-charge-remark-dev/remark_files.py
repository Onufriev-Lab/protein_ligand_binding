import os
import pandas as pd

# === CONFIGURABLE PATHS ===
csv_path = "./pdb_codes.csv"
input_dir = "/home/johann/vt_ml/pdb txt and LIGAND files"
output_dir = "./Remarked-files"

# Create output directory if it doesn't exist
os.makedirs(output_dir, exist_ok=True)

# === LOAD CSV ===
df = pd.read_csv(csv_path)
if "PDB code" not in df.columns or "ligand_charge" not in df.columns:
    raise ValueError("CSV must contain 'PDB code' and 'ligand_charge' columns")

# === PROCESS EACH ENTRY ===
for _, row in df.iterrows():
    pdb_code = str(row["PDB code"]).strip()
    charge = row["ligand_charge"]

    if pd.isna(charge):
        print(f"[{pdb_code}] Charge is empty. Skipping.")
        continue

    input_file = os.path.join(input_dir, f"{pdb_code}_LIGAND.pdb")
    output_file = os.path.join(output_dir, f"{pdb_code}_LIGAND_remarked.pdb")

    if not os.path.exists(input_file):
        print(f"[{pdb_code}] File not found: {input_file}")
        continue

    try:
        with open(input_file, "r") as f:
            lines = f.readlines()

        remark_line = f"REMARK LIGAND_NET_CHARGE {charge} END_LIGAND_INFO\n"
        new_lines = [remark_line] + lines

        with open(output_file, "w") as f:
            f.writelines(new_lines)

        print(f"[{pdb_code}] Remarked file written to {output_file}")
    except Exception as e:
        print(f"[{pdb_code}] Error processing file: {e}")
