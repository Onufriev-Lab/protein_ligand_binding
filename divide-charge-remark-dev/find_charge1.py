import os
import subprocess
import pandas as pd

# === CONFIGURABLE VARIABLES ===
csv_file = "pdb_codes.csv"
pdb_column = "PDB code"
ligand_pdb_name = "_ligand_only.pdb"
ligand_mol2_name = "_ligand_only.mol2"
charge_file_name = "_charge.txt"
target_pH = 7.0

# === FUNCTION TO PARSE CHARGES FROM MOL2 ===
def extract_total_charge(mol2_path):
    total_charge = 0.0
    try:
        with open(mol2_path, 'r') as f:
            for line in f:
                if line.lower().startswith("charge"):
                    parts = line.split()
                    if len(parts) >= 2:
                        try:
                            total_charge += float(parts[1])
                        except ValueError:
                            continue
    except FileNotFoundError:
        print(f"[ERROR] .mol2 file not found: {mol2_path}")
        return None
    return total_charge

# === MAIN PROCESSING LOOP ===
df = pd.read_csv(csv_file)

# Add a new column for ligand charges if not already present
if "ligand_charge" not in df.columns:
    df["ligand_charge"] = None

for idx, row in df.iterrows():
    pdb_code = str(row[pdb_column]).strip()
    pdb_dir = f'./divided-files/'
    lig_file = pdb_code + ligand_pdb_name
    pdb_path = os.path.join(pdb_dir, lig_file)
    lig_mol2_file = pdb_code + ligand_mol2_name
    mol2_path = os.path.join(pdb_dir, lig_mol2_file)
    temp_mol2_path = os.path.join(pdb_dir, f"temp{ligand_mol2_name}")
    charge_path = os.path.join(pdb_dir, charge_file_name)

    if not os.path.exists(pdb_path):
        print(f"[{pdb_code}] LIGAND.pdb not found. Skipping.")
        print(f"pdb_path = {pdb_path}")
        continue

    # Run obabel to convert to mol2 with specified pH
    try:
        subprocess.run([
            "obabel", pdb_path,
            "-O", temp_mol2_path,
            "-d"
        ], check=True)
        print(f"[{pdb_code}] Successfully converted to mol2 w/o H.")
    except subprocess.CalledProcessError as e:
        print(f"[{pdb_code}] Error running obabel: {e}")
        continue

    try:
        subprocess.run([
            "obabel", temp_mol2_path,
            "-O", mol2_path,
            f"-p", str(target_pH)
        ], check=True)
        print(f"[{pdb_code}] Successfully converted to mol2 w/ H.")
    except subprocess.CalledProcessError as e:
        print(f"[{pdb_code}] Error running obabel: {e}")
        continue

    # Extract total charge
    total_charge = extract_total_charge(mol2_path)
    if total_charge is None:
        print(f"[{pdb_code}] Failed to extract charge.")
        continue

    # Write charge to file
    try:
        with open(charge_path, 'w') as f:
            f.write(f"LIGAND = {total_charge:.4f}\n")
        print(f"[{pdb_code}] Charge written to {charge_path}")
    except Exception as e:
        print(f"[{pdb_code}] Failed to write charge.txt: {e}")
        continue

    # Save charge to DataFrame
    df.at[idx, "ligand_charge"] = total_charge

# Save the updated CSV
df.to_csv(csv_file, index=False)
print(f"Updated CSV with ligand charges saved to {csv_file}")
