import subprocess
import os
import pandas as pd

# === CONFIGURABLE VARIABLES ===
csv_file = "refined_gold_data_only_full.csv"
pdb_column = "PDB code"
charge_column = "lig_charge"

def run_antechamber_and_parmchk2(pdb_code, charge):
    pdb_dir = f"./divided-files/{pdb_code}"
    input_pdb = os.path.join(pdb_dir, "LIGAND.pdb")
    output_mol2 = os.path.join(pdb_dir, "LIGAND.mol2")
    output_frcmod = os.path.join(pdb_dir, "LIGAND.frcmod")

    # Run antechamber
    cmd_antechamber = [
        "antechamber",
        "-i", input_pdb,
        "-fi", "pdb",
        "-o", output_mol2,
        "-fo", "mol2",
        "-c", "bcc",
        "-s", "2",
        "-nc", str(charge),
        "-pl", "15"
    ]

    try:
        print(f"[{pdb_code}] Running antechamber...")
        subprocess.run(cmd_antechamber, check=True)
        print(f"[{pdb_code}] Antechamber completed.")
    except subprocess.CalledProcessError as e:
        print(f"[{pdb_code}] Antechamber failed: {e}")
        return

    if not os.path.isfile(output_mol2):
        print(f"[{pdb_code}] Error: {output_mol2} was not created.")
        return

    # Run parmchk2
    cmd_parmchk2 = [
        "parmchk2",
        "-i", output_mol2,
        "-f", "mol2",
        "-o", output_frcmod
    ]

    try:
        print(f"[{pdb_code}] Running parmchk2...")
        subprocess.run(cmd_parmchk2, check=True)
        print(f"[{pdb_code}] parmchk2 completed.")
    except subprocess.CalledProcessError as e:
        print(f"[{pdb_code}] parmchk2 failed: {e}")

if __name__ == "__main__":
    try:
        df = pd.read_csv(csv_file)
    except FileNotFoundError:
        print(f"Error: CSV file '{csv_file}' not found.")
        exit(1)

    for idx, row in df.iterrows():
        pdb_code = str(row[pdb_column]).strip()
        charge = row.get(charge_column)

        try:
            charge = float(charge)
        except (ValueError, TypeError):
            print(f"[{pdb_code}] Invalid charge '{charge}'; skipping.")
            continue

        run_antechamber_and_parmchk2(pdb_code, charge)
