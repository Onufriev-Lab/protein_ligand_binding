import subprocess
import sys
import os

def run_antechamber_and_parmchk2(pdb_code, charge):
    pdb_dir = f"./{pdb_code}"
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
        print(f"Running antechamber for {pdb_code}...")
        subprocess.run(cmd_antechamber, check=True)
        print("Antechamber completed.")
    except subprocess.CalledProcessError as e:
        print(f"Antechamber failed for {pdb_code}: {e}")
        return

    # Check that the .mol2 file exists before running parmchk2
    if not os.path.isfile(output_mol2):
        print(f"Error: {output_mol2} was not created.")
        return

    # Run parmchk2
    cmd_parmchk2 = [
        "parmchk2",
        "-i", output_mol2,
        "-f", "mol2",
        "-o", output_frcmod
    ]

    try:
        print("Running parmchk2...")
        subprocess.run(cmd_parmchk2, check=True)
        print("parmchk2 completed.")
    except subprocess.CalledProcessError as e:
        print(f"parmchk2 failed for {pdb_code}: {e}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python run_antechamber.py <pdb_code> <charge>")
        sys.exit(1)

    pdb_code = sys.argv[1]
    charge = sys.argv[2]

    run_antechamber_and_parmchk2(pdb_code, charge)
