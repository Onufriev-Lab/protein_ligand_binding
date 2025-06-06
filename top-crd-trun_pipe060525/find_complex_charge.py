import sys
import os
import re
import math

# Charged amino acid residues and their formal charges
charged_residues = {
    "ARG": +1,
    "HIP": +1,
    "HIE": 0,
    "HID": 0,
    "LYS": +1,
    "LYN": 0,
    "ASP": -1,
    "ASH": 0,
    "GLU": -1,
    "GLH": 0
}

def analyze_charges(pdb_code):
    path = f"./divided-files/{pdb_code}/"
    pdb_file = os.path.join(path, f"{pdb_code}t.pdb")
    charge_file = os.path.join(path, "charge.txt")
    top_file = os.path.join(path, f"{pdb_code}.top")
    output_file = os.path.join(path, "complex_charges.txt")

    # Check file existence
    for fname in [pdb_file, charge_file, top_file]:
        if not os.path.exists(fname):
            print(f"Error: Missing required file: {fname}")
            return

    # --- Read ligand charge ---
    with open(charge_file) as f:
        line = f.readline()
        match = re.match(r"\s*LIGAND\s*=\s*([-+]?[0-9]*\.?[0-9]+)", line)
        if not match:
            print("Error: Could not parse ligand charge from line:", line.strip())
            return
        ligand_charge = float(match.group(1))

    # --- Count protein charge ---
    charged_residues_list = []
    total_protein_charge = 0
    seen_residues = set()

    with open(pdb_file) as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                resname = line[17:20].strip()
                resid = line[22:26].strip()
                chain = line[21].strip()
                unique_id = (resname, resid, chain)
                if unique_id in seen_residues:
                    continue
                seen_residues.add(unique_id)
                if resname in charged_residues:
                    charge = charged_residues[resname]
                    total_protein_charge += charge
                    charged_residues_list.append(f"{resname} {chain}{resid}  Charge: {charge}")

    # --- Read number of water molecules from line 10 of .top file ---
    with open(top_file) as f:
        lines = f.readlines()
        try:
            water_line = lines[9]  # Line 10 = index 9
            num_waters = int(water_line.strip().split()[0])
        except (IndexError, ValueError):
            print("Error: Could not parse number of water molecules from .top file.")
            return

    # --- Calculate No ---
    No = (num_waters * 0.15) / 56

    # --- Total complex charge ---
    total_charge = total_protein_charge + ligand_charge

    # --- Calculate ion counts using corrected formulas ---
    half_q_over_no = total_charge / (2 * No)
    sqrt_term = math.sqrt(1 + (half_q_over_no) ** 2)
    cl_minus = No * sqrt_term + (total_charge / 2)
    na_plus = No * sqrt_term - (total_charge / 2)


    # --- Print results ---
    print(f"Total protein charge: {total_protein_charge}")
    print(f"Ligand charge: {ligand_charge}")
    print(f"Total complex charge: {total_charge}")
    print(f"Number of water molecules: {num_waters}")
    print(f"No (waters × 0.15 / 56): {No:.4f}")
    print(f"Chloride ions (Cl-): {cl_minus:.2f}")
    print(f"Sodium ions (Na+): {na_plus:.2f}")

    # --- Write output file ---
    with open(output_file, "w") as out:
        out.write("Charged residues:\n")
        for entry in charged_residues_list:
            out.write(entry + "\n")
        out.write(f"\nLigand charge from charge.txt: {ligand_charge}\n")
        out.write(f"Total protein charge: {total_protein_charge}\n")
        out.write(f"Total complex charge: {total_charge}\n")
        out.write(f"Number of water molecules from {pdb_code}.top: {num_waters}\n")
        out.write(f"No (waters × 0.15 / 56): {No:.4f}\n")
        out.write(f"Chloride ions (Cl-): {cl_minus:.2f}\n")
        out.write(f"Sodium ions (Na+): {na_plus:.2f}\n")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python charge_analysis.py <pdb_code>")
    else:
        analyze_charges(sys.argv[1])
