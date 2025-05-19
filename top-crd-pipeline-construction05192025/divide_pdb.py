import os
import pandas as pd

def split_pdb_file(pdb_code, source_path, output_dir):
    with open(source_path, 'r') as f:
        lines = f.readlines()

    # Find indices of all TER lines
    ter_indices = [i for i, line in enumerate(lines) if line.startswith("TER")]

    if len(ter_indices) < 2:
        raise ValueError(f"File {source_path} does not contain at least two TER lines.")

    second_last_ter_index = ter_indices[-2]

    # Prepare the output directory
    pdb_output_dir = os.path.join(output_dir, pdb_code)
    os.makedirs(pdb_output_dir, exist_ok=True)

    # Write COMPLEX.pdb (entire source file)
    with open(os.path.join(pdb_output_dir, "COMPLEX.pdb"), 'w') as f:
        f.writelines(lines)

    # Write LIGAND.pdb: first line + lines after the second-last TER
    ligand_lines = [lines[0]] + lines[second_last_ter_index + 1:]
    with open(os.path.join(pdb_output_dir, "LIGAND.pdb"), 'w') as f:
        f.writelines(ligand_lines)

    # Write PROTEIN.pdb: up to and including the second-last TER + END
    protein_lines = lines[:second_last_ter_index + 1] + ["END\n"]
    with open(os.path.join(pdb_output_dir, "PROTEIN.pdb"), 'w') as f:
        f.writelines(protein_lines)

def process_pdbs(csv_file, output_dir="."):
    df = pd.read_csv(csv_file)

    for _, row in df.iterrows():
        pdb_code = row["PDB code"]
        source_path = f"{pdb_code}-7.0-trunc_cap_cleaned.pdb"
        if os.path.exists(source_path):
            try:
                split_pdb_file(pdb_code, source_path, output_dir)
            except Exception as e:
                print(f"Error processing {pdb_code}: {e}")
        else:
            print(f"File not found: {source_path}")

# Example usage:
# Make sure the CSV is in the working directory
process_pdbs("refined_gold_data_only.csv")
