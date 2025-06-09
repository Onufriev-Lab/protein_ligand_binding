import os
import csv
import pandas as pd

# Set the path to your directory
directory = "/home/johann/vt_ml/pdb txt and LIGAND files"
output_txt = "pdb_codes_list.txt"
output_csv = "pdb_codes.csv"

# Collect pdb_codes and ligand names
records = []

for filename in os.listdir(directory):
    if filename.endswith("_LIGAND.pdb"):
        pdb_code = filename.replace("_LIGAND.pdb", "")
        ligand_name = "Not found"

        file_path = os.path.join(directory, filename)
        try:
            with open(file_path, 'r') as file:
                for line in file:
                    if line.startswith("LIGAND"):
                        ligand_name = line[17:20].strip()
                        break
        except Exception as e:
            ligand_name = f"Error: {e}"

        records.append((pdb_code, ligand_name))

# Sort records by pdb_code
records.sort(key=lambda x: x[0])

# Write Python list to txt
with open(output_txt, "w") as f:
    f.write("pdb_codes = [\n")
    for code, _ in records:
        f.write(f"    '{code}',\n")
    f.write("]\n")
    f.write(f"\n# Total files: {len(records)}\n")

# Write CSV
with open(output_csv, "w", newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(['PDB code', 'Ligand Name'])  # header
    for code, ligand in records:
        writer.writerow([code, ligand])

