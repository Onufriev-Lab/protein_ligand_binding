import os
import csv

# Set the path to your directory
directory = "/home/johann/vt_ml/pdb txt and LIGAND files"
output_txt = "pdb_codes_list.txt"
output_csv = "pdb_codes.csv"

# Collect pdb_codes
pdb_codes = []
for filename in os.listdir(directory):
    if filename.endswith("_LIGAND.pdb"):
        pdb_code = filename.replace("_LIGAND.pdb", "")
        pdb_codes.append(pdb_code)

# Sort for consistency
pdb_codes.sort()

# Write Python list to txt
with open(output_txt, "w") as f:
    f.write("pdb_codes = [\n")
    for code in pdb_codes:
        f.write(f"    '{code}',\n")
    f.write("]\n")
    f.write(f"\n# Total files: {len(pdb_codes)}\n")

# Write CSV
with open(output_csv, "w", newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(['PDB code'])  # header
    for code in pdb_codes:
        writer.writerow([code])

print(f"Found {len(pdb_codes)} matching files.")
print(f"Output written to:\n  - {output_txt}\n  - {output_csv}")
