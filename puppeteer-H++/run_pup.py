import pandas as pd
import subprocess

# Path to your CSV file
csv_file = "./pdb_codes.csv"
pdb_column = "PDB code"

# Path to your JavaScript file
js_script = "pupH++1.js"

# Read CSV
df = pd.read_csv(csv_file)

# Drop missing PDB codes and clean
pdb_codes = df[pdb_column].dropna().astype(str).str.strip()

# Loop and run the JS script
for pdb_code in pdb_codes:
    print(f"Running: node {js_script} {pdb_code}")
    try:
        result = subprocess.run(
            ["node", js_script, pdb_code],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        print(f"[{pdb_code}] Success\n{result.stdout}")
    except subprocess.CalledProcessError as e:
        print(f"[{pdb_code}] Failed\n{e.stderr}")
