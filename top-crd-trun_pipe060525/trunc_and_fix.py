import subprocess
import pandas as pd
import os
import shutil

# === CONFIGURABLE VARIABLES ===
csv_path = "refined_gold_data_only_full.csv"
pdb_column = "PDB code"
ligand_column = "Ligand Name"
radius_column = "trunc_rad"
output_dir = "trunc-cap-cleaned-files"

def run_pipeline(pdb_code, ligand_name, radius, output_dir):
    # Step 1: Run the bash script
    bash_script = './TruncAndCapLinux2.sh'
    try:
        subprocess.run([bash_script, pdb_code, ligand_name, str(radius)], check=True)
    except subprocess.CalledProcessError as e:
        print(f"[{pdb_code}] Error running bash script: {e}")
        return

    # Step 2: Rename 'trunc_cap.pdb'
    original_output = 'trunc_cap.pdb'
    renamed_output = f"{pdb_code}-{radius}-trunc_cap.pdb"

    if not os.path.exists(original_output):
        print(f"[{pdb_code}] Error: Expected output file '{original_output}' not found.")
        return

    try:
        shutil.move(original_output, renamed_output)
        print(f"[{pdb_code}] Renamed '{original_output}' to '{renamed_output}'")
    except Exception as e:
        print(f"[{pdb_code}] Error renaming file: {e}")
        return

    # Step 3: Run the Python post-processing script
    try:
        subprocess.run(['python3', 'fix_trunc_capped2.py', renamed_output], check=True)
    except subprocess.CalledProcessError as e:
        print(f"[{pdb_code}] Error running Python script: {e}")
        return

    # Step 4: Move cleaned file into subdirectory
    cleaned_file = renamed_output.replace('.pdb', '_cleaned.pdb')
    final_destination = os.path.join(output_dir, cleaned_file)
    if os.path.exists(cleaned_file):
        try:
            shutil.move(cleaned_file, final_destination)
            print(f"[{pdb_code}] Moved cleaned file to '{final_destination}'")
        except Exception as e:
            print(f"[{pdb_code}] Error moving cleaned file: {e}")
    else:
        print(f"[{pdb_code}] Cleaned file not found: {cleaned_file}")

if __name__ == "__main__":
    # Create the output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    # Load input CSV
    try:
        df = pd.read_csv(csv_path)
    except FileNotFoundError:
        print(f"Error: CSV file '{csv_path}' not found.")
        exit(1)

    # Iterate through the DataFrame
    for idx, row in df.iterrows():
        pdb_code = str(row[pdb_column]).strip()
        ligand_name = str(row[ligand_column]).strip()
        try:
            radius = float(row[radius_column])
        except (ValueError, KeyError):
            print(f"[Row {idx}] Invalid or missing radius; skipping.")
            continue

        print(f"\nProcessing: PDB={pdb_code}, Ligand={ligand_name}, Radius={radius}")
        run_pipeline(pdb_code, ligand_name, radius, output_dir)
