import os
import sys
import shutil
import pandas as pd

def create_bash_script(pdb_code, gpu_id, subdir_path):
    script_content = f"""#!/bin/bash

export CUDA_VISIBLE_DEVICES={gpu_id}

MODEL_NAME="{pdb_code}"
#min
pmemd -O -i min-imp.in -o "$MODEL_NAME"_min.out -c "$MODEL_NAME".crd -p "$MODEL_NAME".top -r "$MODEL_NAME"_min.rst -ref "$MODEL_NAME".crd
#heat1
pmemd.cuda -O -i heat-imp.in -o "$MODEL_NAME"_heat1.out -c "$MODEL_NAME"_min.rst -p "$MODEL_NAME".top -r "$MODEL_NAME"_heat1.rst -x "$MODEL_NAME"_heat1.mdcrd -ref "$MODEL_NAME"_min.rst 
#equilibrate 1
pmemd.cuda -O -i density-imp.in -o "$MODEL_NAME"_equil1.out -c "$MODEL_NAME"_heat1.rst -p "$MODEL_NAME".top -r "$MODEL_NAME"_equil1.rst -x "$MODEL_NAME"_equil1.mdcrd -ref "$MODEL_NAME"_heat1.rst
#equilibrate 2
pmemd.cuda -O -i equil-imp.in -o "$MODEL_NAME"_equil2.out -c "$MODEL_NAME"_equil1.rst -p "$MODEL_NAME".top -r "$MODEL_NAME"_equil2.rst -x "$MODEL_NAME"_equil2.mdcrd -ref "$MODEL_NAME"_equil1.rst

#md
pmemd.cuda -O -i md-imp.in -o "$MODEL_NAME"_md1.out -c "$MODEL_NAME"_equil2.rst -p "$MODEL_NAME".top -r "$MODEL_NAME"_md1.rst -x "$MODEL_NAME"_md1.mdcrd
"""
    script_path = os.path.join(subdir_path, f"{pdb_code}_run.sh")
    with open(script_path, 'w') as f:
        f.write(script_content)
    os.chmod(script_path, 0o755)  # Make it executable

def main(input_dir, csv_file, gpu_id):
    # Read CSV
    df = pd.read_csv(csv_file)
    if "PDB codes" not in df.columns:
        print("Error: CSV file must contain a column named 'PDB codes'")
        sys.exit(1)

    for pdb_code in df["PDB codes"].dropna().astype(str):
        subdir_path = os.path.join(input_dir, pdb_code)
        if not os.path.isdir(subdir_path):
            print(f"Warning: Subdirectory not found for {pdb_code}, skipping.")
            continue

        src_top = os.path.join(subdir_path, f"{pdb_code}_hmr.top")
        dst_top = os.path.join(subdir_path, f"{pdb_code}.top")
        if not os.path.isfile(src_top):
            print(f"Warning: .top file not found for {pdb_code}, skipping.")
            continue

        shutil.copyfile(src_top, dst_top)
        create_bash_script(pdb_code, gpu_id, subdir_path)
        print(f"Processed {pdb_code}")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script.py <input_dir> <csv_file> <gpu_id>")
        sys.exit(1)
    input_dir = sys.argv[1]
    csv_file = sys.argv[2]
    gpu_id = sys.argv[3]
    main(input_dir, csv_file, gpu_id)
