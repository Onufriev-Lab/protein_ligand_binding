import os
import re
import glob
import argparse

def scan_pdb_directories(target_dir):
    pdb_dir_pattern = re.compile(r'^[a-zA-Z0-9]{4}$')  # 4-character alphanumeric PDB codes

    total_dirs = 0
    with_file = 0
    without_file = 0

    for entry in os.listdir(target_dir):
        full_path = os.path.join(target_dir, entry)
        if os.path.isdir(full_path) and pdb_dir_pattern.fullmatch(entry):
            total_dirs += 1
            md1_files = glob.glob(os.path.join(full_path, '*_md1.out'))
            if md1_files:
                with_file += 1
            else:
                without_file += 1

    print(f"Total directories scanned: {total_dirs}")
    print(f"Directories with *_md1.out: {with_file}")
    print(f"Directories without *_md1.out: {without_file}")

if __name__ == "__main__": #usage:  python count_MD_progress.py dir_with_MD_directories
    parser = argparse.ArgumentParser(description="Scan PDB directories for *_md1.out files.")
    parser.add_argument("target_directory", help="Path to the directory containing 4-letter PDB subdirectories.")
    args = parser.parse_args()

    scan_pdb_directories(args.target_directory)
