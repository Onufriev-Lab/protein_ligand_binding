import os
import argparse
import pandas as pd

def get_unique_pdb_codes(directory):
    filenames = os.listdir(directory)
    pdb_codes = set()

    for fname in filenames:
        if fname.endswith('.top') or fname.endswith('.crd'):
            pdb_code = os.path.splitext(fname)[0]
            pdb_codes.add(pdb_code)

    return sorted(pdb_codes)

def save_to_csv(pdb_codes, output_file="all_pdb_codes.csv"):
    df = pd.DataFrame(pdb_codes, columns=["PDB codes"])
    df.to_csv(output_file, index=False)

def main():
    parser = argparse.ArgumentParser(description="Save unique PDB codes to a CSV file")
    parser.add_argument("directory", help="Directory to scan for .top and .crd files")
    args = parser.parse_args()

    pdb_codes = get_unique_pdb_codes(args.directory)
    save_to_csv(pdb_codes)

if __name__ == "__main__":
    main()
