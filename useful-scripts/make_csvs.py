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

def split_into_groups(pdb_codes, n):
    total = len(pdb_codes)
    size = total // n
    groups = []

    for i in range(n):
        start = i * size
        end = (i + 1) * size if i != n - 1 else total  # last group takes remainder
        groups.append(pdb_codes[start:end])

    return groups

def save_groups_to_csv(groups):
    for i, group in enumerate(groups):
        df = pd.DataFrame(group, columns=["PDB codes"])
        df.to_csv(f"group_{i}.csv", index=False)

def main():
    parser = argparse.ArgumentParser(description="Split PDB codes into n CSV groups")
    parser.add_argument("directory", help="Directory to scan for .top and .crd files")
    parser.add_argument("n", type=int, help="Number of groups to create")
    args = parser.parse_args()

    pdb_codes = get_unique_pdb_codes(args.directory)
    groups = split_into_groups(pdb_codes, args.n)
    save_groups_to_csv(groups)

if __name__ == "__main__":
    main()
