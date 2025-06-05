import sys

def compare_residue_names(pdb_file1, pdb_file2):
    with open(pdb_file1, 'r') as f1, open(pdb_file2, 'r') as f2:
        lines1 = f1.readlines()
        lines2 = f2.readlines()

    for i, (line1, line2) in enumerate(zip(lines1, lines2), start=1):
        if line1.startswith(('ATOM', 'HETATM')) and line2.startswith(('ATOM', 'HETATM')):
            resname1 = line1[17:20].strip()
            resname2 = line2[17:20].strip()

            # Stop comparing once both files have reached water residues
            if resname1 in ('HOH', 'WAT') and resname2 in ('HOH', 'WAT'):
                break

            if resname1 != resname2:
                print(f"Line {i}: ResName mismatch")
                print(f"File 1: {line1.rstrip()}")
                print(f"File 2: {line2.rstrip()}\n")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python compare_residue_names.py file1.pdb file2.pdb")
        sys.exit(1)
    compare_residue_names(sys.argv[1], sys.argv[2])