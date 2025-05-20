import sys

def rename_numeric_ligand_to_LIG(pdb_in, pdb_out):
    """
    Replaces 3-character ligand names that are all digits with 'LIG' in a PDB file.
    Only affects fields in columns 18–20 (residue name).
    """
    with open(pdb_in, 'r') as infile, open(pdb_out, 'w') as outfile:
        for line in infile:
            if line.startswith(('HETATM', 'ATOM', 'TER')):
                resname = line[17:20]
                if resname.isdigit():
                    line = line[:17] + 'LIG' + line[20:]
            outfile.write(line)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python rename_ligand.py input.pdb output.pdb")
        sys.exit(1)
    
    input_pdb = sys.argv[1]
    output_pdb = sys.argv[2]
    rename_numeric_ligand_to_LIG(input_pdb, output_pdb)

