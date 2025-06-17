import sys
from collections import defaultdict

CAP_ATOMS = {
    'NME': {'N', 'H'},
    'ACE': {'C', 'O'}  # This varies, so be cautious
}

# Full residue signatures using common heavy atoms
RESIDUE_SIGNATURES = {
    'ALA': {'CA', 'CB'},
    'ARG': {'CA', 'CB', 'CG', 'CD', 'NE', 'CZ', 'NH1', 'NH2'},
    'ASN': {'CA', 'CB', 'CG', 'OD1', 'ND2'},
    'ASP': {'CA', 'CB', 'CG', 'OD1', 'OD2'},
    'CYS': {'CA', 'CB', 'SG'},
    'GLN': {'CA', 'CB', 'CG', 'CD', 'OE1', 'NE2'},
    'GLU': {'CA', 'CB', 'CG', 'CD', 'OE1', 'OE2'},
    'GLY': {'CA'},  # Only backbone atom
    'HIS': {'CA', 'CB', 'CG', 'ND1', 'CD2', 'CE1', 'NE2'},
    'ILE': {'CA', 'CB', 'CG1', 'CG2', 'CD1'},
    'LEU': {'CA', 'CB', 'CG', 'CD1', 'CD2'},
    'LYS': {'CA', 'CB', 'CG', 'CD', 'CE', 'NZ'},
    'MET': {'CA', 'CB', 'CG', 'SD', 'CE'},
    'PHE': {'CA', 'CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'},
    'PRO': {'CA', 'CB', 'CG', 'CD'},
    'SER': {'CA', 'CB', 'OG'},
    'THR': {'CA', 'CB', 'OG1', 'CG2'},
    'TRP': {'CA', 'CB', 'CG', 'CD1', 'CD2', 'NE1', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2'},
    'TYR': {'CA', 'CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ', 'OH'},
    'VAL': {'CA', 'CB', 'CG1', 'CG2'}
}


def guess_residue_name(atom_names):
    for resname, signature in RESIDUE_SIGNATURES.items():
        if signature == atom_names:
            return resname
    else:
        print("Residue name not found")
    return None

def read_pdb(pdb_lines):
    residues = defaultdict(list)
    for line in pdb_lines:
        if line.startswith(('ATOM', 'HETATM')):
            res_id = (line[17:20].strip(), int(line[22:26]), line[21])  # (resname, resnum, chain)
            residues[res_id].append(line)
    return residues

def fix_mislabeled_caps(residues):
    fixed_lines = []
    res_keys = sorted(residues.keys(), key=lambda x: (x[2], x[1]))  # Sort by chain and resnum

    for i, key in enumerate(res_keys):
        resname, resnum, chain = key
        lines = residues[key]
        atom_names = {
            line[12:16].strip()
            for line in lines
            if line[76:78].strip() != 'H' and line[12:16].strip() not in {'N', 'C', 'O'}
            }

        print ("atom_names = ", atom_names)


        if resname in CAP_ATOMS:
            is_cap = atom_names.issubset(CAP_ATOMS[resname])
            if not is_cap:
                new_name = guess_residue_name(atom_names)
                if new_name:
                    print(f"Fixing mislabeled residue at {resname} {resnum}{chain} -> {new_name}")
                    new_lines = []
                    for line in lines:
                        new_line = line[:17] + new_name.rjust(3) + line[20:]
                        new_lines.append(new_line)
                    fixed_lines.extend(new_lines)
                    continue  # Skip original lines
        fixed_lines.extend(lines)

    return fixed_lines

def main(pdb_file):
    with open(pdb_file) as f:
        lines = f.readlines()

    residues = read_pdb(lines)
    fixed_lines = fix_mislabeled_caps(residues)

    with open(pdb_file.replace(".pdb", "_fixed.pdb"), "w") as f:
        for line in fixed_lines:
            f.write(line)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python fix_mislabeled_AA.py <input.pdb>")
    else:
        main(sys.argv[1])
