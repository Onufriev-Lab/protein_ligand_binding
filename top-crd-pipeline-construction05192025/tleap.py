import os
import sys
import subprocess

def create_tleap_input(pdb_code):
    directory = os.path.join('.', pdb_code)
    os.makedirs(directory, exist_ok=True)

    tleap_lines = [
        "source leaprc.protein.ff19SB",
        "source leaprc.gaff",
        "set default PBradii mbondi3",
        f"loadamberparams ./{pdb_code}/LIGAND.frcmod",
        f"LIG = loadmol2 ./{pdb_code}/LIGAND.mol2",
        f"PROT = loadpdb ./{pdb_code}/PROTEIN.pdb",
        "COMPLEX = combine {PROT LIG}",
        f"saveamberparm COMPLEX ./{pdb_code}/{pdb_code}.top ./{pdb_code}/{pdb_code}.crd",
        f"savepdb COMPLEX ./{pdb_code}/{pdb_code}.pdb",
        "quit"
    ]

    tleap_path = os.path.join(directory, "tleap.in")
    with open(tleap_path, 'w') as f:
        f.write('\n'.join(tleap_lines))

    print(f"tleap.in file created at: {tleap_path}")

    # Run tleap command
    try:
        subprocess.run(['tleap', '-f', tleap_path], check=True)
        print("tleap completed successfully.")
    except subprocess.CalledProcessError as e:
        print(f"Error running tleap: {e}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python create_tleap_input.py <pdb_code>")
        sys.exit(1)
    
    pdb_code = sys.argv[1]
    create_tleap_input(pdb_code)
