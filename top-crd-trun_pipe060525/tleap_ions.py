import os
import sys
import subprocess

def create_tleap_input(pdb_code, dist, Cl, Na):


    tleap_lines = [
        "source leaprc.protein.ff19SB",
        "source leaprc.gaff",
        f"loadamberparams ./divided-files/{pdb_code}/LIGAND.frcmod",
        f"LIG = loadmol2 ./divided-files/{pdb_code}/LIGAND.mol2",
        f"PROT = loadpdb ./divided-files/{pdb_code}/PROTEIN.pdb",
        "COMPLEX = combine {PROT LIG}",
        "source leaprc.water.opc",
        f"solvateOct COMPLEX OPCBOX {dist}",
        f"addions COMPLEX Na+ {Na}",
        f"addions COMPLEX Cl- {Cl}",
        f"saveamberparm COMPLEX ./divided-files/{pdb_code}/{pdb_code}.top ./divided-files/{pdb_code}/{pdb_code}.crd",
        f"savepdb COMPLEX ./divided-files/{pdb_code}/{pdb_code}ti.pdb",
        "quit"
    ]

    tleap_path = f"./divided-files/{pdb_code}/tleap.in"
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
    if len(sys.argv) != 5:
        print("Usage: python create_tleap_input.py <pdb_code>")
        sys.exit(1)
    
    pdb_code = sys.argv[1]
    dist = sys.argv[2]
    Cl = sys.argv[3]
    Na = sys.argv[4]
    create_tleap_input(pdb_code, dist, Cl, Na)
