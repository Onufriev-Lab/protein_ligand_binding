import subprocess
from pathlib import Path

# Mapping from pdb_code to Z
Z_LOOKUP = {
    '4mre': 151,
    '6r9u': 165,
    '6oqc': 153
}

# Root directory
base_path = Path.home() / 'full-imp-pipeline-work/gold-training-062725-heat-fix'

# Traverse to each mbondi directory
for x_dir in base_path.iterdir():
    if not x_dir.is_dir():
        continue

    for y_dir in x_dir.iterdir():
        mbondi_dir = y_dir / 'mbondi'
        if not mbondi_dir.is_dir():
            continue

        # Find the .top file to get pdb_code
        top_files = list(mbondi_dir.glob('*.top'))
        if not top_files:
            print(f"WARNING: No .top file found in {mbondi_dir}")
            continue

        pdb_code = top_files[0].stem[:4].lower()
        if pdb_code not in Z_LOOKUP:
            print(f"WARNING: Unknown PDB code '{pdb_code}' in {mbondi_dir}")
            continue

        Z = Z_LOOKUP[pdb_code]
        W = Z - 1

        # Commands to run
        commands = [
            ['ante-MMPBSA.py', '-p', f'{pdb_code}.top', '-c', 'complex.prmtop', '-s', f':WAT'],
            ['ante-MMPBSA.py', '-p', 'complex.prmtop', '-c', 'protein.prmtop', '-s', f':{Z}'],
            ['ante-MMPBSA.py', '-p', 'complex.prmtop', '-c', 'ligand.prmtop', '-s', f':1-{W}'],
        ]

        # Run each command in mbondi_dir
        for cmd in commands:
            print(f"Running in {mbondi_dir}: {' '.join(cmd)}")
            try:
                subprocess.run(cmd, cwd=mbondi_dir, check=True)
            except subprocess.CalledProcessError as e:
                print(f"ERROR running command: {' '.join(cmd)}\n{e}")
