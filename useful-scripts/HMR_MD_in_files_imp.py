import subprocess
import glob
import re
import shutil
import os
import sys
import pandas as pd

# -----------------------
# Input templates
# -----------------------
input_files = {
    'min-imp.in': """minimise complex
&cntrl
  imin=1,        ! do minimisation
  maxcyc=1000,   ! total cycles
  ncyc=500,      ! switch from steepest descent to conjugate‑gradient
  cut=999.0,     ! nonbonded cutoff (large for GB)

  igb=8,         ! GB model
  ntb=0,         ! no periodic box
  ntc=1,         ! no SHAKE (all bonds flexible)
  ntf=1,         ! force evaluation without SHAKE
  ntpr=100,      ! print energy every 100 steps
  ntr=1,         ! positional restraints on
  restraintmask='(:1-{last_residue}@CA,N,C|:{ligand_residue})',  ! same mask as explicit
  restraint_wt=10.0,                        ! same force constant
/
""",
    'heat-imp.in': """heat complex
&cntrl
  imin=0,        ! MD heating
  irest=0,       ! start from scratch
  ntx=1,         ! coordinates only, no velocities
  nstlim=14286,  ! total steps (≈50 ps @3.5 fs)
  dt=0.0035,     ! timestep (ps)
  
  ntc=2,         ! SHAKE on bonds involving H
  ntf=2,         ! no forces on bonds with SHAKE
  
  ig=-1,         ! random seed
  ntt=3,         ! Langevin thermostat
  gamma_ln=2.0,  ! collision freq (ps**-1)
  tempi=0.0,     ! start temp (K)
  temp0=300.0,   ! target temp (K)
  
  ntpr=286,      ! print to mdout every 286 steps
  ntwx=286,      ! write coords every 286 steps
  
  ntb=0,         ! no periodic box
  igb=8,         ! GB model
  saltcon=0.145, ! 0.145 M salt
  
  cut=999.0,     ! nonbonded cutoff for GB
  
  ntr=1,         ! positional restraints on
  restraintmask='(:1-{last_residue}@CA,N,C|:{ligand_residue})',
  restraint_wt=2.0,
  
  nmropt=1,      ! enable NMR-style &wt blocks
/
&wt TYPE='TEMP0', istep1=0,  istep2=14286,  value1=0.1,  value2=300.0 /
&wt TYPE='END'                         /
""",
    'density-imp.in': """equilibrate complex
&cntrl
  imin=0,        ! MD on
  irest=1,       ! restart (read velocities)
  ntx=5,         ! read coordinates + velocities

  nstlim=1428,  ! total steps
  dt=0.0035,     ! timestep (ps)

  ntc=2,         ! SHAKE on bonds with H
  ntf=2,         ! no force eval on those bonds

  ig=-1,         ! random seed

  ntt=3,         ! Langevin thermostat
  gamma_ln=2.0,  ! collision freq (ps⁻¹)
  temp0=300.0,   ! target temp (K)

  ntpr=2,      ! print to mdout every 286 steps
  ntwx=2,      ! write coordinates every 286 steps

  ntb=0,         ! no periodic box
  igb=8,         ! GB model
  saltcon=0.145, ! 0.145 M salt

  cut=999.0,     ! GB nonbonded cutoff


  ntr=1,                             ! restraints on
  restraintmask='(:1-{last_residue}@CA,N,C|:{ligand_residue})',
  restraint_wt=2.0,                  ! same as explicit
/
""",
    'equil-imp.in': """equilibrate complex
&cntrl
  imin=0,        ! MD on
  irest=1,       ! restart (read velocities)
  ntx=5,         ! read coordinates + velocities

  nstlim=2857142,  ! total steps
  dt=0.0035,     ! timestep (ps)

  ntc=2,         ! SHAKE on bonds with H
  ntf=2,         ! no force eval on those bonds

  ig=-1,         ! random seed

  ntt=3,         ! Langevin thermostat
  gamma_ln=2.0,  ! collision freq (ps⁻¹)
  temp0=300.0,   ! target temp (K)

  ntpr=2857,      ! print to mdout every 286 steps
  ntwx=2857,      ! write coordinates every 286 steps

  ntb=0,         ! no periodic box
  igb=8,         ! GB model
  saltcon=0.145, ! 0.145 M salt

  cut=999.0,     ! GB nonbonded cutoff


  ntr=1,                             ! restraints on
  restraintmask='(:1-{last_residue}@CA,N,C|:{ligand_residue})',
  restraint_wt=0.02,                  ! same as explicit
/
/
""",
    'md-imp.in': """production complex
&cntrl
  imin=0,        ! MD on
  irest=1,       ! restart (read velocities)
  ntx=5,         ! read coordinates + velocities

  nstlim=28571429,  ! total steps
  dt=0.0035,     ! timestep (ps)

  ntc=2,         ! SHAKE on bonds with H
  ntf=2,         ! no force eval on those bonds

  ig=-1,         ! random seed

  ntt=3,         ! Langevin thermostat
  gamma_ln=0.01,  ! collision freq (ps⁻¹)
  temp0=300.0,   ! target temp (K)

  ntpr=28571,      ! print to mdout every 286 steps
  ntwx=28571,      ! write coordinates every 286 steps

  ntb=0,         ! no periodic box
  igb=8,         ! GB model
  saltcon=0.145, ! 0.145 M salt

  cut=999.0,     ! GB nonbonded cutoff
     ! GB screening cutoff

  ntr=1,                             ! restraints on
  restraintmask='(:1-{last_residue}@CA,N,C|:{ligand_residue})',
  restraint_wt=0.02,                  ! same as explicit
/
"""

}

# -----------------------
# Helper Functions
# -----------------------
def get_last_residue_number(pdb_filename):
    with open(pdb_filename, 'r') as file:
        for line in file:
            if line.startswith('TER'):
                parts = line.split()
                try:
                    return int(parts[3])
                except (IndexError, ValueError):
                    continue
    return None

def create_md_input_files(last_residue, input_files, out_dir):
    for filename, content in input_files.items():
        updated_content = content.format(last_residue=last_residue, ligand_residue=last_residue + 1)
        with open(os.path.join(out_dir, filename), 'w') as file:
            file.write(updated_content)
            print(f"File written: {os.path.join(out_dir, filename)}")

def create_pdb(top_file, crd_file, output_pdb):
    command = f"ambpdb -p {top_file} < {crd_file} > {output_pdb}"
    subprocess.run(command, shell=True)
    print(f"PDB file created: {output_pdb}")

def perform_hmr(top_file, modified_top_file):
    parmed_commands = f"""
hmassrepartition dowater
outparm {modified_top_file}
quit
"""
    with open("parmed.in", "w") as f:
        f.write(parmed_commands)
    subprocess.run(f"parmed -O -p {top_file} -i parmed.in", shell=True)
    os.remove("parmed.in")
    print(f"HMR topology written: {modified_top_file}")

# -----------------------
# Main Driver
# -----------------------
def main():
    if len(sys.argv) != 3:
        print("Usage: python script.py <input_directory> <pdb_list.csv>")
        sys.exit(1)

    input_dir = os.path.abspath(sys.argv[1])
    csv_file = os.path.abspath(sys.argv[2])

    # Load CSV and extract list of PDB codes
    try:
        df = pd.read_csv(csv_file)
        pdb_codes = df['PDB codes'].astype(str).tolist()
    except Exception as e:
        print(f"Error reading CSV file: {e}")
        sys.exit(1)

    for pdb_id in pdb_codes:
        top_pattern = os.path.join(input_dir, f"{pdb_id}.top")
        crd_pattern = os.path.join(input_dir, f"{pdb_id}.crd")
        top_files = glob.glob(top_pattern)
        crd_files = glob.glob(crd_pattern)

        if not top_files or not crd_files:
            print(f"Skipping {pdb_id}: Missing .top or .crd file")
            continue

        top_file = top_files[0]
        crd_file = crd_files[0]

        output_dir = os.path.join(os.getcwd(), pdb_id)
        os.makedirs(output_dir, exist_ok=True)

        output_pdb = os.path.join(output_dir, f"{pdb_id}.pdb")
        modified_top = os.path.join(output_dir, f"{pdb_id}_hmr.top")
        copied_crd = os.path.join(output_dir, f"{pdb_id}.crd")
        refc_file = os.path.join(output_dir, "refc")

        try:
            print(f"\nProcessing {pdb_id}...")

            create_pdb(top_file, crd_file, output_pdb)
            shutil.copy(crd_file, copied_crd)
            shutil.copy(crd_file, refc_file)
            perform_hmr(top_file, modified_top)

            last_residue = get_last_residue_number(output_pdb)
            if last_residue:
                create_md_input_files(last_residue, input_files, output_dir)
            else:
                print(f"Could not find TER line in {output_pdb}")

        except Exception as e:
            print(f"Error processing {pdb_id}: {e}")

if __name__ == "__main__":
    main()
