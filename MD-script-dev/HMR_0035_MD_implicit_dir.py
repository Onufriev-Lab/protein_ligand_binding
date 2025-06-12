import subprocess
import glob
import re
import shutil
import os
import sys

# Input file templates (keep your full input_files dict here)
input_files = {
    'min.in': """minimise complex
&cntrl
  imin=1,        ! do minimisation
  maxcyc=1000,   ! total cycles
  ncyc=500,      ! switch from steepest descent to conjugate‑gradient
  cut=999.0,     ! nonbonded cutoff (large for GB)
  rgbmax=10.0,  ! GB cutoff
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
    'heat.in': """heat complex
&cntrl
  imin=0,        ! MD heating
  irest=0,       ! start from scratch
  ntx=1,         ! coordinates only, no velocities
  nstlim=14286,  ! total steps (≈50 ps @3.5 fs)
  dt=0.0035,     ! timestep (ps)
  
  ntc=2,         ! SHAKE on bonds involving H
  ntf=2,         ! no forces on bonds with SHAKE
  
  ig=-1,         ! random seed
  ntt=1,         ! Berendsen thermostat
  tautp=1.0,     ! coupling time (ps)
  tempi=0.0,     ! start temp (K)
  temp0=300.0,   ! target temp (K)
  
  ntpr=286,      ! print to mdout every 286 steps
  ntwx=286,      ! write coords every 286 steps
  
  ntb=0,         ! no periodic box
  igb=8,         ! GB model
  saltcon=0.145, ! 0.145 M salt
  
  cut=999.0,     ! nonbonded cutoff for GB
  rgbmax=10.0,  ! GB screening cutoff
  
  ntr=1,         ! positional restraints on
  restraintmask='(:1-78@CA,N,C|:79)',
  restraint_wt=2.0,
  
  nmropt=1,      ! enable NMR-style &wt blocks
/
&wt TYPE='TEMP0', istep1=0,  istep2=14286,  value1=0.1,  value2=300.0 /
&wt TYPE='END'                         /&cntrl
  imin=0,        ! MD heating
  irest=0,       ! start from scratch
  ntx=1,         ! coordinates only, no velocities
  nstlim=14286,  ! total steps (≈50 ps @3.5 fs)
  dt=0.0035,     ! timestep (ps)
  
  ntc=2,         ! SHAKE on bonds involving H
  ntf=2,         ! no forces on bonds with SHAKE
  
  ig=-1,         ! random seed
  ntt=1,         ! Berendsen thermostat
  tautp=1.0,     ! coupling time (ps)
  tempi=0.0,     ! start temp (K)
  temp0=300.0,   ! target temp (K)
  
  ntpr=286,      ! print to mdout every 286 steps
  ntwx=286,      ! write coords every 286 steps
  
  ntb=0,         ! no periodic box
  igb=8,         ! GB model
  saltcon=0.145, ! 0.145 M salt
  
  cut=999.0,     ! nonbonded cutoff for GB
  rgbmax=10.0,  ! GB screening cutoff
  
  ntr=1,         ! positional restraints on
  restraintmask='(:1-{last_residue}@CA,N,C|:{ligand_residue})',
  restraint_wt=2.0,
  
  nmropt=1,      ! enable NMR-style &wt blocks
/
&wt TYPE='TEMP0', istep1=0,  istep2=14286,  value1=0.1,  value2=300.0 /
&wt TYPE='END'                         /
""",
    'density.in': """equilibrate complex
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
  rgbmax=10.0,  ! GB screening cutoff

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
  rgbmax=10.0,  ! GB screening cutoff

  ntr=1,                             ! restraints on
  restraintmask=':1-{last_residue}@CA,N,C',
  restraint_wt=2.0,                  ! same as explicit
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
  rgbmax=10.0,  ! GB screening cutoff

  ntr=1,                             ! restraints on
  restraintmask=':1-{last_residue}@CA,N,C',
  restraint_wt=2.0,                  ! same as explicit
/
"""

}

def find_all_input_file_pairs(input_dir):
    top_files = glob.glob(os.path.join(input_dir, "*.top"))
    crd_files = glob.glob(os.path.join(input_dir, "*.crd"))

    if not top_files or not crd_files:
        raise FileNotFoundError("No .top or .crd files found in the specified directory.")

    tops_by_id = {}
    crds_by_id = {}

    # Index .top files
    for top in top_files:
        match = re.search(r"pH7_([a-zA-Z0-9]+)\.top$", os.path.basename(top))
        if match:
            pdb_id = match.group(1)
            tops_by_id[pdb_id] = top

    # Index .crd files
    for crd in crd_files:
        match = re.search(r"pH7_([a-zA-Z0-9]+)\.crd$", os.path.basename(crd))
        if match:
            pdb_id = match.group(1)
            crds_by_id[pdb_id] = crd

    # Return list of matched pairs
    matched_pairs = []
    for pdb_id in sorted(set(tops_by_id.keys()) & set(crds_by_id.keys())):
        matched_pairs.append((pdb_id, tops_by_id[pdb_id], crds_by_id[pdb_id]))

    if not matched_pairs:
        raise ValueError("No matching pairs of .top and .crd files found.")

    return matched_pairs


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
    print(f"Modified topology file created: {modified_top_file}")

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

def main():
    if len(sys.argv) != 2:
        print("Usage: python script.py <input_directory>")
        sys.exit(1)

    input_dir = os.path.abspath(sys.argv[1])

    try:
        file_pairs = find_all_input_file_pairs(input_dir)

        for pdb_id, top_file, crd_file in file_pairs:
            print(f"Processing {pdb_id}...")

            output_dir = os.path.join(os.getcwd(), pdb_id)
            os.makedirs(output_dir, exist_ok=True)

            output_pdb = os.path.join(output_dir, f"{pdb_id}.pdb")
            modified_top = os.path.join(output_dir, f"{pdb_id}_hmr.top")
            copied_crd = os.path.join(output_dir, f"{pdb_id}.crd")
            refc_file = os.path.join(output_dir, "refc")

            # Create PDB
            create_pdb(top_file, crd_file, output_pdb)

            # Copy CRD
            shutil.copy(crd_file, copied_crd)
            shutil.copy(crd_file, refc_file)

            # HMR
            perform_hmr(top_file, modified_top)

            # Create input files
            last_residue = get_last_residue_number(output_pdb)
            if last_residue:
                create_md_input_files(last_residue, input_files, output_dir)
            else:
                print(f"Warning: Could not determine last residue from {output_pdb}")

    except Exception as e:
        print(f"Error: {e}")


if __name__ == "__main__":
    main()
