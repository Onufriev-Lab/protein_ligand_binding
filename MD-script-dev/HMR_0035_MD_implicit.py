import subprocess
import glob
import re
import shutil

def find_input_files():
    # Find files matching the pattern for topology and coordinate files
    top_files = glob.glob("0.15_80_10_pH7_*.top")
    crd_files = glob.glob("0.15_80_10_pH7_*.crd")
    
    if not top_files or not crd_files:
        raise FileNotFoundError("Topology or coordinate files not found. Please ensure they follow the naming convention.")
    
    # Extract PDB ID from one of the file names using a more specific regex
    pdb_id_match = re.search(r"pH7_([a-zA-Z0-9]+)\.", top_files[0])
    if not pdb_id_match:
        raise ValueError("PDB ID could not be extracted from the file name.")
    
    pdb_id = pdb_id_match.group(1)
    top_file = top_files[0]
    crd_file = crd_files[0]
    
    return pdb_id, top_file, crd_file

def create_pdb(top_file, crd_file, output_pdb):
    # Command to convert top and crd files to a PDB file
    command = f"ambpdb -p {top_file} < {crd_file} > {output_pdb}"
    subprocess.run(command, shell=True)
    print(f"PDB file created: {output_pdb}")

def perform_hmr(top_file, modified_top_file):
    # Content of the parmed input file
    parmed_commands = f"""
hmassrepartition dowater
outparm {modified_top_file}
quit
"""
    # Write the parmed.in file
    with open("parmed.in", "w") as f:
        f.write(parmed_commands)

    # Command to run parmed for HMR
    command = f"parmed -O -p {top_file} -i parmed.in"
    subprocess.run(command, shell=True)
    print(f"Modified topology file created: {modified_top_file}")

def main():
    try:
        # Find input files and extract the PDB ID
        pdb_id, top_file, crd_file = find_input_files()
        crd_files = glob.glob('0.15_80_10_pH7_*.crd')
        # Define output file names
        output_pdb = f"{pdb_id}.pdb"
        modified_top_file = f"{pdb_id}.top"

        # Create PDB file
        create_pdb(top_file, crd_file, output_pdb)

        new_crd_name = f"{pdb_id}.crd"
        refc_name = "refc"

        for crd_file in crd_files:
            shutil.copy(crd_file, new_crd_name)
            print(f"CRD file copied to {new_crd_name}")
            shutil.copy(crd_file, refc_name)
            print(f"CRD file copied to {refc_name}")

        # Perform HMR
        perform_hmr(top_file, modified_top_file)
    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    main()
 

################

def get_last_residue_number(pdb_filename):
    last_residue_number = None
    with open(pdb_filename, 'r') as file:
        for line in file:
            if line.startswith('TER'):
                parts = line.split()
                last_residue_number = int(parts[3])  # Extract the residue number, which is typically the fourth element in the TER line
                break
    return last_residue_number

def create_md_input_files(last_residue, input_files):
    for filename, content in input_files.items():
        with open(filename, 'w') as file:
            updated_content = content.format(last_residue=last_residue, ligand_residue=last_residue + 1)
            file.write(updated_content)
            print(f'File {filename} has been written with updated restraints.')

# Define the templates for input files with placeholders for dynamic content
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

pdb_files = glob.glob('./*.pdb')

if pdb_files:
    for pdb_filename in pdb_files:
        last_residue = get_last_residue_number(pdb_filename)
        if last_residue:
            create_md_input_files(last_residue, input_files)
        else:
            print(f"Could not determine the last residue number from the PDB file: {pdb_filename}")
else:
    print("No PDB files found in the current directory.")
