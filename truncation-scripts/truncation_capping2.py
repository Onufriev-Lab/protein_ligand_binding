import os
import subprocess
from pymol import cmd

def truncate_and_cap(pdb_path, ligand_resname, cutoff_distance, output_prefix="truncated"):
    base = os.path.splitext(os.path.basename(pdb_path))[0]
    output_dir = os.path.dirname(os.path.abspath(pdb_path))
    truncated_pdb = os.path.join(output_dir, f"{output_prefix}.pdb")
    capped_pdb = os.path.join(output_dir, f"{output_prefix}_capped.pdb")

    # -------- Step 1: Truncate in PyMOL --------
    cmd.reinitialize()
    cmd.load(pdb_path, "complex")

    cmd.select("ligand", f"resn {ligand_resname}")
    if cmd.count_atoms("ligand") == 0:
        raise ValueError(f"Ligand {ligand_resname} not found in {pdb_path}")

    cmd.select("protein_trunc", f"byres (complex and polymer.protein within {cutoff_distance} of ligand)")
    cmd.select("final", "ligand or protein_trunc")

    if cmd.count_atoms("final") == 0:
        raise RuntimeError("No atoms selected for final structure")

    cmd.save(truncated_pdb, "final")
    print(f"[PyMOL] Saved truncated PDB: {truncated_pdb}")

    # -------- Step 2: Add ACE/NME caps in tleap --------
    tleap_input = f"""source leaprc.protein.ff19SB
mol = loadpdb "{truncated_pdb}"
set mol head ACE
set mol tail NME
savepdb mol "{capped_pdb}"
quit
"""
    
    leapin_path = os.path.join(output_dir, "cap_leap.in")
    with open(leapin_path, "w") as f:
        f.write(tleap_input)

    print(f"[tleap] Running tleap...")
    result = subprocess.run(["tleap", "-f", leapin_path], capture_output=True, text=True)

    if result.returncode != 0:
        print("[tleap ERROR]")
        print(result.stderr)
        raise RuntimeError("tleap failed")
    else:
        print("[tleap] Capping complete.")
        print(result.stdout)

    if not os.path.isfile(capped_pdb):
        raise FileNotFoundError("Capped PDB not generated.")

    print(f"[Output] Final capped PDB: {capped_pdb}")


if __name__ == "__main__":
    # Example usage
    truncate_and_cap(
        pdb_path="3nq9.pdb",         # Replace with your input PDB file
        ligand_resname="OCA",        # Replace with your ligand resname
        cutoff_distance=8.0,         # Ångström
        output_prefix="3nq9_trunc"
    )
