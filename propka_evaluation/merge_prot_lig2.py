import sys
import os

# Launch PyMOL in headless mode
import pymol
pymol.pymol_argv = ['pymol', '-qc']  # quiet + no GUI
pymol.finish_launching()

from pymol import cmd

def align_and_merge(pdb_file, mol2_file, ligand_resn, output_pdb):
    cmd.reinitialize()

    # Load files
    print(f"Loading files...", flush=True)
    cmd.load(pdb_file, "complex")
    cmd.load(mol2_file, "ligand_mol2")

    # Select the ligand from the complex using resn
    ligand_sel = f"complex and resn {ligand_resn} and not solvent"
    cmd.select("ligand_pdb", ligand_sel)
    if cmd.count_atoms("ligand_pdb") == 0:
        print(f"[ERROR] No atoms found for selection: {ligand_sel}", flush=True)
        return

    # Align mol2 ligand to pdb ligand
    print(f"Aligning ligand...", flush=True)
    cmd.align("ligand_mol2", "ligand_pdb")

    # Select protein (exclude water and ligand)
    protein_sel = f"complex and polymer.protein and not resn {ligand_resn}"
    cmd.select("protein_only", protein_sel)
    if cmd.count_atoms("protein_only") == 0:
        print(f"[ERROR] No atoms found for protein selection: {protein_sel}", flush=True)
        return

    # Create separate objects for protein and ligand, then merge
    print("Creating combined object with correct order...", flush=True)
    cmd.create("prot_obj", "protein_only")
    cmd.create("lig_obj", "ligand_mol2")
    cmd.create("combined", "prot_obj")
    cmd.load("combined", "ligand_mol2")

    # Save the result
    print(f"Saving combined object to {output_pdb}...", flush=True)
    cmd.save(output_pdb, "combined")
    print("Done!", flush=True)

    # Cleanup
    cmd.delete("all")



if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python merge_prot_lig.py <complex.pdb> <ligand.mol2> <LIGAND_RESN> <output.pdb>")
        sys.exit(1)

    pdb_file = sys.argv[1]
    mol2_file = sys.argv[2]
    ligand_resn = sys.argv[3]
    output_pdb = sys.argv[4]

    # Expand tilde to absolute path
    pdb_file = os.path.expanduser(pdb_file)
    mol2_file = os.path.expanduser(mol2_file)
    output_pdb = os.path.expanduser(output_pdb)

    align_and_merge(pdb_file, mol2_file, ligand_resn, output_pdb)
