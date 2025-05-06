from typing import List
import os

def truncate_and_cap_proteins(pdb_codes: List[str], ligand_codes: List[str], distances: List[float]):
    from pymol import cmd

    for pdb, ligand_code, cutoff in zip(pdb_codes, ligand_codes, distances):
        cmd.reinitialize()
        pdb_file = f"{pdb}.pdb"
        cmd.load(pdb_file, "complex")

        # Select ligand by residue name
        cmd.select("ligand", f"resn {ligand_code}")

        # Select all protein residues with at least one atom within cutoff of ligand
        cmd.select("close_residues", f"byres (polymer.protein within {cutoff} of ligand)")

        # Save selection as a new object
        cmd.create("truncated", "close_residues or ligand")

        # Cap N-terminus and C-terminus
        # WARNING: This assumes linear peptide chains and that ACE and NME are already available
        # You may need to add parameters for accurate modeling after export (e.g., in Amber)
        cmd.remove("not (truncated)")
        cmd.sort()

        # Identify terminal residues
        cmd.iterate("truncated and name N", "stored.nterm = resi", space=globals())
        cmd.iterate("truncated and name OXT", "stored.cterm = resi", space=globals())
        nterm = globals().get("stored.nterm", None)
        cterm = globals().get("stored.cterm", None)

        if nterm:
            cmd.fab("NME", "nme_cap")
            cmd.align("nme_cap and name C", f"truncated and resi {nterm} and name N")
            cmd.fuse("nme_cap", f"truncated and resi {nterm} and name N", mode=1)
            cmd.delete("nme_cap")

        if cterm:
            cmd.fab("ACE", "ace_cap")
            cmd.align("ace_cap and name N", f"truncated and resi {cterm} and name OXT")
            cmd.fuse("ace_cap", f"truncated and resi {cterm} and name OXT", mode=1)
            cmd.delete("ace_cap")

        # Save result
        output_file = f"{pdb}_truncated_capped.pdb"
        cmd.save(output_file, "truncated")
        print(f"Saved: {output_file}")

# Example usage (provide real PDB files beforehand):
truncate_and_cap_proteins(["3nq9", "4gid"], ["OCA", "0GH"], [8.0, 8.0])
