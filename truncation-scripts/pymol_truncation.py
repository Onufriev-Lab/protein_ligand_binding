from typing import List
import os
from pymol import cmd

def truncate_proteins_by_distance(
    pdb_codes: List[str],
    ligand_resnames: List[str],
    distances: List[float]
):
    # Get the directory containing the script
    script_dir = os.path.dirname(os.path.abspath(__file__))

    for pdb_code, lig_code, dist in zip(pdb_codes, ligand_resnames, distances):
        cmd.reinitialize()

        # Build input file path and load
        pdb_path = os.path.join(script_dir, f"{pdb_code}.pdb")
        cmd.load(pdb_path, "complex")

        # Select the ligand
        cmd.select("ligand", f"resn {lig_code}")

        # Select protein residues within cutoff distance of ligand
        cmd.select("truncation", f"byres (polymer.protein within {dist} of ligand)")

        # Include the ligand in the output
        cmd.select("final", "truncation or ligand")

        # Save to output file
        output_path = os.path.join(script_dir, f"{pdb_code}_truncated.pdb")
        cmd.save(output_path, "final")
        print(f"Saved: {output_path}")

if __name__ == "__main__":
    truncate_proteins_by_distance(
        pdb_codes=["3nq9", "4gid"],
        ligand_resnames=["OCA", "0GH"],
        distances=[8.0, 8.0]
    )

