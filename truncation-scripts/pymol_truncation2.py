from typing import List
import os
from pymol import cmd

def truncate_proteins_by_distance(
    pdb_codes: List[str],
    ligand_resnames: List[str],
    distances: List[float]
):
    script_dir = os.path.dirname(os.path.abspath(__file__))

    for pdb_code, lig_code, dist in zip(pdb_codes, ligand_resnames, distances):
        cmd.reinitialize()

        pdb_path = os.path.join(script_dir, f"{pdb_code}.pdb")
        if not os.path.isfile(pdb_path):
            print(f"File not found: {pdb_path}")
            continue

        cmd.load(pdb_path, "complex")

        # Select the ligand by residue name
        ligand_sel = f"complex and resn {lig_code}"
        cmd.select("ligand", ligand_sel)
        if cmd.count_atoms("ligand") == 0:
            print(f"Warning: No ligand atoms found for resn {lig_code} in {pdb_code}")
            continue

        # Select full protein residues within distance of ligand
        trunc_sel = f"byres (complex and polymer.protein within {dist} of ligand)"
        cmd.select("truncation", trunc_sel)
        if cmd.count_atoms("truncation") == 0:
            print(f"Warning: No protein atoms found within {dist} Å of ligand in {pdb_code}")
            continue

        # Split discontinuous fragments into separate objects and assign unique chains
        cmd.create("fragments", "truncation")
        num_frags = cmd.count_fragments("fragments")

        frag_objects = []
        for i in range(num_frags):
            frag_name = f"frag_{i}"
            cmd.fragment("fragments", frag_name, i + 1)  # i+1 since fragments are 1-indexed
            cmd.alter(frag_name, f'chain="{chr(65 + i % 26)}"')  # A-Z cycling
            frag_objects.append(frag_name)

        # Create final object
        cmd.create("ligand_obj", "ligand")
        for frag in frag_objects:
            cmd.create("ligand_obj", f"{frag}", 1, 1)

        output_path = os.path.join(script_dir, f"{pdb_code}_truncated.pdb")
        if cmd.count_atoms("ligand_obj") > 0:
            cmd.save(output_path, "ligand_obj")
            print(f"Saved: {output_path}")
        else:
            print(f"Warning: Empty selection for {pdb_code}. No file saved.")


if __name__ == "__main__":
    truncate_proteins_by_distance(
        pdb_codes=["3nq9", "4gid"],
        ligand_resnames=["OCA", "0GH"],
        distances=[8.0, 8.0]
    )
