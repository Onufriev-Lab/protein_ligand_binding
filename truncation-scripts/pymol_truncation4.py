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

        # Select ligand by residue name
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

        # Extract to new object and break into fragments
        cmd.create("trunc_obj", "truncation")
        cmd.fragment("trunc_obj")  # creates trunc_obj_0001, trunc_obj_0002, ...

        # Assign chain IDs to fragments
        frag_names = cmd.get_object_list("trunc_obj_*")
        chain_id = 65  # ASCII 'A'
        for frag in frag_names:
            chain_letter = chr(chain_id)
            cmd.alter(frag, f'chain="{chain_letter}"')
            chain_id += 1
            if chain_id > 90:  # Past 'Z'
                print("Warning: More than 26 fragments. Reusing chain IDs.")
                chain_id = 65

        # Create final combined object
        cmd.create("ligand_obj", "ligand")
        for frag in frag_names:
            cmd.create("ligand_obj", frag, 1, 1)

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
