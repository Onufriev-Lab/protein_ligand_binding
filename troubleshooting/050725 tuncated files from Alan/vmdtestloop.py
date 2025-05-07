from vmd import molecule, atomsel
import subprocess
import sys
import os
import pandas as pd


#Dependency: vmdpython- https://vmd.robinbetz.com/ 
def process_protein(protein, ligand):
    
    #change "fullsetligands/" the folder where your pdbs are
    molecule.load("pdb", "fullsetligands/" + protein + ".pdb")
    big_sel = atomsel(f'protein or resname {ligand}')
    lig_sel = atomsel(f'resname {ligand}')
    ligand_in_protein_sasa = big_sel.sasa(srad=1.4, restrict=lig_sel)
    
    orig_sasa = ligand_in_protein_sasa
    print(ligand_in_protein_sasa)
    new_sasa = None
    radius = 0
    sasa_distances = []
    atom_numbers = []
    while new_sasa != orig_sasa and radius < 20:
        radius +=1
        #this runs my truncation script, "TruncAndCapLinux2.sh"
        subprocess.run(["bash", "TruncAndCapLinux2.sh", protein, ligand, str(radius)], capture_output=True, text=True)
        print(f"{protein}-{radius}.pdb")
        mol = molecule.load("pdb", "trunc_cap.pdb")
        big_sel = atomsel(f'protein or resname {ligand}')
        lig_sel = atomsel(f'resname {ligand}')
        ligand_in_protein_sasa = big_sel.sasa(srad=1.4, restrict=lig_sel)
        new_sasa = ligand_in_protein_sasa
        atom_numbers.append(molecule.numatoms(mol))
        sasa_distances.append(new_sasa)
        print(new_sasa)
        molecule.delete(mol)
    line_count = 0
    
    #counts number of residues in result
    with open('cpptraj_reslistnohoh.txt', 'r') as file:
        for line in file:
            line_count += 1
    file.close()
    print(radius, orig_sasa, sasa_distances, atom_numbers, line_count)
    return (radius, orig_sasa, sasa_distances, atom_numbers, line_count)

if __name__ == '__main__':
    # Define the dataset as a list of (protein, ligand) pairs.
    
    #change this to whatever your csv of pdbs is
    df = pd.read_csv('pdb_refined_w_AA_lengths2_Monomers_noIons_passed_first_screen.csv',
                    usecols=["PDB code", "Ligand Name"])

    # List to store each protein-ligand processing result
    results = []
    
    # Iterate over each row in the dataframe
    for idx, row in df.iterrows():
        protein = row["PDB code"]
        ligand = row["Ligand Name"]
        print(f"Processing protein {protein} with ligand {ligand} (row {idx})")
        try:
            result = process_protein(protein, ligand)
        except Exception as e:
            print(f"Error processing protein {protein} with ligand {ligand}: {e}")
            result = "error"
        results.append((protein, ligand, result[0], result[1],result[2],result[3],result[4]))
    
    # Convert the list of dictionaries into a DataFrame
    results_df = pd.DataFrame(results)
    
    # Write the results DataFrame to a new CSV file
    results_df.to_csv("processed_results.csv", index=False)
    
    
  
