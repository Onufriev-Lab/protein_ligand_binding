def clean_pdb(input_pdb, output_pdb):
    with open(input_pdb, 'r') as infile:
        lines = infile.readlines()

    # Step 1: Fix 'ATOMNME' lines
    for i, line in enumerate(lines):
        if line.startswith("ATOMNME"):
            lines[i] = "ATOM" + "   " + line[10:]

    # Step 2: Find indices of all 'TER' lines
    ter_indices = [i for i, line in enumerate(lines) if line.startswith("TER")]

    if len(ter_indices) >= 3:
        target_ter_index = ter_indices[-3]  # third-to-last TER line
        # Remove the 3 lines *after* this TER line
        del lines[target_ter_index + 1: target_ter_index + 4]

    # Step 3: Write output
    with open(output_pdb, 'w') as outfile:
        outfile.writelines(lines)

if __name__ == "__main__":
    clean_pdb('3nq9-8-trunc_cap.pdb', '3nq9-8-trunc_cap_cleaned.pdb')