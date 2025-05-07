import re

def clean_pdb(input_pdb, output_pdb):
    
    with open(input_pdb, 'r') as infile:
        lines = infile.readlines()

    # Step 1: Fix lines starting with 'ATOM' followed by one or more 'NME' within columns 1–17
    for i, line in enumerate(lines):
        prefix = line[0:17]
        if line.startswith("ATOM") and re.match(r'^ATOM(?:NME)+', prefix):
            # Replace all leading NMEs in columns 1–17 with 'ATOM   '
            new_prefix = re.sub(r'^ATOM(?:NME)+', 'ATOM', prefix)
            # Replace the original line with corrected prefix + remainder
            lines[i] = new_prefix + line[17:]


    # Step 2: Find all 'TER' line indices
    ter_indices = [i for i, line in enumerate(lines) if line.startswith("TER")]

    # Step 3: Apply condition and remove lines if met
    if len(ter_indices) >= 3:
        third_last_ter_idx = ter_indices[-3]
        second_last_ter_idx = ter_indices[-2]

        third_last_ter_line = lines[third_last_ter_idx]
        second_last_ter_line = lines[second_last_ter_idx]

        nme_condition = third_last_ter_line[17:20] == "NME"
        ace_condition = second_last_ter_line[17:20] == "ACE"

        if nme_condition and ace_condition:
            del lines[third_last_ter_idx + 1 : third_last_ter_idx + 4]

    # Step 4: Write cleaned file
    with open(output_pdb, 'w') as outfile:
        outfile.writelines(lines)

if __name__ == "__main__":
    clean_pdb('3nq9-8-trunc_cap.pdb', '3nq9-8-trunc_cap_cleaned.pdb')