This file contains information about the files in the useful-scripts directory

TruncAndCapLinux2.sh is a shell script that truncates proteins (by residue) and caps them based on distance from the ligand (any atom of the ligand, not center).  It takes three arguments in the following order:  PDB code, ligand name (three letter code), radius to truncate (in Angstroms).  The files associated with the PDB code need to be in the same folder as the script. This script outputs a pdb file with the following name pdbcode-radius-trunc_cap.pdb (example: 3ary-7-trunc_cap.pdb).  WARNING: This script does not currently work correctly if the ligand name is a number.

fix_trun_capped2.py is a python script that fixes the output from TruncAndCapLinux2.sh.  It removes extraneouly added NMEs and ACEs.  Input file format is pdbcode-radius-trunc_cap.pdb as descdribed above.  Output file format is pdbcode-radius-trunc_cap_cleaned.pdb.

trunc_cap_fix.py is a python script that takes a .csv file and extracts pdb codes, ligand names, and radii from columns to run TruncAndCapLinux2.sh and fix_trun_capped2.py to produce pdbcode-radius-trunc_cap_cleaned.pdb files in a batchwise manner.  The name of the column containing the pdb codes needs to be "PDB code".  The name of the ligand name column needs to be "Ligand Name".  The name of the radii column needs to be "trunc_rad".

one_csv_to_rule_them_all.py is a python script that collects unique .top/.crd filenames (without the extension), and creates a .csv file with a column named "PDB codes" that contains them.  
Usage:  python one_csv_to_rule_them_all.py directory_containing_top_crd_files

make_csvs.py is a python script that collects unique .top/.crd filenames (without the extension) and divides them evenly between .csv file in a column named "PDB codes".  THe number of .csv files is specified in the input.  The output files are named group_n.csv where n starts counting at 0 to easily match GPU ids.
Usage:  python make_csvs.py directory_containing_top_crd_files integer_number_of_csv_files

fix_mislabeled_AA is a python script that renames amino acids that have been named ACE or NME by the truncation script.  It does not change ACE and NME on the ends of peptides.
Usage:  python fix_mislabeled_AA.py trunc_pdb_file
