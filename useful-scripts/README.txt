This file contains information about the files in the useful-scripts directory

TruncAndCapLinux2.sh is a shell script that truncates proteins (by residue) and caps them based on distance from the ligand (any atom of the ligand, not center).  It takes three arguments in the following order:  PDB code, ligand name (three letter code), radius to truncate (in Angstroms).  The files associated with the PDB code need to be in the same folder as the script. This script outputs a pdb file with the following name pdbcode-radius-trunc_cap.pdb (example: 3ary-7-trunc_cap.pdb).

fix_trun_capped2.py is a python script that fixes the output from TruncAndCapLinux2.sh.  It removes extraneouly added NMEs and ACEs.  Input file format is pdbcode-radius-trunc_cap.pdb as descdribed above.  Output file format is pdbcode-radius-trunc_cap_cleaned.pdb.

