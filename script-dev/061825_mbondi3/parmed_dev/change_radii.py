import os
import sys
import subprocess

def change_radii(top_file, modified_top_file):
    parmed_commands = f"""
changeRadii mbondi3
outparm {modified_top_file}
quit
"""
    with open("parmed.in", "w") as f:
        f.write(parmed_commands)

    subprocess.run(f"parmed -O -p {top_file} -i parmed.in", shell=True)
    os.remove("parmed.in")
    print(f"radii changed to mbondi3: {modified_top_file}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python change_radii_script.py /path/to/file.top")
        sys.exit(1)

    top_file = sys.argv[1]

    if not os.path.isfile(top_file):
        print(f"Error: File '{top_file}' not found.")
        sys.exit(1)

    # Extract base filename and define output path in current directory
    base_name = os.path.basename(top_file)
    modified_top_file = os.path.join(os.getcwd(), base_name)

    change_radii(top_file, modified_top_file)
