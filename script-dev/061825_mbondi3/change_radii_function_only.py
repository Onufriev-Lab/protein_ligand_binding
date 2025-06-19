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