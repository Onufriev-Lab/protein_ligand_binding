import os
import shutil
from pathlib import Path

# Define the root directory
base_path = Path.home() / 'full-imp-pipeline-work/gold-training-062725-heat-fix'

# Loop through X/Y/mbondi
for x_dir in base_path.iterdir():
    if not x_dir.is_dir():
        continue

    for y_dir in x_dir.iterdir():
        if not y_dir.is_dir():
            continue

        pdb_code = y_dir.name  # Now matches Y
        top_file = Path.cwd() / f"{pdb_code}.top"
        if not top_file.exists():
            print(f"WARNING: {top_file.name} not found in current directory. Skipping.")
            continue

        mbondi_dir = y_dir / 'mbondi'
        if mbondi_dir.is_dir():
            dest_file = mbondi_dir / f"{pdb_code}.top"
            shutil.copy2(top_file, dest_file)
            print(f"Copied {top_file.name} to {dest_file}")
        else:
            print(f"WARNING: Directory not found: {mbondi_dir}")
