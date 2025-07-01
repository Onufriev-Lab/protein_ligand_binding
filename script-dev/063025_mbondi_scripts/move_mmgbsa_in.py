import shutil
from pathlib import Path

# Source file to copy
source_file = Path.cwd() / 'mmgbsa.in'
if not source_file.exists():
    raise FileNotFoundError(f"{source_file} does not exist in the current directory.")

# Define base path
base_path = Path.home() / 'full-imp-pipeline-work/gold-training-062725-heat-fix'

# Traverse and copy mmgbsa.in into each mbondi directory
for x_dir in base_path.iterdir():
    if not x_dir.is_dir():
        continue

    for y_dir in x_dir.iterdir():
        mbondi_dir = y_dir / 'mbondi'
        if mbondi_dir.is_dir():
            dest_file = mbondi_dir / 'mmgbsa.in'
            shutil.copy2(source_file, dest_file)
            print(f"Copied mmgbsa.in to {dest_file}")
        else:
            print(f"WARNING: mbondi directory does not exist: {mbondi_dir}")
