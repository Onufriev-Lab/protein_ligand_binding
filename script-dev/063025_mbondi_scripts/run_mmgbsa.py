import subprocess
from pathlib import Path
import glob
from concurrent.futures import ProcessPoolExecutor, as_completed

# Define base path
base_path = Path.home() / 'full-imp-pipeline-work/gold-training-062725-heat-fix'

# Gather all mbondi directories with corresponding mdcrd file
tasks = []
for x_dir in base_path.iterdir():
    if not x_dir.is_dir():
        continue

    for y_dir in x_dir.iterdir():
        mbondi_dir = y_dir / 'mbondi'
        if not mbondi_dir.is_dir():
            continue

        # Find the trajectory file in the parent directory
        traj_files = glob.glob(str(y_dir / '*_md1.mdcrd'))
        if not traj_files:
            print(f"WARNING: No *_md1.mdcrd found in {y_dir}")
            continue

        traj_file = Path(traj_files[0]).name
        tasks.append((mbondi_dir, traj_file))


def run_mmgbsa(mbondi_dir, traj_file):
    cmd = [
        'MMPBSA.py', '-O',
        '-i', 'mmgbsa.in',
        '-cp', 'complex.prmtop',
        '-rp', 'protein.prmtop',
        '-lp', 'ligand.prmtop',
        '-y', f'../{traj_file}'
    ]

    try:
        subprocess.run(cmd, cwd=mbondi_dir, check=True)
        return f"✅ Success in {mbondi_dir}"
    except subprocess.CalledProcessError as e:
        return f"❌ Error in {mbondi_dir}: {e}"


# Run in parallel
max_workers = min(12, len(tasks))  # Adjust based on available CPU cores
with ProcessPoolExecutor(max_workers=max_workers) as executor:
    futures = [executor.submit(run_mmgbsa, mbondi_dir, traj_file) for mbondi_dir, traj_file in tasks]
    for future in as_completed(futures):
        print(future.result())
