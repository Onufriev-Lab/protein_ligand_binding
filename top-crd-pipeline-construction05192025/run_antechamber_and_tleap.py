import subprocess
import sys

def run_pipeline(pdb_code, charge):
    # Step 1: Run antechamber.py
    try:
        subprocess.run(
            ['python', 'antechamber.py', pdb_code, str(charge)],
            check=True
        )
        print(f"antechamber.py completed for {pdb_code}")
    except subprocess.CalledProcessError as e:
        print(f"Error running antechamber.py: {e}")
        return

    # Step 2: Run tleap.py
    try:
        subprocess.run(
            ['python', 'tleap.py', pdb_code],
            check=True
        )
        print(f"tleap.py completed for {pdb_code}")
    except subprocess.CalledProcessError as e:
        print(f"Error running tleap.py: {e}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python run_pipeline.py <pdb_code> <charge>")
        sys.exit(1)

    pdb_code = sys.argv[1]
    charge = sys.argv[2]

    run_pipeline(pdb_code, charge)
