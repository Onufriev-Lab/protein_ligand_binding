import os
import shutil
import glob
import sys

def collect_and_rename_files(root_dir, prefix="0.15_80_10_pH7_"):
    # Get absolute path to the root directory
    root_dir = os.path.abspath(root_dir)

    # Determine the new target directory at the same level
    parent_dir = os.path.dirname(root_dir)
    target_dir = os.path.join(parent_dir, f"{prefix.strip('_')}_collection")
    os.makedirs(target_dir, exist_ok=True)

    # Loop through all subdirectories (one level deep)
    subdirs = [d for d in glob.glob(os.path.join(root_dir, '*')) if os.path.isdir(d)]

    for subdir in subdirs:
        # Look for .top and .crd files
        for pattern in ('*.top', '*.crd'):
            files = glob.glob(os.path.join(subdir, pattern))
            for filepath in files:
                filename = os.path.basename(filepath)
                new_filename = f"{prefix}{filename}"
                new_filepath = os.path.join(target_dir, new_filename)

                print(f"Copying {filepath} to {new_filepath}")
                shutil.copy2(filepath, new_filepath)

    print(f"\nAll matching files copied to: {target_dir}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python collect_files.py <root_directory>")
        sys.exit(1)

    root_dir = sys.argv[1]
    collect_and_rename_files(root_dir)
