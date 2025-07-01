import os
from pathlib import Path

# Replace these with your actual values
x_values = ['mb3', 'mb3-gbsa3', 'mb2', 'mb2-gbsa3']  # Replace with actual names
y_values = ['4mre', '6r9u', '6oqc']        # Replace with actual names

base_path = Path.home() / 'full-imp-pipeline-work/gold-training-062725-heat-fix'

for x in x_values:
    for y in y_values:
        mbondi_path = base_path / x / y / 'mbondi'
        mbondi_path.mkdir(parents=True, exist_ok=True)
        print(f"Created: {mbondi_path}")
