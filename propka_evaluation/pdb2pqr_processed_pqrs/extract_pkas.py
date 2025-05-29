import pandas as pd
import re


# Placeholder empty DataFrame
df = pd.DataFrame(columns=['amino acid', 'pKa'])

# Dictionary to convert 3-letter amino acid codes to 1-letter codes
three_to_one = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
    'GLX': 'Z', 'ASX': 'B', 'UNK': 'X'
}

# Function to extract amino acid and pKa from a log file
def extract_pka_data(log_text):
    pattern = re.compile(r'^\s*([A-Z]{3})\s+(\d+)\s+[A-Z]\s+([\d\.]+[*]?)')
    rows = []

    for line in log_text.splitlines():
        match = pattern.match(line)
        if match:
            aa3, pos, pka = match.groups()
            aa1 = three_to_one.get(aa3.upper(), '?')
            if '*' in pka:
                pka = pka.replace('*', '')
            rows.append({'amino acid': f'{aa1}{pos}', 'pKa': float(pka)})

    return pd.DataFrame(rows)

# This function expects the content of the log to be passed as a string to `log_text`
# Example usage:
with open('example.log') as f:
    log_content = f.read()
    df = extract_pka_data(log_content)




