import requests
import pandas

excel_data_df = pandas.read_excel('pdb_refined_w_AA_lengths2.xlsx', sheet_name=None, skiprows=1)
sheet_df = excel_data_df["pdb_refined_w_AA_lengths2"]
codes = sheet_df["PDB code"]
initialcount = len(codes.index)
print(sheet_df.head(n=20))

for index, code in enumerate(codes):
    response = requests.get(f"https://data.rcsb.org/rest/v1/core/assembly/{code}/1")
    if response.status_code == 200:
        print(response.json()['rcsb_struct_symmetry'])
        stoichiometry = response.json()['rcsb_struct_symmetry'][0]['stoichiometry']
        
        if stoichiometry[0] != "A1" or len(stoichiometry) != 1 :
            sheet_df.drop(index, inplace=True)
    else:
        sheet_df.drop(index, inplace=True)
finalcount = len(sheet_df.index)
sheet_df.to_csv('pdb_refined_w_AA_lengths2_Monomers.csv')
print(sheet_df.head(n=20))

print(f"initial count: {initialcount} final count: {finalcount}")