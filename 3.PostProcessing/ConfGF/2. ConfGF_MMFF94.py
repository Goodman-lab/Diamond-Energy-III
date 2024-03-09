import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import AddHs
import os

def compute_energy_from_file(file_path):
    data = pd.read_pickle(file_path)

    if not hasattr(data, 'rdmol'):
        return None

    rdmol = data.rdmol
    mol_name = data.smiles.split("_")[0]
    conf_num = data.smiles.split("_")[1].split("H")[0]

    mol_h = AddHs(rdmol, addCoords=True)
    props = AllChem.MMFFGetMoleculeProperties(mol_h)
    ff = AllChem.MMFFGetMoleculeForceField(mol_h, props)
    energy = ff.CalcEnergy()

    return conf_num, energy


all_files = [f for f in os.listdir('.') if os.path.isfile(f) and f.endswith('.pkl')]

# Dictionary to store energy values for each molecule
molecules_dict = {}

for file_name in all_files:
    result = compute_energy_from_file(file_name)

    if result is not None:
        conf_num, energy = result
        molecule_name = file_name.split("_")[1]

        # Check if the molecule_name is already in the dictionary
        if molecule_name not in molecules_dict:
            molecules_dict[molecule_name] = pd.DataFrame(columns=["Conformer_Number", "MMFF94 Energy"])

        # Add energy value to the corresponding molecule entry in the dictionary
        molecules_dict[molecule_name].loc[len(molecules_dict[molecule_name])] = {"Conformer_Number": conf_num,
                                                                                 "MMFF94 Energy": energy}

# Save each molecule's energies to its own .csv file
for molecule_name, energy_df in molecules_dict.items():
    csv_name = molecule_name + ".csv"
    energy_df.to_csv(csv_name, index=False)
