import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import AddHs
from rdkit.Chem import SDWriter
from rdkit.Chem import rdmolops



# Load the data from the pickle file
data = pd.read_pickle('samples_all.pkl')

# Create an empty DataFrame to store molecular name and energy values
energy_df = pd.DataFrame(columns=["Molecular Name", "MMFF94 Energy"])

for item in data:
    rdmol = item.rdmol
    # Compare with each molecule in the test set
    mol_name = item.smiles  # Default to SMILES representation

    # Add hydrogens to the molecule
    mol_h = AddHs(rdmol, addCoords=True)

    # Use MMFF94 force field to get the molecule properties
    props = AllChem.MMFFGetMoleculeProperties(mol_h)

    # Create a force field object for the molecule
    ff = AllChem.MMFFGetMoleculeForceField(mol_h, props)

    # Calculate the energy of the molecule
    energy = ff.CalcEnergy()

    # Append the result to the energy DataFrame
    energy_df = energy_df.append({"Molecular Name": mol_name, "MMFF94 Energy": energy}, ignore_index=True)

    # Save the molecule to an SDF file
    output_filename = f"{mol_name}.sdf"
    writer = SDWriter(output_filename)
    writer.write(mol_h)
    writer.close()

# Save the energy DataFrame to a CSV file
energy_df.to_csv("Geodiff_MMFF94energy.csv", index=False)
