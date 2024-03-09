mol_names=["c4", "c5", "c5n", "c6", "c6n", "c7", "c7c", "c7n", "c8", "c8a", "c8s", "c8si", "c9", "c9n", "c9s",
           "c10n", "c12p", "c12p0", "c12p1", "c15", "c15helix", "c15a", "c15b", "c15c", "c15d", "c15e", "c15f", "c15g", "c16", "c17",
           "c17l", "c18", "c18n", "c19", "c21", "c26", "c20n", "c39n",
           "e1", "e2", "e3", "e4", "e5", "e6", "e7", "e8", "e9", "e10",
           "e11", "e12", "e13", "e14", "e15", "e16", "e17", "e18", "e19", "e20",
           "e21", "e22", "e23", "e24", "e25", "e26", "e27", "e28", "e29", "e30",
           "e31", "e32", "e33", "e34", "e35", "e36", "e37", "e38", "e39", "e40"]

def mol_ID(mol_name):
    return mol_names.index(mol_name)
#print(len(mol_names))



#training_set = set(["c5", "c5n", "c6", "c6n", "c7", "c8", "c8s", "c8si", "c9", "c9n", "c9s",
#           "c10n", "c12p0", "c12p1", "c15", "c15a", "c15b", "c15c", "c15e", "c15f", "c15g", "c16", "c17", "c18", "c19", "c21", "c26",
#           "e1", "e2", "e3", "e4", "e5", "e6", "e8", "e10",
#           "e11", "e13", "e14", "e15", "e16", "e17", "e19",
#           "e21", "e22", "e23", "e24", "e25", "e26", "e28", "e29", "e30",
#           "e32", "e33", "e34", "e35", "e36", "e37"])


#validation_set = set(["c7c", "c15d", "c18n", "e9", "e7", "e18", "e20", "e31"])
#test_set = set(["c4", "c7n", "c8a", "c12p", "c15helix", "c17l", "c20n", "c39n", "e12", "e27", "e38", "e39", "e40"])


# Function to check for overlapping molecules among sets
#def check_unique_sets(*sets):
#    all_mols = set()
#    for s in sets:
#        if len(all_mols.intersection(s)) > 0:
#            overlapping_molecules = all_mols.intersection(s)
#            raise ValueError(f"Overlapping molecules found: {overlapping_molecules}")
#       all_mols.update(s)
# Check if the sets are unique
#check_unique_sets(training_set, validation_set, test_set)


# Generate the IDs based on the sets
#training_set_ID = [mol_ID(mol_name) for mol_name in training_set]
#validation_set_ID = [mol_ID(mol_name) for mol_name in validation_set]
#test_set_ID = [mol_ID(mol_name) for mol_name in test_set]
#print("training_set_ID:", training_set_ID)
#print("validation_set_ID:", validation_set_ID)
#print("test_set_ID:", test_set_ID)
#print("training_set_length:", len(training_set_ID))
#print("validation_set_length:", len(validation_set_ID))
#print("test_set_length:", len(test_set_ID))

test_set = set(["c4", "c7n", "c8a", "c12p", "c15helix", "c17l", "c20n", "c39n", "e12", "e27", "e38", "e39", "e40"])
# Define the test set and their InChI strings
test_set_inchi = {
    "c4":"InChI=1S/C4H10/c1-3-4-2/h3-4H2,1-2H3",
    "c7n":"InChI=1S/C7H16/c1-3-5-7-6-4-2/h3-7H2,1-2H3",
    "c8a":"InChI=1S/C8H18/c1-5-7(3)8(4)6-2/h7-8H,5-6H2,1-4H3/t7-,8+",
    "c12p":"InChI=1S/C12H26/c1-9(11(3,4)5)10(2)12(6,7)8/h9-10H,1-8H3",
    "c15helix":"InChI=1S/C15H32/c1-7-12(3)9-14(5)11-15(6)10-13(4)8-2/h12-15H,7-11H2,1-6H3/t12-,13-,14+,15+/m0/s1",
    "c17l":"InChI=1S/C17H36/c1-3-5-7-9-11-13-15-17-16-14-12-10-8-6-4-2/h3-17H2,1-2H3",
    "c20n":"InChI=1S/C20H42/c1-3-5-7-9-11-13-15-17-19-20-18-16-14-12-10-8-6-4-2/h3-20H2,1-2H3",
    "c39n":"InChI=1S/C39H80/c1-3-5-7-9-11-13-15-17-19-21-23-25-27-29-31-33-35-37-39-38-36-34-32-30-28-26-24-22-20-18-16-14-12-10-8-6-4-2/h3-39H2,1-2H3",
    "e12":"InChI=1S/C6H12O3/c1-5-3-8-4-6(2-7)9-5/h5-7H,2-4H2,1H3/t5-,6-/m1/s1",
    "e27":"InChI=1S/C9H18O2/c1-3-5-8-6-4-7-9(10-2)11-8/h8-9H,3-7H2,1-2H3/t8-,9-/m0/s1",
    "e38":"InChI=1S/C6H12O2/c1-5-2-3-6(7)4-8-5/h5-7H,2-4H2,1H3/t5-,6+/m1/s1",
    "e39":"InChI=1S/C26H44O10/c1-17-2-8-23(29-12-17)34-19-4-10-25(31-14-19)36-21-6-11-26(32-16-21)35-20-5-9-24(30-15-20)33-18-3-7-22(27)28-13-18/h17-27H,2-16H2,1H3/t17-,18-,19-,20-,21-,22-,23+,24+,25+,26+/m0/s1",
    "e40":"InChI=1S/C12H22O11/c13-1-3-5(15)6(16)9(19)12(22-3)23-10-4(2-14)21-11(20)8(18)7(10)17/h3-20H,1-2H2/t3-,4-,5+,6+,7-,8-,9-,10-,11-,12+/m1/s1",
}
print("len_test_set:", len(test_set))

import os
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import SDWriter



# Load the data
data = pd.read_pickle('ConfVAE_output.pkl')

# Convert InChI to molecules
name_to_mol = {name: Chem.MolFromInchi(inchi) for name, inchi in test_set_inchi.items()}

# Create a dictionary to keep track of conformer counts for each molecule
conformer_count = {name: 0 for name in name_to_mol.keys()}

for mol in data:
    mol_name = None
    for name, test_mol in name_to_mol.items():
        if Chem.MolToSmiles(mol, isomericSmiles=False) == Chem.MolToSmiles(test_mol, isomericSmiles=False):
            mol_name = name
            conformer_count[name] += 1
            break

    if mol_name is None:
        mol_name = f"unknown_{Chem.MolToSmiles(mol)}"

    energy_df = pd.DataFrame(columns=["Conformer Order", "MMFF94 Energy"])
    mol_h = AllChem.AddHs(mol, addCoords=True)
    AllChem.MMFFOptimizeMolecule(mol_h, maxIters=1000, ignoreInterfragInteractions=False)
    ff = AllChem.MMFFGetMoleculeForceField(mol_h, AllChem.MMFFGetMoleculeProperties(mol_h))
    energy = ff.CalcEnergy()

    new_row = pd.DataFrame({"Conformer Order": [conformer_count[mol_name]], "MMFF94 Energy": [energy]})
    energy_df = pd.concat([energy_df, new_row], ignore_index=True)

    output_filename = f"{mol_name}_{conformer_count[mol_name]}_MMFF94_ConfVAEoutput.sdf"
    writer = SDWriter(output_filename)
    writer.write(mol_h)
    writer.close()

    csv_filename = f"{mol_name}_MMFF94.csv"
    # Append the energy for the conformer to the CSV instead of overwriting
    if not os.path.exists(csv_filename):
        energy_df.to_csv(csv_filename, index=False, mode='a', header=True)
    else:
        energy_df.to_csv(csv_filename, index=False, mode='a', header=False)

print("SDF files and CSV files have been created.")
