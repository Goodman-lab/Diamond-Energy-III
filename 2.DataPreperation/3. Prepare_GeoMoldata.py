mol_names=["c4", "c5", "c5n", "c6", "c6n", "c7", "c7c", "c7n", "c8", "c8a", "c8s", "c8si", "c9", "c9n", "c9s",
           "c10n", "c12p", "c12p0", "c12p1", "c15", "c15helix", "c15a", "c15b", "c15c", "c15d", "c15e", "c15f", "c15g", "c16", "c17",
           "c17l", "c18", "c18n", "c19", "c21", "c26", "c20n", "c39n",
           "e1", "e2", "e3", "e4", "e5", "e6", "e7", "e8", "e9", "e10",
           "e11", "e12", "e13", "e14", "e15", "e16", "e17", "e18", "e19", "e20",
           "e21", "e22", "e23", "e24", "e25", "e26", "e27", "e28", "e29", "e30",
           "e31", "e32", "e33", "e34", "e35", "e36", "e37", "e38", "e39", "e40"]

def mol_ID(mol_name):
    return mol_names.index(mol_name)



training_set = set(["c5", "c5n", "c6", "c6n", "c7", "c8", "c8s", "c8si", "c9", "c9n", "c9s",
           "c10n", "c12p0", "c12p1", "c15", "c15a", "c15b", "c15c", "c15e", "c15f", "c15g", "c16", "c17", "c18", "c19", "c21", "c26",
           "e1", "e2", "e3", "e4", "e5", "e6", "e8", "e10",
           "e11", "e13", "e14", "e15", "e16", "e17", "e19",
           "e21", "e22", "e23", "e24", "e25", "e26", "e28", "e29", "e30",
           "e32", "e33", "e34", "e35", "e36", "e37"])
validation_set = set(["c7c", "c15d", "c18n", "e9", "e7", "e18", "e20", "e31"])
test_set = set(["c4", "c7n", "c8a", "c12p", "c15helix", "c17l", "c20n", "c39n", "e12", "e27", "e38", "e39", "e40"])


# Function to check for overlapping molecules among sets
def check_unique_sets(*sets):
    all_mols = set()
    for s in sets:
        if len(all_mols.intersection(s)) > 0:
            overlapping_molecules = all_mols.intersection(s)
            raise ValueError(f"Overlapping molecules found: {overlapping_molecules}")
        all_mols.update(s)
# Check if the sets are unique
check_unique_sets(training_set, validation_set, test_set)


# Generate the IDs based on the sets
training_set_ID = [mol_ID(mol_name) for mol_name in training_set]
validation_set_ID = [mol_ID(mol_name) for mol_name in validation_set]
test_set_ID = [mol_ID(mol_name) for mol_name in test_set]
print("training_set_ID:", training_set_ID)
print("validation_set_ID:", validation_set_ID)
print("test_set_ID:", test_set_ID)
print("training_set_length:", len(training_set_ID))
print("validation_set_length:", len(validation_set_ID))
print("test_set_length:", len(test_set_ID))



#geom_id=120834290
#totalconfs=?
#uniqueconfs=?
#conformerweights=[0.50002, 0.49998]
#rd_mol=<rdkit.Chem.rdchem.Mol object at 0x7f634784d900>
#smiles=rd_mol.rdkit_chem_to_smile("rd_mol")
#{'conformers': [{'geom_id': 120834290, 'set': 1, 'degeneracy': 2, 'totalenergy': 0.0,
#                 'relativeenergy': 0.0, 'boltzmannweight': 1.0, 'conformerweights': [0.50002, 0.49998],
#                 'rd_mol': <rdkit.Chem.rdchem.Mol object at 0x7f634784d900>}],
#                 'totalconfs': 2, 'temperature': 298.15, 'uniqueconfs': 1,
#                 'lowestenergy': 0.0 , 'poplowestpct': 100.0, 'ensembleenergy': 0.0,
#                 'ensembleentropy': 0.0, 'ensemblefreeenergy': 0.0, 'charge': 0,
#                 'smiles': 'O=[N+]([O-])c1ncco1'}


import os
import pickle
from rdkit import Chem
import numpy as np


# Save the training, validation, and test set IDs
# Create NumPy array with the desired data
split0_data = np.array([np.array(training_set_ID, dtype=np.int),
                        np.array(validation_set_ID, dtype=np.int),
                        np.array(test_set_ID, dtype=np.int)])
# Save the NumPy array into a file named 'split0.npy'
np.save('split0.npy', split0_data)
print("Saved training, validation, and test set IDs into 'split0.npy'")



SDF_EXTENSION = 'H.sdf'
base_path = os.getcwd()

def read_sdf_files_from_dir(dir_path):
    mol_list = []
    for file_name in os.listdir(dir_path):
        if file_name.endswith(SDF_EXTENSION):
            sdf_path = os.path.join(dir_path, file_name)
            supplier = Chem.SDMolSupplier(sdf_path)
            mol_list += [mol for mol in supplier if mol is not None]
    return mol_list

def merge_conformers(mol_list):
    if not mol_list:
        return None
    merged_mol = Chem.Mol(mol_list[0])
    for mol in mol_list[1:]:
        Chem.CombineMols(merged_mol, mol)
    return merged_mol

def create_molecule_dict(merged_mol, geom_id, len_mol_list):
    return {
        'geom_id': geom_id,
        'set': 1,
        'degeneracy': 2,
        'totalenergy': 0.0,
        'relativeenergy': 0.0,
        'boltzmannweight': 1.0,
        'conformerweights': [1.0] * len_mol_list,
        'rd_mol': merged_mol,
    }

def main():
    mol_dict = {name: mol_ID(name) for name in mol_names}

    for mol_name, geom_id in mol_dict.items():
        dir_path = os.path.join(base_path, mol_name)

        if os.path.exists(dir_path):
            mol_list = read_sdf_files_from_dir(dir_path)

            if mol_list:
                # Merge conformers into one molecule
                merged_mol = merge_conformers(mol_list)

                if merged_mol:
                    mol_data_dict = {'conformers': [create_molecule_dict(merged_mol, geom_id, len(mol_list))]}

                    # Update additional attributes
                    mol_data_dict.update({
                        'totalconfs': len(mol_list),
                        'uniqueconfs': len(mol_list),  # Only one merged conformer
                        'temperature': 298.15,
                        'lowestenergy': 0.0,
                        'poplowestpct': 100.0,
                        'ensembleenergy': 0.0,
                        'ensembleentropy': 0.0,
                        'ensemblefreeenergy': 0.0,
                        'charge': 0,
                        'smiles': Chem.MolToSmiles(merged_mol),
                    })

                    pickle_file_name = f"{mol_data_dict['smiles']}.pkl"

                    if not os.path.exists(pickle_file_name):
                        with open(pickle_file_name, 'wb') as f:
                            pickle.dump(mol_data_dict, f)
                        print(f"Serialized 1 conformer for {mol_name} into {pickle_file_name}")
                    else:
                        print(f"File {pickle_file_name} already exists. Skipping.")

if __name__ == "__main__":
    main()
