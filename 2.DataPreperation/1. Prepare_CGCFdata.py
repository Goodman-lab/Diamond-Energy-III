mol_names=["c4", "c5", "c5n", "c6", "c6n", "c7", "c7c", "c7n", "c8", "c8a", "c8s", "c8si", "c9", "c9n", "c9s",
           "c10n", "c12p", "c12p0", "c12p1", "c15", "c15helix", "c15a", "c15b", "c15c", "c15d", "c15e", "c15f", "c15g", "c16", "c17",
           "c17l", "c18", "c18n", "c19", "c21", "c26", "c20n", "c39n",
           "e1", "e2", "e3", "e4", "e5", "e6", "e7", "e8", "e9", "e10",
           "e11", "e12", "e13", "e14", "e15", "e16", "e17", "e18", "e19", "e20",
           "e21", "e22", "e23", "e24", "e25", "e26", "e27", "e28", "e29", "e30",
           "e31", "e32", "e33", "e34", "e35", "e36", "e37", "e38", "e39", "e40"]

def mol_ID(mol_name):
    return mol_names.index(mol_name)

print("len_mol_names:", len(mol_names))


training_set = set(["c5", "c5n", "c6", "c6n", "c7", "c8", "c8s", "c8si", "c9", "c9n", "c9s",
           "c10n", "c12p0", "c12p1", "c15", "c15a", "c15b", "c15c", "c15e", "c15f", "c15g", "c16", "c17", "c18", "c19", "c21", "c26",
           "e1", "e2", "e3", "e4", "e5", "e6", "e8", "e10",
           "e11", "e13", "e14", "e15", "e16", "e17", "e19",
           "e21", "e22", "e23", "e24", "e25", "e26", "e28", "e29", "e30",
           "e32", "e33", "e34", "e35", "e36", "e37"])
validation_set = set(["c7c", "c15d", "c18n", "e9", "e7", "e18", "e20", "e31"])
test_set = set(["c4", "c7n", "c8a", "c12p", "c15helix", "c17l", "c20n", "c39n", "e12", "e27", "e38", "e39", "e40"])


print("len_training_set:", len(training_set))
print("len_validation_set:", len(validation_set))
print("len_test_set:", len(test_set))



import os
import pickle
import pandas as pd
from rdkit import Chem


def categorize_molecule(subdir_name, mol, training, validation, test):
    if subdir_name in training_set:
        training.append(mol)
    if subdir_name in validation_set:
        validation.append(mol)
    if subdir_name in test_set:
        test.append(mol)


def update_dataframe(subdir_name, count, data_frame):
    new_idx = len(data_frame)
    if subdir_name in training_set:
        data_frame.loc[new_idx] = [subdir_name, count]
        new_idx += 1
    if subdir_name in validation_set:
        data_frame.loc[new_idx] = [subdir_name, count]
        new_idx += 1
    if subdir_name in test_set:
        data_frame.loc[new_idx] = [subdir_name, count]




training_data = []
validation_data = []
test_data = []


df_data = pd.DataFrame(columns=['Molecule', 'Conformer_Count'])

for subdir, _, files in os.walk('.'):
    subdir_name = os.path.basename(subdir)
    current_mol_object_count = 0  # Reset the counter for each subdirectory
    if subdir_name in training_set or subdir_name in validation_set or subdir_name in test_set:
        for file in files:
            if file.endswith("H.sdf"):
                file_path = os.path.join(subdir, file)
                suppl = Chem.SDMolSupplier(file_path)
                for mol in suppl:
                    if mol is not None:
                        categorize_molecule(subdir_name, mol, training_data, validation_data, test_data)
                        current_mol_object_count += 1  # Update the counter for the molecules in this subdirectory

        # Update the DataFrame only once per subdirectory
        if current_mol_object_count > 0:
            new_idx = len(df_data)
            df_data.loc[new_idx] = [subdir_name, current_mol_object_count]


df_data.reset_index(drop=True, inplace=True)

df_train = df_data[df_data['Molecule'].isin(training_set)]
df_val = df_data[df_data['Molecule'].isin(validation_set)]
df_test = df_data[df_data['Molecule'].isin(test_set)]

df_train.to_csv('training_data.csv', index=False)
df_val.to_csv('validation_data.csv', index=False)
df_test.to_csv('test_set_data.csv', index=False)

with open("train.pkl", "wb") as f:
    pickle.dump(training_data, f)
with open("val.pkl", "wb") as f:
    pickle.dump(validation_data, f)
with open("test.pkl", "wb") as f:
    pickle.dump(test_data, f)
