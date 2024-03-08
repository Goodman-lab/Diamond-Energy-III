mol_names=["c4", "c5", "c5n", "c6", "c6n", "c7", "c7c", "c7n", "c8", "c8a", "c8s", "c8si", "c9", "c9n", "c9s",
           "c10n", "c12p", "c12p0", "c12p1", "c15", "c15helix", "c15a", "c15b", "c15c", "c15d", "c15e", "c15f", "c15g", "c16", "c17",
           "c17l", "c18", "c18n", "c19", "c21", "c26", "c20n", "c39n",
           "e1", "e2", "e3", "e4", "e5", "e6", "e7", "e8", "e9", "e10",
           "e11", "e12", "e13", "e14", "e15", "e16", "e17", "e18", "e19", "e20",
           "e21", "e22", "e23", "e24", "e25", "e26", "e27", "e28", "e29", "e30",
           "e31", "e32", "e33", "e34", "e35", "e36", "e37", "e38", "e39", "e40"]

def mol_ID(mol_name):
    return mol_names.index(mol_name)
print(len(mol_names))



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



import os
import pickle
import copy
from rdkit import Chem
import torch
from rdkit.Chem import HybridizationType
from torch_scatter import scatter
from rdkit.Chem.rdchem import BondType as BT
from torch_geometric.data import Data

BOND_TYPES = {t: i for i, t in enumerate(BT.names.values())}
BOND_NAMES = {i: t for i, t in enumerate(BT.names.keys())}

SDF_EXTENSION = 'H.sdf'
base_path = os.getcwd()

def read_sdf_files_from_dir(dir_path):
    mol_list = []
    for file_name in os.listdir(dir_path):
        if file_name.endswith(SDF_EXTENSION):
            sdf_path = os.path.join(dir_path, file_name)
            supplier = Chem.SDMolSupplier(sdf_path)
            print(f"Reading file: {file_name}")
            mol_list += [mol for mol in supplier if mol is not None]
    return mol_list

def rdmol_to_data(mol: Chem.Mol, mol_name, smiles=None, data_cls=dict):
    assert mol.GetNumConformers() == 1
    N = mol.GetNumAtoms()

    pos = torch.tensor(mol.GetConformer(0).GetPositions(), dtype=torch.float32)

    atomic_number = []
    aromatic = []
    sp = []
    sp2 = []
    sp3 = []
    num_hs = []
    for atom in mol.GetAtoms():
        atomic_number.append(atom.GetAtomicNum())
        aromatic.append(1 if atom.GetIsAromatic() else 0)
        hybridization = atom.GetHybridization()
        sp.append(1 if hybridization == HybridizationType.SP else 0)
        sp2.append(1 if hybridization == HybridizationType.SP2 else 0)
        sp3.append(1 if hybridization == HybridizationType.SP3 else 0)

    z = torch.tensor(atomic_number, dtype=torch.long)

    row, col, edge_type = [], [], []
    for bond in mol.GetBonds():
        start, end = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        row += [start, end]
        col += [end, start]
        edge_type += 2 * [BOND_TYPES[bond.GetBondType()]]

    edge_index = torch.tensor([row, col], dtype=torch.long)
    edge_type = torch.tensor(edge_type)

    perm = (edge_index[0] * N + edge_index[1]).argsort()
    edge_index = edge_index[:, perm]
    edge_type = edge_type[perm]

    row, col = edge_index
    hs = (z == 1).to(torch.float32)

    num_hs = scatter(hs[row], col, dim_size=N, reduce='sum').tolist()


    if smiles is None:
        #smiles = Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)
        smiles=mol_name


    boltzmannweight = torch.tensor([1], dtype=torch.float32)
    totalenergy = torch.tensor([1], dtype=torch.float32)

    data = Data(atom_type=z, pos=pos, edge_index=edge_index, edge_type=edge_type,
                rdmol=copy.deepcopy(mol), smiles=smiles,
                boltzmannweight=boltzmannweight, totalenergy=totalenergy)

    # data.nx = to_networkx(data, to_undirected=True)

    return data

def main():
    mol_dict = {name: mol_ID(name) for name in mol_names}

    train_data, val_data, test_data = [], [], []

    for mol_name, geom_id in mol_dict.items():
        dir_path = os.path.join(base_path, mol_name)

        if os.path.exists(dir_path):
            print(f"Processing directory for molecule: {mol_name}")
            mol_list = read_sdf_files_from_dir(dir_path)

            for mol in mol_list:
                data_obj = rdmol_to_data(mol, mol_name)  # Convert mol into Data format

                # Modified this part to allow a molecule to belong to multiple sets
                #if mol_name in training_set:
                #    train_data.append(data_obj)
                #if mol_name in validation_set:
                #    val_data.append(data_obj)
                if mol_name in mol_names:
                    test_data.append(data_obj)

    # Save the datasets to .pkl files
    #with open('train.pkl', 'wb') as f:
    #    pickle.dump(train_data, f)
    #with open('val.pkl', 'wb') as f:
    #    pickle.dump(val_data, f)
    with open('test.pkl', 'wb') as f:
        pickle.dump(test_data, f)

    print("Data serialized into test.pkl.")
    #("Data serialized into train.pkl, val.pkl, and test.pkl.")



if __name__ == "__main__":
    main()
