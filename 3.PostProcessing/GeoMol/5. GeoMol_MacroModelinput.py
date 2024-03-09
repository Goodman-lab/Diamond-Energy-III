
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


import pandas as pd
from rdkit import Chem
from rdkit.Chem import MolFromInchi as InchiToMol
from rdkit.Chem import AddHs
from rdkit.Chem import AllChem


# Create a reverse mapping from InChI to name
inchi_name_mapping = {v: k for k, v in test_set_inchi.items()}

# Convert the InChI strings to SMILES
inchi_smiles_mapping = {}
for inchi, name in inchi_name_mapping.items():
    mol = InchiToMol(inchi)
    if mol:
        smiles = Chem.MolToSmiles(mol)
        inchi_smiles_mapping[smiles] = name
    else:
        print(f"Failed to convert InChI to Mol: {inchi}")

# Load data from pickle
data = pd.read_pickle('qm9_test_mols.pkl')

# Process the data
for smiles, mols in data.items():
    mol_name = inchi_smiles_mapping.get(smiles)
    if mol_name:
        conf_id = 0
        for mol in mols:
            mol_h = AddHs(mol, addCoords=True)
            # Save the molecule with hydrogens to an SDF file
            writer = Chem.SDWriter(f'{mol_name}_{conf_id}.sdf')
            writer.write(mol_h)
            writer.close()

            conf_id += 1

print("Molecule conformers have been saved as .sdf files!")

