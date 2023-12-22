from rdkit import Chem
from rdkit.Chem import AllChem
import os
import numpy as np
from scipy.stats import wasserstein_distance



def read_sdf_to_mol(file_path):
    supplier = Chem.SDMolSupplier(file_path)
    molecules = []
    for mol in supplier:
        if mol is not None:
            # Add hydrogens
            mol_with_h = Chem.AddHs(mol, addCoords=True)
            # Generate 3D coordinates for the hydrogens, aligning them with the existing coordinates
            #AllChem.EmbedMolecule(mol_with_h, AllChem.ETKDG())
            #AllChem.MMFFOptimizeMolecule(mol_with_h)
            molecules.append(mol_with_h)
        else:
            print(f"Warning: Invalid molecule in {file_path}")
    return molecules


def find_dihedral_indices(molecule):
    dihedrals = []
    for atom in molecule.GetAtoms():
        if len(atom.GetNeighbors()) < 2:
            continue  # Skip atoms with less than two neighbors
        for neighbor1 in atom.GetNeighbors():
            for neighbor2 in neighbor1.GetNeighbors():
                if neighbor2.GetIdx() == atom.GetIdx():
                    continue  # Skip the back-reference to the original atom
                for neighbor3 in neighbor2.GetNeighbors():
                    if neighbor3.GetIdx() == neighbor1.GetIdx():
                        continue  # Skip the back-reference to neighbor1
                    dihedrals.append((atom.GetIdx(), neighbor1.GetIdx(), neighbor2.GetIdx(), neighbor3.GetIdx()))
    return dihedrals


def calculate_dihedral_angles(molecule):
    dihedrals = []
    dihedral_indices = find_dihedral_indices(molecule)
    for indices in dihedral_indices:
        angle = AllChem.GetDihedralDeg(molecule.GetConformer(), *indices)
        dihedrals.append(angle)
    if not dihedrals:
        print(f"Warning: No dihedral angles found for a molecule.")
    return dihedrals


# Function to create a distribution of dihedral angles
def create_angle_distribution(angles, num_bins=360):
    hist, bin_edges = np.histogram(angles, bins=num_bins, range=(-180, 180), density=True)
    return hist, bin_edges

# Function to compare angle distributions
def compare_angle_distributions(model_distribution, benchmark_distribution):
    distance = wasserstein_distance(model_distribution, benchmark_distribution)
    return distance

# Function to calculate dihedral descriptors for a list of molecules
def calculate_dihedral_descriptors(molecules):
    all_dihedral_distributions = []
    for molecule in molecules:
        if molecule is None:
            continue  # Skip invalid molecules
        dihedral_angles = calculate_dihedral_angles(molecule)
        if not dihedral_angles:  # Check for empty list
            continue
        angle_distribution, _ = create_angle_distribution(dihedral_angles)
        if np.isnan(angle_distribution).any() or np.isinf(angle_distribution).any():
            continue  # Skip distributions with NaN or infinity
        if np.sum(angle_distribution) == 0:
            continue  # Skip if the histogram is all zeros
        all_dihedral_distributions.append(angle_distribution)
    if not all_dihedral_distributions:
        print("Warning: No valid dihedral distributions for all molecules.")
        return None
    return np.array(all_dihedral_distributions)




# Function to recursively find all SDF files within a directory and its subdirectories
def find_sdf_files(directory, sdf_files):
    for entry in os.listdir(directory):
        path = os.path.join(directory, entry)
        if os.path.isdir(path):
            # Recursive call for subdirectories
            find_sdf_files(path, sdf_files)
        elif entry.endswith('.sdf') and not entry.endswith('H.sdf'):  # Change here to read only *.sdf files
            # Extract molecule name from file name, assuming it's the first part before '_'
            molecule_name = entry.split('_')[0]
            if molecule_name not in sdf_files:
                sdf_files[molecule_name] = []
            sdf_files[molecule_name].append(path)



# Main function to compute and compare dihedral distributions
def compute_and_compare(models_dirs, benchmark_dir):
    # Find all SDF files in the benchmark directory, including subdirectories
    benchmark_sdf_files = {}
    find_sdf_files(benchmark_dir, benchmark_sdf_files)

    # Read and store benchmark molecules and calculate their dihedral distributions
    benchmark_distributions = {}
    for molecule_name, file_paths in benchmark_sdf_files.items():
        benchmark_molecules = []
        for file_path in file_paths:
            benchmark_molecules.extend(read_sdf_to_mol(file_path))
        descriptors = calculate_dihedral_descriptors(benchmark_molecules)
        if descriptors is not None and descriptors.size > 0:
            benchmark_distributions[molecule_name] = np.mean(descriptors, axis=0)
        else:
            print(f"Warning: No valid dihedral descriptors for {molecule_name} in benchmark.")

    # Go through each model directory and compare with benchmark
    for model_dir in models_dirs:
        model_distributions = {}
        for file_name in os.listdir(model_dir):
            if file_name.endswith('.sdf'):
                molecule_name = file_name.split('_')[0]
                if molecule_name not in benchmark_distributions:
                    continue
                model_molecules = read_sdf_to_mol(os.path.join(model_dir, file_name))
                descriptors = calculate_dihedral_descriptors(model_molecules)
                if descriptors is not None and descriptors.size > 0:
                    model_distributions[molecule_name] = np.mean(descriptors, axis=0)
                else:
                    print(f"Warning: No valid dihedral descriptors for {molecule_name} in {model_dir}.")

        # Compare distributions for common molecules
        for molecule_name, model_distribution in model_distributions.items():
            if molecule_name in benchmark_distributions and model_distribution is not None:
                benchmark_distribution = benchmark_distributions[molecule_name]
                if benchmark_distribution is not None:
                    distance = compare_angle_distributions(model_distribution, benchmark_distribution)
                    print(f"Distribution comparison for {molecule_name} in {model_dir}: {distance}")
                else:
                    print(f"Warning: Missing valid benchmark distribution for {molecule_name}.")





# List of model directories
models_dirs = ['CGCF_closeness', 'GeoMol_closeness', 'ConfVAE_closeness']
benchmark_dir = 'Diamond_MM3_closeness'

# Compute and compare
compute_and_compare(models_dirs, benchmark_dir)
