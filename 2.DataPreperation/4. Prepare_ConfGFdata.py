
import os
import shutil
from zipfile import ZipFile

# Initialize molecule names and IDs
mol_names = [
    "c4", "c5", "c5n", "c6", "c6n", "c7", "c7c", "c7n", "c8", "c8a", "c8s", "c8si",
    "c9", "c9n", "c9s", "c10n", "c12p", "c12p0", "c12p1", "c15", "c15helix",
    "c15a", "c15b", "c15c", "c15d", "c15e", "c15f", "c15g", "c16", "c17",
    "c17l", "c18", "c18n", "c19", "c21", "c26",
    "e1", "e2", "e3", "e4", "e5", "e6", "e7", "e8", "e9",
    "e10", "e11", "e12", "e13", "e14", "e15", "e16", "e17", "e18", "e19",
    "e20", "e21", "e22", "e23", "e24", "e25", "e26", "e27", "e28", "e29",
    "e30", "e31", "e32", "e33", "e34", "e35", "e36", "e37"
]
mol_ID = [
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14,
    15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29,
    30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40,
    41, 42, 43, 44, 45, 46, 47, 48, 49, 50,
    51, 52, 53, 54, 55, 56, 57, 58, 59, 60,
    61, 62, 63, 64, 65, 66, 67, 68, 69, 70
]

# Create a new directory to store the SDF files
target_dir = "ConfGF_diamonddata"
if not os.path.exists(target_dir):
    os.makedirs(target_dir)

# Loop through each molecule name and ID
for mol in mol_names:
    # Create a path to the current molecule subdirectory
    mol_subdir = os.path.join(os.getcwd(), mol)

    # Check if the subdirectory exists
    if os.path.exists(mol_subdir):
        # Loop through all the files in the current molecule subdirectory
        for filename in os.listdir(mol_subdir):
            # Check if the file ends with "_diamond.sdf"
            if filename.endswith("_diamond.sdf"):
                # Generate the full file path
                src_file_path = os.path.join(mol_subdir, filename)

                # Generate the target file path
                target_file_path = os.path.join(target_dir, filename)

                # Copy the file
                shutil.copy(src_file_path, target_file_path)

# Create a ZIP file
zip_filename = "ConfGF_diamonddata.zip"
with ZipFile(zip_filename, 'w') as zipf:
    for root, dirs, files in os.walk(target_dir):
        for filename in files:
            zipf.write(os.path.join(root, filename), os.path.relpath(os.path.join(root, filename), target_dir))

print(f"Copying complete. The files have been zipped into {zip_filename}")
