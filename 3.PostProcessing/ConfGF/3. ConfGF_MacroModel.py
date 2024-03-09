import os
import pandas as pd
import rdkit
from rdkit import Chem
from rdkit.Chem import AddHs
from rdkit.Chem import rdmolfiles

def count_files_in_directory(directory, extension):
    return len([f for f in os.listdir(directory) if f.endswith(extension)])

def run_command(command):
    os.system(command)


def generate_minifile_for_diamond(file_name_without_ext):
    content = f"""{file_name_without_ext}.mae
{file_name_without_ext}_MM3-out.mae
 MMOD       0      1      0      0     0.0000     0.0000     0.0000     0.0000
 DEBG    1003      0      0      0     0.0000     0.0000     0.0000     0.0000
 FFLD       2      1      0      0     1.0000     0.0000     0.0000     0.0000
 BDCO       0      0      0      0    41.5692 99999.0000     0.0000     0.0000
 BGIN       0      0      0      0     0.0000     0.0000     0.0000     0.0000
 READ      -1      0      0      0     0.0000     0.0000     0.0000     0.0000
 ELST      -1      0      0      0     0.0000     0.0000     0.0000     0.0000
 WRIT       0      0      0      0     0.0000     0.0000     0.0000     0.0000
 END       0      0      0      0     0.0000     0.0000     0.0000     0.0000
"""

    with open(f"{file_name_without_ext}_MM3.com", "w+") as f:
        f.write(content)


def frequency_test(file_name_without_ext):
    with open(f"{file_name_without_ext}_MM3_fretest.com", "w+") as f:
        f.write(f"{file_name_without_ext}_MM3-out.mae\n")
        f.write(f"{file_name_without_ext}_MM3_fretest-out.mae\n")
        f.write("\n")
        f.write(" DEBG     211      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
        f.write(" FFLD       2      1      0      0     1.0000     0.0000     0.0000     0.0000\n")
        f.write(" BDCO       0      0      0      0    41.5692 99999.0000     0.0000     0.0000\n")
        f.write(" READ       0      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
        f.write(" MINI       4      0    250      0     0.0000     0.0000     0.1000     0.0000\n")
        f.write(" ELST       0      1      0      0     0.0000     0.0000     0.0000     0.0000\n")
        f.write(" MTST       0      0      0      0     2.0000     0.0000     0.0000     0.0000\n")
        f.write(" RRHO       0      0      1      0   300.0000     1.0000     2.0000     0.0000\n")
    f.close()


def generate_OPLS4minifile_for_diamond(file_name_without_ext):
    content = f"""{file_name_without_ext}.mae
{file_name_without_ext}_OPLS4-out.mae
 MMOD       0      1      0      0     0.0000     0.0000     0.0000     0.0000
 DEBG    1003      0      0      0     0.0000     0.0000     0.0000     0.0000
 FFLD      16      1      0      0     1.0000     0.0000     0.0000     0.0000
 BDCO       0      0      0      0    41.5692 99999.0000     0.0000     0.0000
 BGIN       0      0      0      0     0.0000     0.0000     0.0000     0.0000
 READ      -1      0      0      0     0.0000     0.0000     0.0000     0.0000
 ELST      -1      0      0      0     0.0000     0.0000     0.0000     0.0000
 WRIT       0      0      0      0     0.0000     0.0000     0.0000     0.0000
 END       0      0      0      0     0.0000     0.0000     0.0000     0.0000
"""

    with open(f"{file_name_without_ext}_OPLS4.com", "w+") as f:
        f.write(content)

def OPLS4frequency_test(file_name_without_ext):
    with open(f"{file_name_without_ext}_OPLS4_fretest.com", "w+") as f:
        f.write(f"{file_name_without_ext}_OPLS4-out.mae\n")
        f.write(f"{file_name_without_ext}_OPLS4_fretest-out.mae\n")
        f.write("\n")
        f.write(" DEBG     211      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
        f.write(" FFLD      16      1      0      0     1.0000     0.0000     0.0000     0.0000\n")
        f.write(" BDCO       0      0      0      0    41.5692 99999.0000     0.0000     0.0000\n")
        f.write(" READ       0      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
        f.write(" MINI       4      0    250      0     0.0000     0.0000     0.1000     0.0000\n")
        f.write(" ELST       0      1      0      0     0.0000     0.0000     0.0000     0.0000\n")
        f.write(" MTST       0      0      0      0     2.0000     0.0000     0.0000     0.0000\n")
        f.write(" RRHO       0      0      1      0   300.0000     1.0000     2.0000     0.0000\n")
    f.close()


def list_files_in_directory(directory, extension, exclude_substring=None):
    if exclude_substring:
        return [f for f in os.listdir(directory) if f.endswith(extension) and exclude_substring not in f]
    else:
        return [f for f in os.listdir(directory) if f.endswith(extension)]




def getsdf_from_file(file_path):
    try:
        data = pd.read_pickle(file_path)

        if not hasattr(data, 'rdmol'):
            return None

        rdmol = data.rdmol
        file_name = data.smiles
        mol_h = AddHs(rdmol, addCoords=True)

        # Use SDWriter to write the molecule to a file
        writer = rdmolfiles.SDWriter(file_name)
        writer.write(mol_h)
        writer.close()

        return file_name
    except Exception as e:
        print(f"Error processing {file_path}. Error: {e}")
        return None



def main():
    # List all .pkl files in the current directory
    all_files = [f for f in os.listdir('.') if os.path.isfile(f) and f.endswith('.pkl')]
    for file in all_files:
        getsdf_from_file(file)

    # Convert SDF to MAE
    run_command("bash sdftomae.sh")

    # Get current directory and list the MAE files
    current_directory = os.getcwd()
    mae_files = list_files_in_directory(current_directory, ".mae", "-out")  # Exclude files with "_out" substring

    # Generate .com files for ConfGF output
    for file_name in mae_files:
        file_name_without_ext = os.path.splitext(file_name)[0]
        generate_minifile_for_diamond(file_name_without_ext)
        generate_OPLS4minifile_for_diamond(file_name_without_ext)

    # Run batch MM3 energy calculation
    for file_name in mae_files:
        file_name_without_ext = os.path.splitext(file_name)[0]
        command = f"bmin -WAIT -SUBHOST localhost:3 {file_name_without_ext}_MM3"
        os.system(command)

    # Run batch MM3 frequency test
    for file_name in mae_files:
        file_name_without_ext = os.path.splitext(file_name)[0]
        frequency_test(file_name_without_ext)
        command = f"bmin -WAIT -SUBHOST localhost:3 {file_name_without_ext}_MM3_fretest"
        os.system(command)

    # Run batch OPLS4 energy calculation
    for file_name in mae_files:
        file_name_without_ext = os.path.splitext(file_name)[0]
        command = f"bmin -WAIT -SUBHOST localhost:3 {file_name_without_ext}_OPLS4"
        os.system(command)

    # Run batch OPLS4 frequency test
    for file_name in mae_files:
        file_name_without_ext = os.path.splitext(file_name)[0]
        OPLS4frequency_test(file_name_without_ext)
        command = f"bmin -WAIT -SUBHOST localhost:3 {file_name_without_ext}_OPLS4_fretest"
        os.system(command)


if __name__ == "__main__":
    main()


