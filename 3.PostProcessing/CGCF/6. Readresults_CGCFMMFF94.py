import os
import pandas as pd
import glob

# Assuming the CSV files are in the current directory, we'll search for all files matching the pattern
csv_files = glob.glob('*_MMFF94.csv')

# Initialize a dictionary to hold the lowest energy values for each molecule
lowest_energies = {}

# Loop through each CSV file
for file in csv_files:
    # Read the CSV file into a DataFrame
    df = pd.read_csv(file)

    # Extract the molecule name from the filename (remove the '_MMFF94.csv' part)
    molecule_name = file.replace('_MMFF94.csv', '')

    # Find the lowest energy value in the 'MMFF94 Energy' column
    lowest_energy = df['MMFF94 Energy'].min()

    # Store the lowest energy value in the dictionary using the molecule name as the key
    lowest_energies[molecule_name] = lowest_energy

# Convert the dictionary to a DataFrame
lowest_energies_df = pd.DataFrame(list(lowest_energies.items()), columns=['Molecule', 'Lowest MMFF94 Energy'])

# Save the DataFrame to a new CSV file
lowest_energies_df.to_csv('lowest_energy.csv', index=False)

# Output the first few rows of the DataFrame to verify
lowest_energies_df.head()

