import re
import pandas as pd
import glob

# Define a pattern to extract the energy value from the log files
energy_pattern = re.compile(r"Total Energy =\s+([\d.]+) kcal/mol")

# Initialize a dictionary to hold the energy values and conformer counts for each molecule
energies_mm3 = {}

# Assuming the log files are in the current directory, list all log files
log_files = glob.glob('*_ConfVAEoutput_MM3_fretest.log')

# Loop through each log file
for log_file in log_files:
    # Extract the molecule name from the filename
    molecule_name = log_file.split('_')[0]

    # Read the content from the file
    with open(log_file, 'r') as file:
        content = file.read()

    # Find all energy values in the content
    energy_values = energy_pattern.findall(content)
    if energy_values:
        # Convert energy values to float
        energy_values = [float(value) for value in energy_values]

        # Store all energy values for the molecule
        if molecule_name not in energies_mm3:
            energies_mm3[molecule_name] = energy_values
        else:
            energies_mm3[molecule_name].extend(energy_values)

# Prepare data for DataFrame
data = []

# Calculate conformer counts within energy ranges
for molecule, energies in energies_mm3.items():
    min_energy = min(energies)
    count_3kcal = sum(1 for e in energies if e <= min_energy + 3)
    count_5kcal = sum(1 for e in energies if e <= min_energy + 5)
    data.append([molecule, min_energy, count_3kcal, count_5kcal])

# Convert the data to a DataFrame
df = pd.DataFrame(data, columns=['Molecule', 'Lowest OPLS4 Energy', 'Conformers within 3 kcal/mol', 'Conformers within 5 kcal/mol'])

# Save the DataFrame to a new CSV file
df.to_csv('MM3_energy_and_conformers.csv', index=False)

# Output the DataFrame to verify
print(df.head())

