import pandas as pd
import glob

# List of CSV files to read from current working directory
csv_files = glob.glob("*_MMFF94.csv")

# Initialize an empty DataFrame to hold the final results
final_data = pd.DataFrame(columns=['Molecular Name', 'Lowest Energy (MMFF94-kcal/mol)',
                                   'Conformers within 3 kcal/mol', 'Conformers within 5 kcal/mol'])

# Loop through each CSV file to gather data
for csv_file in csv_files:
    # Read the CSV into a DataFrame
    df = pd.read_csv(csv_file)

    # Get the lowest energy value
    min_energy = df['MMFF94 Energy'].min()

    # Count how many conformers are within 3 kcal/mol and 5 kcal/mol of the lowest energy value
    within_3_kcal = len(df[df['MMFF94 Energy'] <= min_energy + 3])
    within_5_kcal = len(df[df['MMFF94 Energy'] <= min_energy + 5])

    # Extract the molecule name from the CSV filename
    mol_name = csv_file.split('_MMFF94.csv')[0]

    # Create a new DataFrame row
    new_row = pd.DataFrame([{
        'Molecular Name': mol_name,
        'Lowest Energy (MMFF94-kcal/mol)': min_energy,
        'Conformers within 3 kcal/mol': within_3_kcal,
        'Conformers within 5 kcal/mol': within_5_kcal
    }])

    # Concatenate the new row to the existing DataFrame
    final_data = pd.concat([final_data, new_row], ignore_index=True)

# Save the final DataFrame to a new CSV file
final_data.to_csv('ConfGF_final_results.csv', index=False)
