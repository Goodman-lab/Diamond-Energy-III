import pandas as pd
from tqdm import tqdm
import matplotlib.pyplot as plt
import seaborn as sns

tqdm.pandas()

# Load the files
file_names = ['CGCF_Diamond_MMFF94.csv', 'CGCF_ISO17_MMFF94.csv', 'CGCF_qm9_MMFF94.csv', 'CGCF_drugs_MMFF94.csv']
dataframes = {}

# Load each file, sort by 'Molecule', and rename the 'Lowest MMFF94 Energy' column to the file name without '.csv'
for file in tqdm(file_names):
    df = pd.read_csv(file).sort_values(by='Molecule')
    # Create a new column name by removing '.csv' from the file name
    new_column_name = file.replace('.csv', '')
    # Rename the 'Lowest MMFF94 Energy' column
    df.rename(columns={'Lowest MMFF94 Energy': new_column_name}, inplace=True)
    # Add the dataframe to the dictionary with the new column name as the key
    dataframes[new_column_name] = df[['Molecule', new_column_name]]

# Assuming the first dataframe has been correctly processed and has the 'Molecule' column
first_key = next(iter(dataframes))
merged_df = dataframes[first_key]

# Merge the rest of the dataframes on 'Molecule'
for key in dataframes:
    if key != first_key:
        merged_df = merged_df.merge(dataframes[key], on='Molecule', how='outer')

# Show the head of the final dataframe
print(merged_df.head())

# Save the merged dataframe to a CSV file
#merged_df.to_csv('merged_data.csv', index=False)
#print("Merged DataFrame saved to 'merged_data.csv'.")


# Calculate the mean of the energy columns for each method
mean_energies = merged_df.iloc[:, 1:].mean()

# Create a dataframe for the mean energies
mean_energies_df = pd.DataFrame(mean_energies).reset_index()
mean_energies_df.columns = ['Method', 'Mean Energy']

# Show the dataframe
print(mean_energies_df)

# Create a bar plot for the mean energies
plt.figure(figsize=(10, 6))
sns.barplot(x='Mean Energy', y='Method', data=mean_energies_df, palette='viridis')
plt.title('Mean Lowest MMFF94 Energy for Each Method')
plt.xlabel('Mean Energy (kcal/mol)')
plt.ylabel('Method')
plt.tight_layout()  # Adjust the layout to fit everything
plt.show()

