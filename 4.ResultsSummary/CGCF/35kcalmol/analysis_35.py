
import pandas as pd
from tqdm.auto import tqdm

# Enable progress_apply for pandas
tqdm.pandas()

# Load the CSV file into a DataFrame
df = pd.read_csv('merged_conformers_data_GeoMol3kcal_OPLS4_statistic.csv', encoding='ascii')

import matplotlib.pyplot as plt
import numpy as np

# Set the style
plt.style.use('seaborn-darkgrid')

# Create a color palette
palette = plt.get_cmap('Set1')

# Prepare data
x = df['Molecule']
columns = df.columns[1:]

# Plot each column
plt.figure(figsize=(15,10))
for i, column in enumerate(columns):
    plt.plot(x, df[column], marker='o', linestyle='-', color=palette(i), linewidth=1.5, alpha=0.9, label=column)

# Add legend
plt.legend(loc=2, ncol=2)

# Add titles
plt.title('Conformer Within 3kcal/mol of the Lowest Energy of Various Methods', loc='left', fontsize=12, fontweight=0, color='orange')
plt.xlabel('Molecule')
plt.ylabel('Count')

# Show the plot
plt.xticks(rotation=90)
plt.tight_layout()
plt.show()
