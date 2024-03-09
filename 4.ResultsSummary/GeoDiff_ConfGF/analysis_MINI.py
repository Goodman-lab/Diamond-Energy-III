import pandas as pd
from tqdm.auto import tqdm

# Load the CSV file into a DataFrame
file_path = 'Merged_Energy_Data_MINI_OPLS4.csv'
df = pd.read_csv(file_path, encoding='ascii')

# Show the head of the DataFrame
df_head = df.head()
print(df_head)

import matplotlib.pyplot as plt
import numpy as np

# Set the style
plt.style.use('seaborn-darkgrid')

# Create a color palette
palette = plt.get_cmap('Set1')

# Plot multiple lines
plt.figure(figsize=(14, 8))
num=0
for column in df.drop('Molecule', axis=1):
    num+=1
    plt.plot(df['Molecule'], df[column], marker='', color=palette(num), linewidth=1, alpha=0.9, label=column)

# Add legend
plt.legend(loc=2, ncol=2)

# Add titles
plt.title("Energy values by Molecule", loc='left', fontsize=12, fontweight=0, color='orange')
plt.xlabel("Molecule")
plt.ylabel("Energy Value(OPLS4 kcal/mol)")

# Show the plot
plt.xticks(rotation=90)
plt.tight_layout()
plt.show()