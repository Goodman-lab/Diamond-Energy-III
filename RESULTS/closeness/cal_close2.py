import pandas as pd
import numpy as np

# Example distances for each model
model_distances = {
    'CGCF_closeness': [0.004095441595441592, 0.003138888888888889, 0.0018518518518518452, 0.004256669256669259, 0.004722996377672633, 0.004115226337448562, 0.004097024585446216, 0.003437663015127803],
    'GeoMol_closeness': [0.0034861111111111104, 0.003933625892635758, 0.0032524778299426187, 0.004115226337448562, 0.004352626007302265, 0.0018518518518518452],
    'ConfVAE_closeness': [0.00447608279742572, 0.004097024585446216, 0.00348611111111111, 0.0026969222743870613, 0.004256669256669257, 0.004166666666666669, 0.0039886039886039846, 0.001388888888888884]
}

# Prepare data for CSV
data_for_csv = {
    "Statistic": ["Average Distance", "Median Distance", "Standard Deviation", "Max Distance", "Min Distance"]
}

for model, distances in model_distances.items():
    distances_array = np.array(distances)
    data_for_csv[model] = [
        np.mean(distances_array),
        np.median(distances_array),
        np.std(distances_array),
        np.max(distances_array),
        np.min(distances_array)
    ]

# Convert to DataFrame
df = pd.DataFrame(data_for_csv)

# Save to CSV
df.to_csv("model_statistical_comparison.csv", index=False)

print(df)

