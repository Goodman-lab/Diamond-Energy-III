import numpy as np
import pandas as pd
from scipy.stats import kruskal

# Corrected definitions of distances
cgcf_distances = [0.004095441595441592, 0.003138888888888889, 0.0018518518518518452, 0.004256669256669259, 0.004722996377672633, 0.004115226337448562, 0.004097024585446216, 0.003437663015127803]
geomol_distances = [0.0034861111111111104, 0.003933625892635758, 0.0032524778299426187, 0.004115226337448562, 0.004352626007302265, 0.0018518518518518452]
confvae_distances = [0.00447608279742572, 0.004097024585446216, 0.00348611111111111, 0.0026969222743870613, 0.004256669256669257, 0.004166666666666669, 0.0039886039886039846, 0.001388888888888884]

# Perform Kruskal-Wallis test
statistic, p_value = kruskal(cgcf_distances, geomol_distances, confvae_distances)
print(f"Kruskal-Wallis test statistic: {statistic}, p-value: {p_value}")

# Example distances for each model
data = {
    'CGCF_closeness': cgcf_distances,
    'GeoMol_closeness': geomol_distances,
    'ConfVAE_closeness': confvae_distances
}

# Handle unequal lengths
max_len = max(len(cgcf_distances), len(geomol_distances), len(confvae_distances))
for key in data.keys():
    data[key] += [np.nan] * (max_len - len(data[key]))

df = pd.DataFrame(data)
df.to_csv("model_comparisons.csv", index=False)

