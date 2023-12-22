import pandas as pd


distances_data = """
e12, CGCF_closeness, 0.004095441595441592
c8a, CGCF_closeness, 0.003138888888888889
c4, CGCF_closeness, 0.0018518518518518452
e27, CGCF_closeness, 0.004256669256669259
c17l, CGCF_closeness, 0.004722996377672633
c12p, CGCF_closeness, 0.004115226337448562
c15helix, CGCF_closeness, 0.004097024585446216
c7n, CGCF_closeness, 0.003437663015127803
c8a, GeoMol_closeness, 0.0034861111111111104
c15helix, GeoMol_closeness, 0.003933625892635758
c7n, GeoMol_closeness, 0.0032524778299426187
c12p, GeoMol_closeness, 0.004115226337448562
c17l, GeoMol_closeness, 0.004352626007302265
c4, GeoMol_closeness, 0.0018518518518518452
c17l, ConfVAE_closeness, 0.00447608279742572
c15helix, ConfVAE_closeness, 0.004097024585446216
c8a, ConfVAE_closeness, 0.00348611111111111
c7n, ConfVAE_closeness, 0.0026969222743870613
e27, ConfVAE_closeness, 0.004256669256669257
c12p, ConfVAE_closeness, 0.004166666666666669
e12, ConfVAE_closeness, 0.0039886039886039846
c4, ConfVAE_closeness, 0.001388888888888884
"""

# Parse the data into a dictionary
data = {}
for line in distances_data.strip().split('\n'):
    parts = line.split(',')
    molecule, model, distance = parts[0].strip(), parts[1].strip(), float(parts[2].strip())

    if molecule not in data:
        data[molecule] = {}
    data[molecule][model] = distance

# Convert to DataFrame
df = pd.DataFrame(data).T
df.reset_index(inplace=True)
df.rename(columns={'index': 'Molecule'}, inplace=True)

# Save to CSV
df.to_csv("model_closeness_comparison.csv", index=False)
