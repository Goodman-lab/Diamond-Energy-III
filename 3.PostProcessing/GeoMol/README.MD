# Diamond Energy Projects ML Conformer Generation

## Overview
Once the ML models are trained with data from Diamond Energy Projects, they can be used to generate low-energy conformers for test molecules. This allows for a comparison of their performance against traditional conformer search methods. The conformers produced by the ML models are further minimized or have their single points calculated using MMFF94, MM3, and OPLS4 force fields.

## Outputs
The direct outputs from the ML models are stored in `.pkl` files. The scripts within this directory are utilized to prepare the output data from the ML models for additional computational tasks.

## Further Information
For additional details and instructions on these ML models, you can refer to the following resources:

- [CGCF](https://github.com/MinkaiXu/CGCF-ConfGen)
- [ConfVAE](https://github.com/MinkaiXu/ConfVAE-ICML21)
- [GeoMol](https://github.com/PattanaikL/GeoMol)
- [ConfGF](https://github.com/DeepGraphLearning/ConfGF)
- [GeoDiff](https://github.com/MinkaiXu/GeoDiff)

## Data and Scripts

1. `GeoMol_Diamond_test_mols.pkl`  
   Raw output data from the GeoMol model, trained on Diamond Energy Projects data, generating structures for test molecules.

2. `GeoMol_drugs_test_mols.pkl`  
   Raw output data from the GeoMol model, trained on a drugs dataset, with structures generated for test molecules using a pretrained model provided by the author.

3. `GeoMol_qm9_test_mols.pkl`  
   Raw output data from the GeoMol model, trained on the qm9 dataset, generating structures for test molecules.

4. `GeoMol_MMFF94.py`  
   A script that processes `.pkl` structures generated by the GeoMol model, using the RDKit MMFF94 force field for minimization, and saves the results, including the lowest energy values and conformer numbers within a 3 to 5 kcal/mol range from the lowest energy, in a `.csv` file.

5. `GeoMol_MacroModelinput.py`  
   This script reads the `.pkl` file outputs from GeoMol and prepares them in `.sdf` format for further processing with the MacroModel suite.

6. `GeoMol_MacroModel.py`  
   Script for minimizing structures generated by the GeoMol model and performing frequency analysis using the OPLS4 or MM3 force field in the MacroModel package.

7. `Sdftomae.sh`  
   Shell script for converting files from `.sdf` to `.mae` format for use in the MacroModel suite.

8. `GeoMol_slurm.sh` 
   Job submit script.

9. `Readresults_GeoMolOPLS4_MM3.py`  
   Script that reads calculation results from the OPLS4 and MM3 force fields, including the lowest energy values and conformer counts within a 3 to 5 kcal/mol range from the lowest energy, and saves them in a `.csv` file.

