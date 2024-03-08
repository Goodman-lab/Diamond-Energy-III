# Project Overview

This project aims to use the low-energy conformers obtained from previous research to test the effectiveness of current leading machine learning (ML)-based conformer search algorithms. The goal is to ascertain how quickly these algorithms can learn the space of diamond lattice-like low-energy conformers. Subsequently, we evaluate whether the well-trained models, informed by the learned low-energy space fitting into a diamond lattice, can guide a rapid conformer search.

## Scripts for Data Preparation

The scripts contained in this directory are used to prepare the data from Diamond Energy Projects for training the following published ML models. The preparation process involves several models, including CGCF, ConfVAE, GeoMol, ConfGF, and GeoDiff. For a more detailed introduction and instructions on how to use these ML models, refer to the links provided below.

### ML Models and Their Corresponding Data Preparation Scripts

1. **CGCF** - [Prepare_CGCFdata.py](https://github.com/MinkaiXu/CGCF-ConfGen)
2. **ConfVAE** - [Prepare_ConfVAEdata.py](https://github.com/MinkaiXu/ConfVAE-ICML21)
3. **GeoMol** - [Prepare_GeoMoldata.py](https://github.com/PattanaikL/GeoMol)
4. **ConfGF** - [Prepare_ConfGFdata.py](https://github.com/DeepGraphLearning/ConfGF)
5. **GeoDiff** - [Prepare_Geodiffdata.py](https://github.com/MinkaiXu/GeoDiff)

## Additional Information

For further study and information on these ML models, please refer to the respective links provided above.
