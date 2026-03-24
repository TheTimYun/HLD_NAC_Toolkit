# README

**TLDR** Development of ML models for the prediction of critical curvature (Cc) of anionic, cationic, non-ionic and extended surfactants. For anionic surfactants it is Linear regression with PaDEL descriptors ($R^{2}$ is 0.58, RMSE for LOO cross-validation is 1.19). For non-ionics - Ridge with RDKit descriptors ($R^{2}$ is 0.68, RMSE for LOO cross-validation is 2.19). For cationics - Linear Regression with United descriptors ($R^{2}$ is 0.73, RMSE for LOO cross-validation is 1.82). For extended - Linear Regression with PaDEL descriptors ($R^{2}$ is 0.55, RMSE for LOO cross-validation is 0.71).


You will need following libraries to make Сс calculation:
1. [RDKit](https://www.rdkit.org/)
2. [DataMol](https://datamol.io/) and [MolFeat](https://molfeat.datamol.io/)
3. [Padelpy](https://github.com/ecrl/padelpy)
4. [SciKit-learn](https://scikit-learn.org/stable/)

Some other libraries are alco encountered in files (like [cirpy](https://cirpy.readthedocs.io/en/latest/) or [mols2grid](https://github.com/cbouy/mols2grid)), although they are not necessary for mail calculation. 

This set of files is created to find our the best model for the prediction of Cc

This folder contains following files : 



1. *1.Descriptors_calculation.ipynb* - describes the process of descriptors preparation and curation. This file has several outputs:
2. *2.ML_Padel_descs.ipynb* describes the process of finding optimal hyperparameters and training shallow ML algorithms for PaDEL descs. 
3. *2.ML_RDkit_descs.ipynb* describes the process of of finding optimal hyperparameters and training shallow ML algorithms for RDKit descs. 
4. *2. ML_United_descs.ipynb* describes the process of of finding optimal hyperparameters and training shallow ML algorithms for united (RDKit + PaDEL) descs. 
5. *2b.M_torch_molecule_extended.ipynb* and *2b.M_torch_molecule_non-ionicipynb* are not very useful. In these files we tried to predict Cc with deep learning models, implemented in [torch-molecule](https://github.com/liugangcode/torch-molecule). The results were not satisfactory, but we left it here. 
6. *3.Creating_models.ipynb* - create final instances of trained models with optimal hyperparameters, found in files of series *2.ML...* for the prediction of Cc. Outputs, used in calculator:
   * *models_for_prediction_anionic.pickle*
   * *models_for_prediction_cationic.pickle*
   * *models_for_prediction_extended.pickle*
   * *models_for_prediction_non_ionic.pickle*
8. *functions.py* - contains several useful functions (calculation of desciptors, embedding of the molecule) to be further used in calculator
9. *calculator_py.py* is used to calculate Cc values in terminal. Just launch it with ```python calculator_py.py```  and follow the prompts



Folder **DATA**
1. *Dataset.xlsx* and *Dataset_without_reference.xlsx* - Excel worksheets contain names of the molecule, SMILES, Cc, upper and lower boundaries and references (first file)
2. *Data_ionic.csv*, *Data_non_ionic.csv*, *Data_extended.csv*, *Data_cationic.csv* - are just parts of the dataset for the corresponding classes of surfactants. For surfactants, that have multiple instances in the initial dataset, values of Cc are averaged;

Folder **Features_PaDEL**:
1. *padel_descks.pickle* -  list of Pandas Dataframes with corresponding descriptors
2. *features_to_drop_padel.pickle* -  list of names of highly correalted features, that are dropped
3. *imputers_for_padel.pickle* - is a list with set of trained imputers for each classes of surfactants for filling the gaps in PaDEL descriptors (padelpy often do not compute descriptors properly, leaving some NaNs)

Folder **Features_RDKit**:
1. *rdkit_descs.pickle* -  list of Pandas Dataframes with corresponding descriptors
2.  *features_to_drop_rdkit.pickle* -  list of names of highly correalted features, that are dropped

Folder **Features_United**:
1.  *united_descs.pickle*-  list of Pandas Dataframes with corresponding descriptors
2.  *features_to_drop_united.pickle* -  list of names of highly correalted features, that are dropped

Folder **Scores**
1. *scores_padel_descs.pickle*, *scores_rdkit_descs.pickle*, *scores_united_descs.pickle*,  *scores_extended_torch.pickle*, *scores_non_ionic_torch.pickle* -  list of dictionaries, that contain the performance of every algorithm for every class of surfactant: $R^2$, RMSE, MAE for cross-validation, train and test sets