# README

**TLDR** Development of ML method for EACN prediction. The best shallow method: SVR ($R^2$ at cross-validation is 0.92, RMSE is 1.40). The best DL method: SGIR, ($R^2$ is 0.79, RMSE is 2.73 for test set). They are realized as IPython and python calculator.

You will need following libraries to make EACN calculation:
1. [RDKit](https://www.rdkit.org/)
2. [DataMol](https://datamol.io/) and [MolFeat](https://molfeat.datamol.io/)
3. [Padelpy](https://github.com/ecrl/padelpy)
4. [SciKit-learn](https://scikit-learn.org/stable/)
5. [Torch_molecule](https://github.com/liugangcode/torch-molecule) with CUDA installed.

Some other libraries are alco encountered in files (like [cirpy](https://cirpy.readthedocs.io/en/latest/) or [mols2grid](https://github.com/cbouy/mols2grid)), although they are not necessary for EACN calculation. 

This set of files is created to find our the best model for the prediction of EACN of different compounds, both non-polar hydrocarbons and different polar oils. Our approach inherits the one described in well done works: [paper 1](https://pubs.acs.org/doi/10.1021/acs.jpca.4c00936?goto=supporting-info) and [paper 2](https://pubs.acs.org/doi/full/10.1021/acsomega.2c04592).

This folder contains following files: 
1. *targets.csv* - EACN values vs names of molecules. The [reference](https://pubs.acs.org/doi/10.1021/acs.jpca.4c00936?goto=supporting-info) 
2. *1.Descriptors_calculation.ipynb* - describes the process of descriptors preparation and curation. This file has several outputs:
    * *X_descriptors_initial.pickle*, *X_padel_initial.pickle* - full set of standard RDKit descriptors and PaDEL descriptors [reference](https://onlinelibrary.wiley.com/doi/full/10.1002/jcc.21707)
    * *X_standard_dropped.pickle*, *X_padel_dropped.pickle* - set of RDKit and PaDEL descriptors after throwing out descriptors with high correlation
    * *fetures_to_drop_Standard.pickle*, *features_to_drop_Padel.pickle* - features to drop when we get rid of descriptors with high correlation
    * *imputer.pickle* - as there is a high chance, that PaDel descriptors are not calculated properly (NaNs can occur and so on), we make as instance of imputer to fill in the gaps in dataset, if they are present
    * *Targets_with_SMILES_curated.csv* - there were some problems during transformation of names into SMILES, therefore we created a "good" dataset where all hard-to-deal molecules are transformed into corresponding SMILES
    * *EACN.pickle* - target values in form of PandasDataframe series
3. *2. ML_standard_hp.ipynb* - describes the process of the search of optimal shale ML algorithm for the EACN prediction and its hyperparameters. This file has one output:
    * *score_hyper_opt.pickle* - a dictionary, that contains the performance of every algorithm: $R^2$, RMSE, MAE for cross-validation, train and test sets
4. *2a.Outliers.ipynb* - this file is created to visualize molecules, that are hard-to-predict
5. *2b.ML_torch_molecule.ipynb* - here we're looking for the best Deep learning ML algorithm for the search of best EACN prediction from [Torch-molecule](https://github.com/liugangcode/torch-molecule) library. It also makes a final torch model for prediction. It has several outputs:
    * *scores_torch.pickle* - dictionary, that contains performance data:  $R^2$, RMSE, MAE for train and test sets.
    * *sgir_final.pt* - trained instance of SGIR Deep learning model for EACN prediction
6. *3b.Final_calculator_SVR.ipynb* - describes the process of training shale algorithm (SVR) prediction model. It is combined with IPython realization of EACN calculator
    * *trained_pipe.pickle* is the only output, it is the instance of trained pipe (scaler + SVR model) for EACN prediction
7. *3b.Final_calculator_torch.ipynb* - IPython realization of EACN calculator with DL method SGIR.
8. *function.py* - containes *calculate_EACN* funcion, which will the be used in framework for surfactant-containing systems.
9. *calculator_py.py* is used to calculate EACN values in terminal. Just launch it with ```python calculator_py.py```  and follow the prompts

*Trash* folder contains my previous unsuccessful efforts to cluster dataset, create my own GNN, old datasets. Stuff is not well-sorted and described, so use it with high caution.
