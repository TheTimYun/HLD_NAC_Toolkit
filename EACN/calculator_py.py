from warnings import filterwarnings
filterwarnings(action = 'ignore')
import datamol as dm
from molfeat_padel.calc import PadelDescriptors
from sklearn.svm import SVR
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
import pandas as pd
import numpy as np
import pickle
from function import calculate_EACN


gpu_available = int(input("Enter 1, if you have GPU card with CUDA installed and you want to get the results of SGIR prediction, enter 0, if you don't"))

                         
if gpu_available:
    from torch_molecule import SGIRMolecularPredictor
    sgir = SGIRMolecularPredictor()
    sgir.load_from_local('sgir_final.pt')


with open('trained_pipe.pickle', 'rb') as inp:
    pipe = pickle.load(inp)

with open('imputer.pickle', 'rb') as inp:
        imputer = pickle.load(inp) 

with open('app_domain.pickle', 'rb') as inp:
        app_domain = pickle.load(inp) 

with open('features_to_drop_Padel.pickle', 'rb') as inp:
        features_to_drop = pickle.load(inp)

smi = input('Please, enter SMILES of your molecule')
if gpu_available:
    svr_result, sgir_result, flag_AD = calculate_EACN(smi, app_domain = app_domain,  imputer = imputer, features_to_drop = features_to_drop, svr_model = pipe, use_sgir = True, sgir_model = sgir)
    if flag_AD:
        print('You're outside the applicability domain')
    print('EACN of the molecule according to SVR is {}, according to SGIR is {}'.format(svr_result , sgir_result))
else:
    svr_result, flag_AD = calculate_EACN(smi, app_domain = app_domain, imputer = imputer, features_to_drop = features_to_drop, svr_model = pipe)
    if flag_AD:
        print('You're outside the applicability domain')
    print('EACN of the molecule according to SVR is {}'.format(svr_result ))
    
    