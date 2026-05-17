import datamol as dm
from molfeat_padel.calc import PadelDescriptors
from sklearn.svm import SVR
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
import pandas as pd
import numpy as np
import pickle
#from torch_molecule import SGIRMolecularPredictor
from rdkit import Chem
from rdkit.Chem import AllChem
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.pipeline import make_pipeline



def embed_optimize(smi):
    try:
        m = Chem.MolFromSmiles(smi)
        m = Chem.AddHs(m)
        params = AllChem.ETKDGv3()
        params.maxIterations = 1000
        params.useRandomCoords = True
        AllChem.EmbedMolecule(m, params)
        AllChem.MMFFOptimizeMolecule(m, maxIters=500)
        return m
    except:
        print(smi)


def mah(array, mean, pseudoinv):
    diff = array - mean
    return np.sqrt(diff@pseudoinv@diff.T)


    


def calculate_EACN(smi,app_domain, imputer = None, features_to_drop = None, svr_model = None, use_sgir = False, sgir_model = None):
    pipe = svr_model
    desc_calc = PadelDescriptors()
    mol = embed_optimize(smi)
    with dm.without_rdkit_log():
        descs = pd.DataFrame(desc_calc(mol), index = desc_calc.columns).transpose()
    
    descs = descs.drop(columns = features_to_drop)
    invxtx = app_domain['invxtx']
    warning = app_domain['warning']
    pipe_appdomain = app_domain['pipe']
    descs_transformed = pipe_appdomain.transform(descs)

    leverage = np.dot(np.dot(descs_transformed, invxtx), np.transpose(descs_transformed))[0][0]
    if leverage>warning:
        flag_AD = True
    else:
        flag_AD = False
    if use_sgir:
        return (pipe.predict(descs)[0], sgir_model.predict([smi])['prediction'][0][0], flag_AD)
    else:
        return (pipe.predict(descs)[0], flag_AD)


    


