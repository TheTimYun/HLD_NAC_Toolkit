from rdkit import Chem
from rdkit.Chem import AllChem
from padelpy import from_smiles
import pandas as pd
import numpy as np
from rdkit import Chem
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import Ridge
from sklearn.impute import KNNImputer
import molfeat
from molfeat.calc.descriptors import RDKitDescriptors2D, RDKitDescriptors3D
from molfeat.trans import MoleculeTransformer
import datamol as dm

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


def calculate_cc(trained_estimator, feats):
    return trained_estimator.predict(feats)[0]


def calculate_padel_descriptors(smi, imputer, features_to_drop, extended = True):
    feats = from_smiles(smi)
    feats = pd.DataFrame(feats, index = [0])
    for column in feats.columns:
        feats[column] = pd.to_numeric(feats[column], errors = 'coerce')
    if extended:
        feats = pd.DataFrame(imputer.transform(feats), columns = feats.columns)
        smarts_PO = 'OC(C)CO'
        smarts_EO = 'O[CH2][CH2]O'
        PO_frag = Chem.MolFromSmarts(smarts_PO)
        EO_frag = Chem.MolFromSmarts(smarts_EO)
        feats['EOs'] = len(Chem.MolFromSmiles(smi).GetSubstructMatches(EO_frag))
        feats['POs'] = len(Chem.MolFromSmiles(smi).GetSubstructMatches(PO_frag))
    if not extended:
        feats = pd.DataFrame(imputer.transform(feats), columns = feats.columns)
    feats.drop(columns = features_to_drop, inplace=True)
    return feats


def calculate_rdkit_descriptors(smi, features_to_drop, if_non_ionic = True):
    mol = Chem.MolFromSmiles(smi)
    optimized_molecule = embed_optimize(smi)
    calc_2d = RDKitDescriptors2D(ignore_descrs=['AvgIpc', 'Ipc'], do_not_standardize=True)
    calc_3d = RDKitDescriptors3D()
    with dm.without_rdkit_log():
        trans_2D = MoleculeTransformer(calc_2d, verbose = True)
        trans_3D = MoleculeTransformer(calc_3d, verbose = True)
        feats_2d = trans_2D(smi)
        feats_3d = trans_3D(optimized_molecule)
    feats = pd.DataFrame(np.concatenate([feats_2d, feats_3d], axis = 1), columns = calc_2d.columns + calc_3d.columns).drop(columns = 'Alerts')
    if if_non_ionic: 
        smarts_PO = 'OC(C)CO'
        smarts_EO = 'O[CH2][CH2]O'
        PO_frag = Chem.MolFromSmarts(smarts_PO)
        EO_frag = Chem.MolFromSmarts(smarts_EO)
        feats['EOs'] = len(mol.GetSubstructMatches(EO_frag))
        feats['POs'] = len(mol.GetSubstructMatches(PO_frag))
    feats.drop(columns = features_to_drop, inplace=True)
    return feats

def calculate_leverage(row, invxtx):
    return np.dot(np.dot(row, invxtx), np.transpose(row))[0][0]
