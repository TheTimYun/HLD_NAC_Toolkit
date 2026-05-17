from warnings import filterwarnings
filterwarnings(action = 'ignore')
from functions import embed_optimize, calculate_cc, calculate_padel_descriptors, calculate_rdkit_descriptors, calculate_leverage
import pickle
import numpy as np
import pandas as pd


#loadign models
with open('models_for_prediction_anionic.pickle', 'rb') as inp:
    models = pickle.load(inp)
    anionic_model_svr = models[0]
    anionic_model_ridge = models[1]

with open('models_for_prediction_extended.pickle', 'rb') as inp:
    extended_model = pickle.load(inp)[0]

with open('models_for_prediction_non_ionic.pickle', 'rb') as inp:
    models = pickle.load(inp)
    non_ionic_model_svr = models[0]
    non_ionic_model_ridge = models[1]

with open('models_for_prediction_cationic.pickle', 'rb') as inp:
    cationic_model = pickle.load(inp)[0]

#Loading features to drop
with open('Features_PaDEL/featured_to_drop_padel.pickle', 'rb') as inp:
    features_to_drop_padel = pickle.load(inp)
    features_to_drop_anionic = features_to_drop_padel['PaDEL_anionic']
    features_to_drop_cationic_padel = features_to_drop_padel['PaDEL_cationic']
    features_to_drop_extended = features_to_drop_padel['PaDEL_extended']

with open('Features_RDKit/features_to_drop_rdkit.pickle', 'rb') as inp:
    features_to_drop_nonionic = pickle.load(inp)['RDKIt_non_ionic']
    
with open('Features_RDKit/features_to_drop_rdkit.pickle', 'rb') as inp:    
    features_to_drop_cationic_rdkit = pickle.load(inp)['RDKit_cationic']

with open('Features_United/features_to_drop_united.pickle', 'rb') as inp:
    features_to_drop_cationic = pickle.load(inp)['United_cationic']


#Loading imputers

with open('Features_PaDEL/imputers_for_padel.pickle', 'rb') as inp:
    imputers  = pickle.load(inp)
    imputer_extended = imputers['PaDEL_extended']
    imputer_anionic = imputers['PaDEL_anionic']
    imputer_cationic = imputers['PaDEL_cationic']

smi = input('Enter SMILES of you molecule ')

choice = input('Enter 0 for anionic surfactants, 1 for non-ionic, enter 2 for extended surfactant, 3 for cationic surfactant ')

if choice ==  '0':
    with open('Appl_domains_data/Anionic.pickle', 'rb') as inp:
        AData = pickle.load(inp)
        warning = AData['warning']
        invxtx = AData['invxtx']   
        pipe = AData['pipe']  
    feats = calculate_padel_descriptors(smi = smi, imputer = imputer_anionic, features_to_drop = features_to_drop_anionic, extended = False)
    leverage = calculate_leverage(feats, invxtx, pipe)
    if leverage > warning:
        print("You're outside the applicability domain! Leverage is {}, whereas warning level is {}".format(leverage, warning))
    print('Cc acc. to SVR is {}, Cc acc. to Ridge is {} '.format(calculate_cc(anionic_model_svr, feats), calculate_cc(anionic_model_ridge, feats) ))

if choice ==  '3':
    with open('Appl_domains_data/Cationic.pickle', 'rb') as inp:
        AData = pickle.load(inp)
        warning = AData['warning']
        invxtx = AData['invxtx']
        pipe = AData['pipe']  
    array_of_cationics = []
    #feats = calculate_padel_descriptors(smi = smi, imputer = imputer_cationic, features_to_drop = features_to_drop_cationic, extended = False)
    for i in range(5):
        feats_padel = calculate_padel_descriptors(smi = smi, imputer = imputer_cationic, features_to_drop = features_to_drop_cationic_padel, extended = False)
        feats_rdkit = calculate_rdkit_descriptors(smi = smi,  features_to_drop = features_to_drop_cationic_rdkit, if_non_ionic = False)
        feats_cationic = pd.concat([feats_padel, feats_rdkit], axis = 1)
        feats_cationic.drop(columns = features_to_drop_cationic, inplace=True)
        array_of_cationics.append(feats_cationic)
    
    average_descriptors = sum(array_of_cationics)/len(array_of_cationics)
    leverage = calculate_leverage(average_descriptors, invxtx, pipe)
    if leverage > warning:
        print("You're outside the applicability domain! Leverage is {}, whereas warning level is {}".format(leverage, warning))
    print('Cc acc. to Ridge is {} '.format(calculate_cc(cationic_model, feats_cationic),  ))

if choice ==  '2':
    with open('Appl_domains_data/Extended.pickle', 'rb') as inp:
        AData = pickle.load(inp)
        warning = AData['warning']
        invxtx = AData['invxtx']
        pipe = AData['pipe']  
    feats = calculate_padel_descriptors(smi = smi, imputer = imputer_extended, features_to_drop = features_to_drop_extended, extended = True)
    leverage = calculate_leverage(feats, invxtx, pipe)
    if leverage > warning:
        print("You're outside the applicability domain! Leverage is {}, whereas warning level is {}".format(leverage, warning))
    print('Cc is ', calculate_cc(extended_model, feats))

if choice == '1':
    array_svr = []
    array_linear = []
    with open('Appl_domains_data/Nonionic.pickle', 'rb') as inp:
        AData = pickle.load(inp)
        warning = AData['warning']
        invxtx = AData['invxtx']
        pipe = AData['pipe']  
    array_of_nonionics = []
    for i in range(5):
        feats = calculate_rdkit_descriptors(smi = smi, features_to_drop = features_to_drop_nonionic)
        array_svr.append(calculate_cc(non_ionic_model_svr, feats))
        array_linear.append(calculate_cc(non_ionic_model_ridge, feats))
        array_of_nonionics.append(feats)
    average_descriptors = sum(array_of_nonionics)/len(array_of_nonionics)
    leverage = calculate_leverage(average_descriptors, invxtx, pipe)
    if leverage > warning:
        print("You're outside the applicability domain! Leverage is {}, whereas warning level is {}".format(leverage, warning))
    print('Cc for SVR is {}, for linear regression is {} '.format(sum(array_svr)/10, sum(array_linear)/10))
    print('Standard deviation for SVR is {}, for linear regression is {} '.format(np.std(array_svr), np.std(array_linear)))

    